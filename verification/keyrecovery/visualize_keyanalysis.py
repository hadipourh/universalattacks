#!/usr/bin/env python3
"""
Visualization tool for Orthros-PRF Key Analysis Results
Analyzes and visualizes relationships between weak keys and good pairs

Usage:
    python visualize_keyanalysis.py <results_file.txt>
    python visualize_keyanalysis.py results_5r-v1.txt --plot-all
"""

import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import argparse

def parse_results_file(filename):
    """Parse the key analysis results file"""
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract attack parameters
    param_match = re.search(r'Parameters: offset=(\d+), active=\[(\d+), (\d+)\], dy1=(0x[0-9a-f]+), dy2=(0x[0-9a-f]+)', content)
    if param_match:
        params = {
            'offset': int(param_match.group(1)),
            'active': [int(param_match.group(2)), int(param_match.group(3))],
            'dy1': param_match.group(4),
            'dy2': param_match.group(5)
        }
    else:
        params = {}
    
    # Extract weak key classes
    classes = []
    pattern = r"Good pairs: \((.*?)\), Key candidates: \{(.*?)\}"
    
    for match in re.finditer(pattern, content, flags=re.IGNORECASE | re.DOTALL):
        good_pairs_str = match.group(1)
        keys_str = match.group(2)

        pair_matches = re.findall(r"\(0x([0-9a-fA-F]{2}),\s*0x([0-9a-fA-F]{2})\)", good_pairs_str)
        key_matches = re.findall(r"left\s*=\s*0x([0-9a-fA-F]{2}),\s*right\s*=\s*0x([0-9a-fA-F]{2})", keys_str)
        if not key_matches:
            # Backward compatibility with older format
            key_matches = re.findall(r"\(0x([0-9a-fA-F]{2}),\s*0x([0-9a-fA-F]{2})\)", keys_str)

        classes.append({
            'num_good_pairs': len(pair_matches),
            'num_keys': len(key_matches),
            'good_pairs_str': good_pairs_str,
            'keys_str': keys_str,
            'good_pairs': [(p[0].lower(), p[1].lower()) for p in pair_matches],
            'keys': [(k[0].lower(), k[1].lower()) for k in key_matches]
        })
    
    # Extract summary statistics
    summary = {}
    summary_match = re.search(r'Total number of keys: (\d+)', content)
    if summary_match:
        summary['total_keys'] = int(summary_match.group(1))
    
    summary_match = re.search(r'Number of weak keys: (\d+)', content)
    if summary_match:
        summary['num_weak_keys'] = int(summary_match.group(1))
    
    summary_match = re.search(r'Number of unique indices.*?: (\d+)', content)
    if summary_match:
        summary['num_classes'] = int(summary_match.group(1))
    
    summary_match = re.search(r'Size of union of all good pairs: (\d+)', content)
    if summary_match:
        summary['union_size'] = int(summary_match.group(1))
    
    return params, classes, summary

def plot_keys_vs_good_pairs_scatter(classes, params, summary, filename_prefix):
    """Scatter plot: Number of keys vs number of good pairs per class"""
    good_pairs = [c['num_good_pairs'] for c in classes]
    num_keys = [c['num_keys'] for c in classes]
    
    plt.figure(figsize=(10, 6))
    plt.scatter(good_pairs, num_keys, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    plt.xlabel('Number of Good Pairs', fontsize=12, fontweight='bold')
    plt.ylabel('Number of Weak Keys', fontsize=12, fontweight='bold')
    plt.title(f'Weak Keys vs Good Pairs Distribution\n{filename_prefix}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add statistics annotation
    textstr = f"Classes: {summary.get('num_classes', 'N/A')}\nTotal Weak Keys: {summary.get('num_weak_keys', 'N/A')}"
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_scatter.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_scatter.pdf', bbox_inches='tight')
    print(f"✓ Saved scatter plot: {filename_prefix}_scatter.png/pdf")
    plt.close()

def plot_good_pairs_distribution(classes, params, summary, filename_prefix):
    """Histogram: Distribution of good pairs per class"""
    good_pairs = [c['num_good_pairs'] for c in classes]
    
    plt.figure(figsize=(10, 6))
    bins = range(min(good_pairs), max(good_pairs) + 2, 2) if good_pairs else 10
    plt.hist(good_pairs, bins=bins, alpha=0.7, edgecolor='black', color='steelblue')
    plt.xlabel('Number of Good Pairs per Class', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency (Number of Classes)', fontsize=12, fontweight='bold')
    plt.title(f'Distribution of Good Pairs\n{filename_prefix}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add statistics
    textstr = f"Mean: {np.mean(good_pairs):.1f}\nMedian: {np.median(good_pairs):.0f}\nStd: {np.std(good_pairs):.1f}"
    plt.text(0.98, 0.98, textstr, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_good_pairs_dist.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_good_pairs_dist.pdf', bbox_inches='tight')
    print(f"✓ Saved good pairs distribution: {filename_prefix}_good_pairs_dist.png/pdf")
    plt.close()

def plot_keys_distribution(classes, params, summary, filename_prefix):
    """Histogram: Distribution of weak keys per class - shows frequency"""
    num_keys = [c['num_keys'] for c in classes]
    
    # Count frequency of each class size
    from collections import Counter
    frequency_dict = Counter(num_keys)
    class_sizes = sorted(frequency_dict.keys())
    frequencies = [frequency_dict[size] for size in class_sizes]
    
    plt.figure(figsize=(12, 7))
    
    # Create bar plot for better readability
    bars = plt.bar(class_sizes, frequencies, alpha=0.8, edgecolor='black', linewidth=1.5, color='coral', width=0.7)
    
    plt.xlabel('Number of Weak Keys per Class', fontsize=13, fontweight='bold')
    plt.ylabel('Frequency (Number of Classes)', fontsize=13, fontweight='bold')
    plt.title(f'Frequency of Weak Keys per Class\n{filename_prefix}', fontsize=15, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    # Set x-axis to show all class sizes
    plt.xticks(class_sizes, fontsize=10)
    
    # Add value labels on top of each bar
    for i, (size, freq) in enumerate(zip(class_sizes, frequencies)):
        plt.text(size, freq, f'{freq}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Add statistics box
    total_classes = len(classes)
    textstr = f"Total Classes: {total_classes}\n"
    textstr += f"Mean: {np.mean(num_keys):.1f}\n"
    textstr += f"Median: {np.median(num_keys):.0f}\n"
    textstr += f"Mode: {max(frequency_dict, key=frequency_dict.get)}\n"
    textstr += f"Std: {np.std(num_keys):.1f}"
    
    plt.text(0.98, 0.98, textstr, transform=plt.gca().transAxes, 
             fontsize=11, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7, edgecolor='black'))
    
    # Add a table with the frequency data
    table_data = []
    for size, freq in zip(class_sizes[:10], frequencies[:10]):  # Top 10
        pct = (freq / total_classes) * 100
        table_data.append([f"{size} keys", f"{freq} classes", f"{pct:.1f}%"])
    
    if len(class_sizes) > 10:
        remaining = sum(frequencies[10:])
        pct = (remaining / total_classes) * 100
        table_data.append(["Others", f"{remaining} classes", f"{pct:.1f}%"])
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_keys_freq.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_keys_freq.pdf', bbox_inches='tight')
    print(f"✓ Saved weak keys frequency: {filename_prefix}_keys_freq.png/pdf")
    plt.close()

def plot_heatmap(classes, params, summary, filename_prefix):
    """2D histogram/heatmap: Density of (good_pairs, num_keys) combinations"""
    good_pairs = np.array([c['num_good_pairs'] for c in classes])
    num_keys = np.array([c['num_keys'] for c in classes])
    
    plt.figure(figsize=(12, 8))
    
    # Create 2D histogram
    hist, xedges, yedges = np.histogram2d(good_pairs, num_keys, bins=[20, 20])
    
    plt.imshow(hist.T, origin='lower', aspect='auto', cmap='YlOrRd', interpolation='nearest',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plt.colorbar(label='Number of Classes')
    
    plt.xlabel('Number of Good Pairs', fontsize=12, fontweight='bold')
    plt.ylabel('Number of Weak Keys', fontsize=12, fontweight='bold')
    plt.title(f'Density Heatmap: Keys vs Good Pairs\n{filename_prefix}', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_heatmap.pdf', bbox_inches='tight')
    print(f"✓ Saved heatmap: {filename_prefix}_heatmap.png/pdf")
    plt.close()

def plot_top_classes_bar(classes, params, summary, filename_prefix, top_n=20):
    """Bar chart: Top N classes by number of weak keys"""
    # Sort by number of keys
    sorted_classes = sorted(classes, key=lambda x: x['num_keys'], reverse=True)[:top_n]
    
    labels = [f"{c['num_good_pairs']} pairs" for c in sorted_classes]
    keys = [c['num_keys'] for c in sorted_classes]
    
    plt.figure(figsize=(12, 8))
    bars = plt.barh(range(len(labels)), keys, color='mediumseagreen', edgecolor='black', linewidth=0.5)
    plt.yticks(range(len(labels)), labels, fontsize=9)
    plt.xlabel('Number of Weak Keys', fontsize=12, fontweight='bold')
    plt.ylabel('Class (by # Good Pairs)', fontsize=12, fontweight='bold')
    plt.title(f'Top {top_n} Classes by Weak Key Count\n{filename_prefix}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='x')
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        width = bar.get_width()
        plt.text(width, bar.get_y() + bar.get_height()/2, f' {int(width)}',
                ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_top_classes.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_top_classes.pdf', bbox_inches='tight')
    print(f"✓ Saved top classes bar chart: {filename_prefix}_top_classes.png/pdf")
    plt.close()

def plot_frequency_table(classes, params, summary, filename_prefix, hide_config=False):
    """Create frequency distribution showing number of weak keys per class"""
    from collections import Counter
    num_keys = [c['num_keys'] for c in classes]
    frequency_dict = Counter(num_keys)
    
    # Sort by class size
    class_sizes = sorted(frequency_dict.keys())
    frequencies = [frequency_dict[size] for size in class_sizes]
    total_classes = len(classes)
    percentages = [(f/total_classes)*100 for f in frequencies]
    
    # Calculate statistics
    mean = sum(num_keys) / len(num_keys)
    median = sorted(num_keys)[len(num_keys)//2]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Bar chart with color matching paper's histogram style
    bars = ax.bar(class_sizes, frequencies, color='#4B6CD6', alpha=0.85, 
                  edgecolor='white', linewidth=1.0, width=0.7)
    
    # Add value labels on bars
    for size, freq, pct in zip(class_sizes, frequencies, percentages):
        ax.text(size, freq + max(frequencies)*0.01,
               f'{freq}\n({pct:.1f}%)',
               ha='center', va='bottom', fontsize=10)
    
    ax.set_xlabel('Number of weak keys in weak-key class $\\mathcal{W}_i$', fontsize=12)
    ax.set_ylabel('Number of weak-key classes', fontsize=12)
    ax.set_xticks(class_sizes)
    ax.set_ylim(bottom=0, top=max(frequencies)*1.18)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.7)
    ax.set_axisbelow(True)
    
    # Add legend with parameters and statistics
    if hide_config:
        legend_text = f"$\\Delta_1={params['dy1']}$, $\\Delta_2={params['dy2']}$\n"
    else:
        legend_text = f"Offset={params['offset']}, Active S-boxes={params['active']}\n"
        legend_text += f"$\\Delta_1={params['dy1']}$, $\\Delta_2={params['dy2']}$\n"
    legend_text += f"$\\mu = {mean:.2f}$, median $= {median}$\n"
    legend_text += f"Total: {total_classes} classes, {sum(num_keys)} weak keys"
    
    ax.text(0.98, 0.97, legend_text,
           transform=ax.transAxes,
           fontsize=11,
           verticalalignment='top',
           horizontalalignment='right',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                    edgecolor='gray', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_frequency_table.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_frequency_table.pdf', bbox_inches='tight')
    print(f"✓ Saved frequency distribution: {filename_prefix}_frequency_table.png/pdf")
    plt.close()

def plot_combined_dashboard(classes, params, summary, filename_prefix):
    """Combined dashboard with multiple subplots"""
    good_pairs = [c['num_good_pairs'] for c in classes]
    num_keys = [c['num_keys'] for c in classes]
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # 1. Scatter plot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(good_pairs, num_keys, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    ax1.set_xlabel('Good Pairs', fontweight='bold')
    ax1.set_ylabel('Weak Keys', fontweight='bold')
    ax1.set_title('Keys vs Good Pairs', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Good pairs distribution
    ax2 = fig.add_subplot(gs[0, 1])
    bins_gp = range(min(good_pairs), max(good_pairs) + 2, 2) if good_pairs else 10
    ax2.hist(good_pairs, bins=bins_gp, alpha=0.7, edgecolor='black', color='steelblue')
    ax2.set_xlabel('Good Pairs per Class', fontweight='bold')
    ax2.set_ylabel('Frequency', fontweight='bold')
    ax2.set_title('Good Pairs Distribution', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Keys distribution
    ax3 = fig.add_subplot(gs[0, 2])
    bins_k = range(min(num_keys), max(num_keys) + 2) if num_keys else 10
    ax3.hist(num_keys, bins=bins_k, alpha=0.7, edgecolor='black', color='coral')
    ax3.set_xlabel('Weak Keys per Class', fontweight='bold')
    ax3.set_ylabel('Frequency', fontweight='bold')
    ax3.set_title('Weak Keys Distribution', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Top 10 classes
    ax4 = fig.add_subplot(gs[1, :2])
    sorted_classes = sorted(classes, key=lambda x: x['num_keys'], reverse=True)[:10]
    labels = [f"{c['num_good_pairs']}gp" for c in sorted_classes]
    keys = [c['num_keys'] for c in sorted_classes]
    bars = ax4.barh(range(len(labels)), keys, color='mediumseagreen', edgecolor='black', linewidth=0.5)
    ax4.set_yticks(range(len(labels)))
    ax4.set_yticklabels(labels)
    ax4.set_xlabel('Weak Keys', fontweight='bold')
    ax4.set_title('Top 10 Classes', fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='x')
    
    # 5. Summary statistics
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis('off')
    stats_text = f"""SUMMARY STATISTICS
    
Total Classes: {summary.get('num_classes', 'N/A')}
Total Weak Keys: {summary.get('num_weak_keys', 'N/A')}
Total Keys: {summary.get('total_keys', 'N/A')}
Union Size: {summary.get('union_size', 'N/A')}

Good Pairs Stats:
  Mean: {np.mean(good_pairs):.1f}
  Median: {np.median(good_pairs):.0f}
  Min: {min(good_pairs)}
  Max: {max(good_pairs)}

Weak Keys Stats:
  Mean: {np.mean(num_keys):.1f}
  Median: {np.median(num_keys):.0f}
  Min: {min(num_keys)}
  Max: {max(num_keys)}
"""
    ax5.text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
             verticalalignment='center',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))
    
    fig.suptitle(f'Key Analysis Dashboard - {filename_prefix}', fontsize=16, fontweight='bold')
    
    plt.savefig(f'{filename_prefix}_dashboard.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_prefix}_dashboard.pdf', bbox_inches='tight')
    print(f"✓ Saved combined dashboard: {filename_prefix}_dashboard.png/pdf")
    plt.close()

def print_statistics(classes, params, summary):
    """Print detailed statistics"""
    good_pairs = [c['num_good_pairs'] for c in classes]
    num_keys = [c['num_keys'] for c in classes]
    
    print("\n" + "="*70)
    print("STATISTICAL ANALYSIS")
    print("="*70)
    
    print(f"\nAttack Parameters:")
    print(f"  Offset: {params.get('offset', 'N/A')}")
    print(f"  Active indices: {params.get('active', 'N/A')}")
    print(f"  dy1: {params.get('dy1', 'N/A')}, dy2: {params.get('dy2', 'N/A')}")
    
    print(f"\nOverall Statistics:")
    print(f"  Total number of classes: {len(classes)}")
    print(f"  Total weak keys: {summary.get('num_weak_keys', 'N/A')}")
    print(f"  Total key space: {summary.get('total_keys', 'N/A')}")
    print(f"  Union of good pairs: {summary.get('union_size', 'N/A')}")
    
    print(f"\nGood Pairs Statistics:")
    print(f"  Mean: {np.mean(good_pairs):.2f}")
    print(f"  Median: {np.median(good_pairs):.0f}")
    print(f"  Std Dev: {np.std(good_pairs):.2f}")
    print(f"  Min: {min(good_pairs)}, Max: {max(good_pairs)}")
    print(f"  Range: {max(good_pairs) - min(good_pairs)}")
    
    print(f"\nWeak Keys per Class Statistics:")
    print(f"  Mean: {np.mean(num_keys):.2f}")
    print(f"  Median: {np.median(num_keys):.0f}")
    print(f"  Std Dev: {np.std(num_keys):.2f}")
    print(f"  Min: {min(num_keys)}, Max: {max(num_keys)}")
    print(f"  Range: {max(num_keys) - min(num_keys)}")
    
    # Class size distribution
    counter = Counter(num_keys)
    print(f"\nClass Size Distribution (top 10):")
    for size, count in counter.most_common(10):
        print(f"  {size} keys: {count} classes")
    
    print("\n" + "="*70)

def extract_example_class(classes, params, num_keys_target=1, output_file=None):
    """Extract and format an example weak-key class for paper inclusion
    
    Args:
        classes: List of weak-key class dictionaries
        params: Attack parameters
        num_keys_target: Target number of weak keys in the example class (default: 1)
        output_file: Optional file to write the example to
    """
    # Find a class with the target number of keys
    example_class = None
    for cls in classes:
        if cls['num_keys'] == num_keys_target:
            example_class = cls
            break
    
    if not example_class:
        print(f"Warning: No class found with exactly {num_keys_target} weak key(s)")
        # Try to find the smallest class
        example_class = min(classes, key=lambda x: x['num_keys'])
        print(f"Using class with {example_class['num_keys']} weak keys instead")
    
    # Parse the good pairs and keys
    good_pairs = example_class.get('good_pairs')
    if not good_pairs:
        good_pairs_pattern = r"\(0x([0-9a-fA-F]{2}),\s*0x([0-9a-fA-F]{2})\)"
        good_pairs = [(p[0].lower(), p[1].lower())
                      for p in re.findall(good_pairs_pattern, example_class['good_pairs_str'])]

    keys = example_class.get('keys')
    if not keys:
        keys_pattern_new = r"left\s*=\s*0x([0-9a-fA-F]{2}),\s*right\s*=\s*0x([0-9a-fA-F]{2})"
        keys = [
            (k[0].lower(), k[1].lower())
            for k in re.findall(keys_pattern_new, example_class['keys_str'])
        ]
        if not keys:
            keys_pattern_old = r"\(0x([0-9a-fA-F]{2}),\s*0x([0-9a-fA-F]{2})\)"
            keys = [
                (k[0].lower(), k[1].lower())
                for k in re.findall(keys_pattern_old, example_class['keys_str'])
            ]
    
    # Format output
    output = []
    output.append("="*70)
    output.append("EXAMPLE WEAK-KEY CLASS FOR PAPER")
    output.append("="*70)
    output.append(f"\nParameters:")
    output.append(f"  Offset: {params.get('offset', 'N/A')}")
    output.append(f"  Active S-boxes: {params.get('active', 'N/A')}")
    output.append(f"  Δ₁: {params.get('dy1', 'N/A')}, Δ₂: {params.get('dy2', 'N/A')}")
    output.append(f"\nWeak-Key Class:")
    output.append(f"  Number of good pairs: {example_class['num_good_pairs']}")
    output.append(f"  Number of weak keys: {example_class['num_keys']}")
    output.append(f"\nGood Pairs:")
    for i, (p1, p2) in enumerate(good_pairs, 1):
        output.append(f"  {i}. (0x{p1}, 0x{p2})")
    output.append(f"\nWeak Keys (K₀, K₁):")
    for i, (k0, k1) in enumerate(keys, 1):
        output.append(f"  {i}. (0x{k0}, 0x{k1})")
    output.append("="*70)
    
    result = "\n".join(output)
    
    # Print to console
    print("\n" + result)
    
    # Optionally write to file
    if output_file:
        with open(output_file, 'w') as f:
            f.write(result + "\n")
        print(f"\n✓ Example saved to: {output_file}")
    
    return example_class, good_pairs, keys

def main():
    parser = argparse.ArgumentParser(description='Visualize key analysis results')
    parser.add_argument('input_file', help='Results file to analyze')
    parser.add_argument('--output-prefix', '-o', help='Output filename prefix (default: derived from input)')
    parser.add_argument('--plot-all', action='store_true', help='Generate all plots')
    parser.add_argument('--scatter', action='store_true', help='Generate scatter plot')
    parser.add_argument('--histogram', action='store_true', help='Generate histograms')
    parser.add_argument('--heatmap', action='store_true', help='Generate heatmap')
    parser.add_argument('--top-classes', action='store_true', help='Generate top classes bar chart')
    parser.add_argument('--frequency', action='store_true', help='Generate frequency table')
    parser.add_argument('--dashboard', action='store_true', help='Generate combined dashboard')
    parser.add_argument('--no-stats', action='store_true', help='Don\'t print statistics')
    parser.add_argument('--hideconfig', action='store_true', help='Hide offset and active S-boxes from legend')
    parser.add_argument('--example', action='store_true', help='Extract example weak-key class for paper')
    parser.add_argument('--example-keys', type=int, default=1, help='Target number of keys for example class (default: 1)')
    parser.add_argument('--example-output', help='File to save example to')
    
    args = parser.parse_args()
    
    # Parse input file
    print(f"\nParsing {args.input_file}...")
    params, classes, summary = parse_results_file(args.input_file)
    print(f"✓ Found {len(classes)} weak key classes")
    
    # Determine output prefix
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        prefix = args.input_file.replace('.txt', '').replace('results_', 'viz_')
    
    # Print statistics
    if not args.no_stats:
        print_statistics(classes, params, summary)
    
    # Generate plots
    print("\nGenerating visualizations...")
    
    if args.plot_all or args.scatter:
        plot_keys_vs_good_pairs_scatter(classes, params, summary, prefix)
    
    if args.plot_all or args.histogram:
        plot_good_pairs_distribution(classes, params, summary, prefix)
        plot_keys_distribution(classes, params, summary, prefix)
    
    if args.plot_all or args.heatmap:
        plot_heatmap(classes, params, summary, prefix)
    
    if args.plot_all or args.top_classes:
        plot_top_classes_bar(classes, params, summary, prefix)
    
    if args.plot_all or args.frequency:
        plot_frequency_table(classes, params, summary, prefix, hide_config=args.hideconfig)
    
    if args.plot_all or args.dashboard:
        plot_combined_dashboard(classes, params, summary, prefix)
    
    # If no specific plot requested, generate frequency table by default
    if not any([args.plot_all, args.scatter, args.histogram, args.heatmap, args.top_classes, args.frequency, args.dashboard]):
        plot_frequency_table(classes, params, summary, prefix, hide_config=args.hideconfig)
    
    # Extract example class if requested
    if args.example:
        extract_example_class(classes, params, num_keys_target=args.example_keys, output_file=args.example_output)
    
    print("\n✓ Visualization complete!")

if __name__ == '__main__':
    main()
