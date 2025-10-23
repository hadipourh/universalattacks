#!/usr/bin/env python3
"""
Visualization & Statistical Classification Tool for Differential-Linear Distinguisher

Features:
  - Dynamic half-normal threshold derivation.
  - Multi-tier key strength classification (WEAK / PROBABLE-WEAK / POSSIBLE-WEAK / STRONG).
  - Optional auto-tuning of classification parameters based on dataset size.
  - Publication-quality multi-panel figure (bar, scatter, histogram, peak spectrum).
  - Comparison, tail exceedance, and timeline visualizations.

Copyright (C) 2025 Hosein Hadipour & Mostafizar Rahman
License: GPLv3 or later
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import re
import math
from pathlib import Path

COLOR_BANDS = {
    'base': "#4B6CD6",          # deep navy for near-random behaviour
    'moderate': '#0EA5E9',      # bright cyan for 2.5σ–3σ region
    'conservative': '#F97316',  # vivid orange for ≥3σ
    'extreme': '#A855F7',       # saturated violet for ≥4σ
    'fwer': '#DC2626'           # bold crimson for FWER exceedance
}
try:
    from threshold_calculator import DynamicCorrelationThresholdCalculator
except ModuleNotFoundError:  # pragma: no cover - fallback path resolution
    import sys
    sys.path.append(str(Path(__file__).parent))
    from threshold_calculator import DynamicCorrelationThresholdCalculator

GOOD_PAIR_COLOR = '#FFD700'  # bright gold for emphasis

class CorrelationVisualizer:
    def __init__(self, data_dir="Data/R4_V1/Output10", makefile_path="makefile"):
        self.data_dir = Path(data_dir)
        self.makefile_path = makefile_path
        self.threshold_calc = DynamicCorrelationThresholdCalculator(makefile_path=makefile_path,
                                                                     data_dir=str(data_dir))
        self.thresholds = self.threshold_calc.get_visualization_thresholds()
        self.sqrt_n = self.threshold_calc.sqrt_n

    def _build_color_map(self, z_arr, good_mask, q_ext):
        colors = []
        for is_good, z in zip(good_mask, z_arr):
            if is_good:
                colors.append(GOOD_PAIR_COLOR)
            elif z >= q_ext:
                colors.append(COLOR_BANDS['fwer'])
            elif z >= 4.0:
                colors.append(COLOR_BANDS['extreme'])
            elif z >= 3.0:
                colors.append(COLOR_BANDS['conservative'])
            elif z >= 2.5:
                colors.append(COLOR_BANDS['moderate'])
            else:
                colors.append(COLOR_BANDS['base'])
        return np.array(colors, dtype=object)

    # ---------------- Utility / Parsing -----------------
    def auto_tune_params(self, s, q_ext):
        """Derive (R,G,R_hi,borderline_factor) heuristically from sample size s and theoretical extreme q_ext."""
        s = max(10, int(s))
        log_s = math.log10(s)
        q_norm = max(0.0, min(1.0, (q_ext - 3.0) / 2.5))
        R = 1.18 + 0.18 * q_norm + 0.04 * max(0.0, log_s - 3.0)
        R = max(1.1, min(R, 1.75))
        R_hi = max(R + 0.55, 1.85)
        G = 0.55 + 0.18 * log_s
        G = max(0.8, min(G, 1.8))
        borderline = 0.72 + 0.06 * max(0.0, 3.3 - log_s)
        borderline = max(0.65, min(borderline, 0.88))
        return {'R': R, 'G': G, 'R_hi': R_hi, 'borderline_factor': borderline}
        
    def parse_filename(self, filename):
        """Parse key information from filename like 'out_11-1_5-3.txt'."""
        basename = filename.split('/')[-1] if '/' in filename else filename

        match = re.match(r'out_(\d+)-(\d+)_(\d+)-(\d+)\.txt', basename)
        if match:
            left_hi = int(match.group(1))
            left_lo = int(match.group(2))
            right_hi = int(match.group(3))
            right_lo = int(match.group(4))
            return {
                'key1_val1': left_hi,
                'key1_val2': left_lo,
                'key2_val1': right_hi,
                'key2_val2': right_lo,
                'k_left': ((left_hi & 0xF) << 4) | (left_lo & 0xF),
                'k_right': ((right_hi & 0xF) << 4) | (right_lo & 0xF)
            }

        return {
            'key1_val1': '?',
            'key1_val2': '?', 
            'key2_val1': '?',
            'key2_val2': '?',
            'k_left': None,
            'k_right': None
        }
    
    def load_correlation_data(self, filename):
        """Load correlation data from output file."""
        filepath = self.data_dir / filename

        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        key_info = self.parse_filename(filename)
        header_pairs = set()
        header_key_left = None
        header_key_right = None
        records = []

        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith('#'):
                    if stripped.startswith('# KEY'):
                        match = re.search(r'left=0x([0-9a-fA-F]{2})\s+right=0x([0-9a-fA-F]{2})', stripped)
                        if match:
                            header_key_left = int(match.group(1), 16)
                            header_key_right = int(match.group(2), 16)
                    elif stripped.startswith('# GOOD_PAIR'):
                        hex_vals = re.findall(r'0x([0-9a-fA-F]{2})', stripped)
                        if len(hex_vals) >= 2:
                            gx = int(hex_vals[0], 16) & 0xFF
                            gy = int(hex_vals[1], 16) & 0xFF
                            header_pairs.add((min(gx, gy), max(gx, gy)))
                    continue

                parts = stripped.split()
                if len(parts) < 5:
                    continue
                try:
                    px_hi = int(parts[0]) & 0xF
                    py_hi = int(parts[1]) & 0xF
                    px_lo = int(parts[2]) & 0xF
                    py_lo = int(parts[3]) & 0xF
                    corr = int(parts[4])
                except ValueError:
                    continue

                x = ((px_hi & 0xF) << 4) | (px_lo & 0xF)
                y = ((py_hi & 0xF) << 4) | (py_lo & 0xF)
                canonical_pair = (min(x, y), max(x, y))
                is_good = canonical_pair in header_pairs if header_pairs else False

                records.append({
                    'line': line_num,
                    'corr': corr,
                    'x': x,
                    'y': y,
                    'canonical_pair': canonical_pair,
                    'is_good': is_good,
                    'px_hi': px_hi,
                    'py_hi': py_hi,
                    'px_lo': px_lo,
                    'py_lo': py_lo
                })

        if header_key_left is not None:
            key_info['k_left'] = header_key_left
        if header_key_right is not None:
            key_info['k_right'] = header_key_right

        return records, header_pairs, key_info
    
    # ---------------- Primary Single-File Visualization -----------------
    def create_correlation_plot(self, filename=None, save_plot=True, show_plot=True,
                                print_report=True, output_label=None):
        """Create comprehensive correlation visualization (single file).

        Steps:
          1. Load file & parse counts.
          2. Optional auto-tune of (R,G,R_hi,borderline).
          3. Classify using max/gap heuristic.
          4. Produce 4-panel figure (bar, peak spectrum, scatter, histogram).
          5. Save & report.
        """
        # 1. Resolve filename
        if filename is None:
            output_files = list(self.data_dir.glob('out_*.txt'))
            if not output_files:
                raise FileNotFoundError(f"No output files found in {self.data_dir}")
            filename = output_files[0].name

        records, header_pairs, key_info = self.load_correlation_data(filename)
        if not records:
            raise ValueError("No correlation data found in file")
        test_numbers = [rec['line'] for rec in records]
        corr_values = [rec['corr'] for rec in records]
        corr_arr = np.array(corr_values, dtype=float)
        sample_mean_corr = float(np.mean(corr_arr)) if corr_arr.size > 0 else 0.0
        sample_std_corr = float(np.std(corr_arr)) if corr_arr.size > 0 else 0.0
        test_numbers_arr = np.array(test_numbers, dtype=float)
        z_arr = corr_arr / self.sqrt_n
        sample_mean_z = sample_mean_corr / self.sqrt_n if self.sqrt_n > 0 else 0.0
        s = len(corr_arr)
        good_mask = np.array([rec['is_good'] for rec in records], dtype=bool)
        matched_good_pairs = len({rec['canonical_pair'] for rec in records if rec['is_good']})
        unique_good_pairs = len(header_pairs) if header_pairs else 0

        alpha = getattr(self, 'cli_alpha', 0.01)
        # Initial theoretical extreme (used only for tuning if requested)
        q_ext_initial = self._norm_ppf(1 - alpha / (2 * s))

        tuned = None
        if getattr(self, 'cli_auto_tune', False):
            tuned = self.auto_tune_params(s, q_ext_initial)
            R = tuned['R']; G = tuned['G']; R_hi = tuned['R_hi']; borderline = tuned['borderline_factor']
        else:
            R = getattr(self, 'cli_R', 3.0)
            G = getattr(self, 'cli_G', 1.2)
            R_hi = getattr(self, 'cli_R_hi', 6.0)
            borderline = getattr(self, 'cli_borderline_factor', 0.85)

        cls_result = self.classify_simple(corr_arr, alpha=alpha, R=R, G=G, R_hi=R_hi, borderline_factor=borderline)
        is_strong_label = cls_result['label'] == 'STRONG'
        highlight_good_mask = good_mask
        display_good_count = matched_good_pairs if matched_good_pairs else unique_good_pairs
        q_ext = cls_result['metrics']['q_ext']
        q_ext_threshold = q_ext * self.sqrt_n
        single_tail = 2 * (1 - 0.5 * (1 + math.erf(q_ext / math.sqrt(2))))
        single_tail = max(0.0, min(1.0, single_tail))
        estimated_fwer = 1 - max(0.0, (1 - single_tail)) ** len(corr_arr)

        fig = plt.figure(figsize=(20,14), dpi=150)
        ax_bar_corr = plt.subplot(2,2,1)
        ax_bar_thresh = plt.subplot(2,2,2)
        ax_scatter = plt.subplot(2,2,3)
        ax_hist = plt.subplot(2,2,4)

        # Color map shared between bar and scatter plots (strong keys omit "good pair" highlight)
        color_map = self._build_color_map(z_arr, highlight_good_mask, q_ext)

        # Derived label for y-axis (attempt power-of-two form)
        n_samples = int(self.sqrt_n ** 2)
        if n_samples > 0 and (n_samples & (n_samples - 1)) == 0:
            samples_label = f'$2^{{{int(math.log2(n_samples))}}} \\cdot | \\mathrm{{Cr}} |$'
        else:
            samples_label = '$N \\cdot | \\mathrm{Cr} |$'

        # Panel 1: Bar plot (thicken bars beyond 3σ for visual emphasis)
        bar_widths = np.full_like(z_arr, 0.9, dtype=float)
        if not is_strong_label:
            bar_widths[z_arr >= 3.0] = 1.3
        bars = ax_bar_corr.bar(
            test_numbers,
            corr_values,
            color=color_map,
            edgecolor='white',
            linewidth=0.4,
            alpha=1.0,
            width=bar_widths,
        )

        weak_mode = np.any(highlight_good_mask)
        for idx, patch in enumerate(bars):
            face = color_map[idx]
            patch.set_facecolor(face)
            patch.set_alpha(1.0)

            if highlight_good_mask[idx]:
                patch.set_edgecolor('#18181B')
                patch.set_linewidth(1.15 if weak_mode else 1.0)
                patch.set_zorder(3)
            else:
                patch.set_edgecolor(face)
                patch.set_linewidth(0.4 if weak_mode else 0.35)
                patch.set_zorder(1.6 if (not is_strong_label) else 1.5)
        ax_bar_corr.axhline(self.thresholds['random'], color='#2F4F4F', linestyle='-', linewidth=1.5, label='Expected mean (model)')
        ax_bar_corr.axhline(sample_mean_corr, color='#111827', linestyle=':', linewidth=1.3,
                            label=f'Sample mean ({sample_mean_corr:.1f})')
        ax_bar_corr.axhline(self.thresholds['moderate'], color=COLOR_BANDS['moderate'], linestyle='--', linewidth=1.6, label='Moderate (2.5σ)')
        ax_bar_corr.axhline(self.thresholds['conservative'], color=COLOR_BANDS['conservative'], linestyle='--', linewidth=1.6, label='Conservative (3σ)')
        ax_bar_corr.axhline(q_ext_threshold, color=COLOR_BANDS['fwer'], linestyle='--', linewidth=2.0,
                            label=f'FWER≈{estimated_fwer:.2e} (z={q_ext:.2f})')
        four_sigma_threshold = 4.0 * self.sqrt_n
        y_lower, y_upper = ax_bar_corr.get_ylim()
        if four_sigma_threshold < y_upper:
            gradient = np.linspace(0.0, 1.0, 256).reshape(-1, 1)
            transparent_red = mcolors.to_rgba(COLOR_BANDS['fwer'], alpha=0.0)
            accent_red = mcolors.to_rgba(COLOR_BANDS['fwer'], alpha=0.55)
            red_gradient = mcolors.LinearSegmentedColormap.from_list(
                'orthros_four_sigma_glow',
                [transparent_red, accent_red]
            )
            x_min = min(test_numbers) - 0.5
            x_max = max(test_numbers) + 0.5
            ax_bar_corr.imshow(
                gradient,
                aspect='auto',
                cmap=red_gradient,
                extent=[x_min, x_max, four_sigma_threshold, y_upper],
                origin='lower',
                interpolation='bicubic',
                zorder=0.5,
                alpha=1.0,
            )
            ax_bar_corr.set_ylim(y_lower, y_upper)
        ax_bar_corr.set_xlabel('Test Number')
        ax_bar_corr.set_ylabel(samples_label)
        ax_bar_corr.set_title('Correlation vs Test Number (Bar)', fontsize=11, fontweight='bold')
        ax_bar_corr.grid(True, alpha=0.25, axis='y')
        legend_handles, legend_labels = ax_bar_corr.get_legend_handles_labels()
        ax_bar_corr.legend(legend_handles, legend_labels, loc='upper left', fontsize=8,
                           frameon=True, facecolor='white', edgecolor='gray', framealpha=0.95)

        # Panel 2: Peak spectrum
        sorted_desc = np.sort(corr_arr)[::-1]
        K = min(50, len(sorted_desc))
        top_vals = sorted_desc[:K]
        ranks = np.arange(1, K+1)
        ax_bar_thresh.plot(ranks, top_vals, marker='o', linestyle='-', linewidth=1.2, markersize=3, color='tab:blue')
        if K > 5:
            with np.errstate(divide='ignore'):
                slope, intercept = np.polyfit(np.log(ranks), np.log(top_vals + 1e-9), 1)
                fitted = np.exp(intercept + slope * np.log(ranks))
            ax_bar_thresh.plot(ranks, fitted, linestyle='--', color='red', alpha=0.7, label=f'decay slope={slope:.2f}')
        ax_bar_thresh.annotate(f"max={int(top_vals[0])}", (1, top_vals[0]), textcoords='offset points', xytext=(6,4), fontsize=8)
        if K >=5:
            ax_bar_thresh.annotate(f"r5={int(top_vals[4])}", (5, top_vals[4]), textcoords='offset points', xytext=(6,4), fontsize=8)
        ax_bar_thresh.set_xlabel('Peak Rank')
        ax_bar_thresh.set_ylabel(samples_label)
        ax_bar_thresh.set_title('Peak Spectrum', fontsize=11, fontweight='bold')
        ax_bar_thresh.grid(True, alpha=0.3)
        ax_bar_thresh.axhline(sample_mean_corr, color='#111827', linestyle=':', linewidth=1.0,
                              label=f'Sample mean ({sample_mean_corr:.1f})')
        ax_bar_thresh.legend(frameon=False, fontsize=8)

        # Panel 3: Scatter
        scatter_sizes = np.full_like(z_arr, 22.0, dtype=float)
        scatter_sizes[z_arr >= 3.0] = 30.0
        scatter_sizes[z_arr >= 4.0] = 40.0
        scatter_sizes[z_arr >= q_ext] = 48.0
        non_good_mask = ~highlight_good_mask
        ax_scatter.scatter(
            test_numbers_arr[non_good_mask],
            corr_arr[non_good_mask],
            c=color_map[non_good_mask],
            s=scatter_sizes[non_good_mask],
            edgecolor='white',
            linewidth=0.45,
            alpha=1.0,
        )
        if np.any(highlight_good_mask):
            ax_scatter.scatter(
                test_numbers_arr[highlight_good_mask],
                corr_arr[highlight_good_mask],
                marker='*',
                s=scatter_sizes[highlight_good_mask] * 4.0,
                facecolors=GOOD_PAIR_COLOR,
                edgecolors='#18181B',
                linewidths=1.2,
                alpha=0.95,
                label='True good pair'
            )
        peak_mask = z_arr >= 3.0
        if np.any(peak_mask):
            highlight_mask = peak_mask & non_good_mask
            if np.any(highlight_mask):
                ax_scatter.scatter(
                    test_numbers_arr[highlight_mask],
                    corr_arr[highlight_mask],
                    facecolors='none',
                    edgecolors='#111827',
                    linewidths=1.0,
                    s=scatter_sizes[highlight_mask] * 1.15,
                    zorder=4,
                )
        running_max = np.maximum.accumulate(corr_arr)
        ax_scatter.plot(test_numbers, running_max, color='#1F2937', linewidth=1.1, alpha=0.75, label='Running max')
        ax_scatter.axhline(self.thresholds['random'], color='#4B5563', linestyle='-', linewidth=1.4, label='Expected mean (model)')
        ax_scatter.axhline(sample_mean_corr, color='#111827', linestyle=':', linewidth=1.2,
                           label=f'Sample mean ({sample_mean_corr:.1f})')
        ax_scatter.axhline(self.thresholds['moderate'], color=COLOR_BANDS['moderate'], linestyle='--', linewidth=1.6, label='Moderate (2.5σ)')
        ax_scatter.axhline(self.thresholds['conservative'], color=COLOR_BANDS['conservative'], linestyle='--', linewidth=1.6, label='Conservative (3σ)')
        ax_scatter.axhline(q_ext_threshold, color=COLOR_BANDS['fwer'], linestyle='--', linewidth=2.1,
                           label=f'FWER≈{estimated_fwer:.2e} (z={q_ext:.2f})')
        k_label = int(getattr(self, 'cli_label_top_k', 1))
        if k_label > 0:
            order = np.argsort(corr_arr)[::-1][:k_label]
            for idxp in order:
                x = test_numbers[idxp]; y = corr_arr[idxp]; z = z_arr[idxp]
                ax_scatter.annotate(f"{int(y)} (z={z:.2f})", (x,y), textcoords='offset points', xytext=(0,8), ha='center', fontsize=8)
        ax_scatter.set_xlabel('Test Number')
        ax_scatter.set_ylabel(samples_label)
        ax_scatter.set_title('Correlation Scatter', fontsize=11, fontweight='bold')
        ax_scatter.grid(True, alpha=0.25)
        ax_scatter.legend(loc='center left', bbox_to_anchor=(1.02,0.5), fontsize=8, frameon=True, facecolor='white', edgecolor='gray', framealpha=0.95)

        # Panel 4: Histogram
        counts, bins, _ = ax_hist.hist(corr_values, bins=30, color='#4B6CD6', alpha=0.7, edgecolor='white')
        ax_hist.axvline(self.thresholds['random'], color='#2F4F4F', linestyle='-', linewidth=2.0, label='Expected mean (model)')
        ax_hist.axvline(sample_mean_corr, color='#111827', linestyle=':', linewidth=1.8,
                        label=f'Sample mean ({sample_mean_corr:.1f})')
        ax_hist.axvline(self.thresholds['moderate'], color='#C06014', linestyle='--', linewidth=2.0, label='Moderate (2.5σ)')
        ax_hist.axvline(self.thresholds['conservative'], color='#8B4513', linestyle='--', linewidth=2.0, label='Conservative (3σ)')
        ax_hist.axvline(q_ext_threshold, color='#5C0000', linestyle='--', linewidth=2.5,
                        label=f'FWER≈{estimated_fwer:.2e} (z={q_ext:.2f})')
        ax_hist.set_xlabel('Correlation Count')
        ax_hist.set_ylabel('Frequency (log scale)')
        ax_hist.set_yscale('log', nonpositive='clip')
        ax_hist.set_title('Distribution', fontsize=11, fontweight='bold')
        ax_hist.grid(True, alpha=0.3)
        # Summary
        above_liberal = sum(1 for c in corr_values if c >= self.thresholds['liberal'])
        above_moderate = sum(1 for c in corr_values if c >= self.thresholds['moderate'])
        above_conservative = sum(1 for c in corr_values if c >= self.thresholds['conservative'])
        total_tests = len(corr_values)
        count_extreme = cls_result['metrics']['count_extreme']
        thresh_extreme = cls_result['metrics']['threshold_extreme']
        sample_pairs = int(round(self.sqrt_n ** 2))
        summary_text = (
            f"Tests={total_tests}\n"
            f"Pairs/Test (N) = {sample_pairs:,} (√N={self.sqrt_n:.0f})\n"
            f"Mean={sample_mean_corr:.1f} (theory {self.thresholds['random']:.1f})\n"
            f"Std={sample_std_corr:.1f}\n"
            f"Max={int(np.max(corr_arr))} (z={np.max(z_arr):.2f})\n"
            f">=z2: {above_liberal} ({above_liberal/total_tests*100:.1f}%)\n"
            f">=z2.5: {above_moderate} ({above_moderate/total_tests*100:.1f}%)\n"
            f">=z3: {above_conservative} ({above_conservative/total_tests*100:.1f}%)\n"
            f"≥z_extreme({thresh_extreme:.2f}): {count_extreme} ({count_extreme/total_tests*100:.1f}%)\n"
            f"True good pairs: {display_good_count}"
        )
        ax_hist.text(0.98,0.98, summary_text, transform=ax_hist.transAxes, ha='right', va='top', fontfamily='monospace', fontsize=9,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, pad=0.5))

        # Layout
        plt.tight_layout(rect=[0,0.03,1,0.96])

        # Save
        if save_plot:
            if output_label:
                base_name = output_label
            else:
                parts = []
                if hasattr(self, 'exp_params') and self.exp_params:
                    ep = self.exp_params
                    if ep.get('trail_name'): parts.append(ep['trail_name'])
                    if ep.get('rounds'): parts.append(f"{ep['rounds']}r")
                    if ep.get('offset'): parts.append(f"o{ep['offset']}")
                    if ep.get('deg'): parts.append(f"deg{ep['deg']}")
                key_nibbles = f"k{key_info.get('key1_val1','X')}-{key_info.get('key1_val2','X')}_{key_info.get('key2_val1','X')}-{key_info.get('key2_val2','X')}"
                parts.append(key_nibbles)
                max_corr = int(np.max(corr_arr)); z_max = np.max(z_arr)
                parts.append(f"{cls_result['label']}_max{max_corr}_z{z_max:.2f}")
                base_name = "_".join(parts)
            pdf_filename = f"{base_name}.pdf"; png_filename = f"{base_name}_analysis.png"
            fig.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=600, metadata={'Creator':'ORTHROS Key Recovery Analysis'})
            fig.savefig(png_filename, format='png', bbox_inches='tight', dpi=600, facecolor='white', edgecolor='none')
            print("\nHigh-quality plots saved:")
            print(f"  PDF (vector): {pdf_filename}")
            print(f"  PNG (600 DPI): {png_filename}")

        # Report
        if print_report:
            if tuned:
                old_R, old_G, old_R_hi, old_border = (getattr(self,'cli_R',None), getattr(self,'cli_G',None), getattr(self,'cli_R_hi',None), getattr(self,'cli_borderline_factor',None))
                self.cli_R, self.cli_G, self.cli_R_hi, self.cli_borderline_factor = R,G,R_hi,borderline
                self.print_analysis_report(corr_values, filename, key_info, good_pair_count=display_good_count)
                self.cli_R, self.cli_G, self.cli_R_hi, self.cli_borderline_factor = old_R, old_G, old_R_hi, old_border
            else:
                self.print_analysis_report(corr_values, filename, key_info, good_pair_count=display_good_count)
        if show_plot:
            plt.show()
        return fig
    
    # ---------------- Reporting -----------------
    def print_analysis_report(self, corr_values, filename, key_info, good_pair_count=0):
        """Print detailed statistical analysis for a single output file.

        Parameters
        ----------
        corr_values : list[int]
            Sequence of absolute difference counts |c0-c1| per test.
        filename : str
            Name of the file analyzed (out_*.txt).
        key_info : dict
            Parsed key nibble values.
        """

        print("DIFFERENTIAL-LINEAR DISTINGUISHER CORRELATION ANALYSIS")
        print("=" * 80)
        print(f"File: {filename}")
        print(
            f"Key Hypothesis: Nibble1=({key_info.get('key1_val1', '?')},{key_info.get('key1_val2', '?')}), "
            f"Nibble2=({key_info.get('key2_val1', '?')},{key_info.get('key2_val2', '?')})"
        )
        if good_pair_count:
            print(f"True good pairs for this key: {good_pair_count}")

        # Basic statistics
        print(f"\nBASIC STATISTICS:")
        print(f"  Total tests: {len(corr_values)}")
        print(f"  Mean correlation: {np.mean(corr_values):.2f}")
        print(f"  Standard deviation: {np.std(corr_values):.2f}")
        print(f"  Min/Max: {np.min(corr_values)} / {np.max(corr_values)}")
        print(f"  Half-normal mean: {self.thresholds['random']:.2f} (sqrt(n)={self.sqrt_n:.2f})")

        # Threshold analysis
        above_liberal = sum(1 for c in corr_values if c >= self.thresholds['liberal'])
        above_moderate = sum(1 for c in corr_values if c >= self.thresholds['moderate'])
        above_conservative = sum(1 for c in corr_values if c >= self.thresholds['conservative'])
        total_tests = len(corr_values)

        print(f"\nTHRESHOLD ANALYSIS:")
        print(
            f"  Above liberal ({self.thresholds['liberal']:.0f}):     {above_liberal:4d} "
            f"({above_liberal/total_tests*100:.1f}%)"
        )
        print(
            f"  Above moderate ({self.thresholds['moderate']:.0f}):    {above_moderate:4d} "
            f"({above_moderate/total_tests*100:.1f}%)"
        )
        print(
            f"  Above conservative ({self.thresholds['conservative']:.0f}): {above_conservative:4d} "
            f"({above_conservative/total_tests*100:.1f}%)"
        )

        # Statistical significance of max correlation
        max_corr = max(corr_values)
        analysis = self.threshold_calc.analyze_correlation_strength(max_corr)
        print(f"\nMAX CORRELATION ANALYSIS (Half-Normal):")
        print(f"  Value: {max_corr}")
        print(f"  z = value/√n = {analysis['z_score']:.2f}")
        print(f"  Half-normal tail p-value: {analysis['p_value_half_normal']:.6e}")

        # Simple classification output
        classification = self.classify_simple(
            np.array(corr_values),
            alpha=getattr(self, 'cli_alpha', 0.01),
            R=getattr(self, 'cli_R', 3.0),
            G=getattr(self, 'cli_G', 1.2),
            R_hi=getattr(self, 'cli_R_hi', 6.0),
            borderline_factor=getattr(self, 'cli_borderline_factor', 0.85)
        )
        print("\nCLASSIFICATION (SIMPLE MAX+GAP RULE):")
        print(f"  Label: {classification['label']}")
        print(f"  z_max = {classification['metrics']['z_max']:.2f}")
        print(f"  z_second = {classification['metrics']['z_second']:.2f}")
        print(f"  q_ext (theoretical extreme) = {classification['metrics']['q_ext']:.2f}")
        print(f"  ratio = z_max / q_ext = {classification['metrics']['ratio']:.2f}")
        print(f"  gap = z_max - z_second = {classification['metrics']['gap']:.2f}")
        print(f"  Peaks ≥ z_extreme (≥{classification['metrics']['threshold_extreme']:.2f}): {classification['metrics']['count_extreme']}")
        print(f"  Peaks ≥ dense threshold (≥{classification['metrics']['threshold_dense']:.2f}): {classification['metrics']['count_dense']}")
        print(f"  Peaks ≥ 3σ: {classification['metrics']['count_three_sigma']}")
        topk = min(5, len(corr_values))
        print(f"  Mean of top-{topk} z-scores: {classification['metrics']['topk_mean']:.2f}")
        print(f"  Parameters: α={classification['params']['alpha']}, R={classification['params']['R']}, G={classification['params']['G']}, R_hi={classification['params']['R_hi']}")
        if classification['params']['borderline_factor'] > 0:
            print(f"  Borderline factor={classification['params']['borderline_factor']:.2f} (target ratio {classification['params']['borderline_target']:.2f})")
        else:
            print(f"  Borderline evaluation disabled (factor ≤ 0)")
        print("  Decision rationale:")
        print(f"    {classification['rationale']}")
        print("=" * 80)

    # ------------------- Classification Framework -------------------
    def _norm_ppf(self, p):
        """Standard-library inverse CDF for N(0,1) using statistics.NormalDist."""
        try:
            from statistics import NormalDist
            return NormalDist().inv_cdf(p)
        except Exception:
            # As a last resort, fall back to numpy/scipy if present, else raise
            try:
                import numpy as np
                import math as _m
                # Use erfc-based binary search (rare path; avoid manual coeffs)
                # Solve Phi(x)=p where Phi(x)=0.5*(1+erf(x/sqrt(2)))
                lo, hi = -10.0, 10.0
                for _ in range(80):
                    mid = 0.5*(lo+hi)
                    cdf = 0.5*(1.0 + _m.erf(mid/_m.sqrt(2.0)))
                    if cdf < p:
                        lo = mid
                    else:
                        hi = mid
                return 0.5*(lo+hi)
            except Exception:
                raise RuntimeError("No standard inverse normal available; requires Python 3.8+ for statistics.NormalDist.")

    def classify_simple(self, corr_arr, alpha=0.01, R=3.0, G=1.2, R_hi=6.0, borderline_factor=0.85):
        """Simple max+gap classifier.

        Statistical model & rationale
        -----------------------------
        Under a *strong* (null) key hypothesis we model each absolute difference
        count |c0-c1| as approximately HalfNormal(sigma = sqrt(n)), i.e. if
        Z_i = |c0-c1|_i / sqrt(n) then Z_i ~ |N(0,1)| independently (idealized).

        Let M_s = max_{1..s} Z_i. For small family-wise error rate alpha, solve
            P(M_s <= q_ext) = 1 - alpha  =>  P(M_s > q_ext) = alpha
        Since P(Z > t) = 2(1 - Phi(t)) for half-normal (two tails of N(0,1)),
            P(M_s > t) ≈ 1 - (1 - 2(1-Phi(t)))^s ≈ 2s(1-Phi(t))  (union bound)
        Set 2 s (1 - Phi(q_ext)) = alpha  ->  Phi(q_ext) = 1 - alpha/(2s)
        Thus q_ext = Phi^{-1}(1 - alpha/(2s)).

                Decision logic (multi-tier):
                    1. Automatic WEAK if ratio = z_max / q_ext ≥ R_hi (extreme outlier/high-tier).
                    2. Else WEAK if ratio ≥ R and isolation gap (z_max - z_second) ≥ G (isolated spike).
                    3. Else WEAK if sustained extreme activity: ≥3 peaks above the extreme threshold or ≥2 peaks with high average.
                    4. Else PROBABLE-WEAK when large-but-crowded outliers or dense high peaks are present.
                    5. Else POSSIBLE-WEAK if ratio > 1 and either (a) ratio exceeds a softened borderline target or (b) multiple ≥3σ peaks exist.
                    6. Else STRONG.

        Parameters
        ----------
        corr_arr : array-like
            Raw |c0-c1| counts.
        alpha : float
            Target family-wise error probability for observing a spurious max.
        R : float
            Multiplicative factor beyond theoretical extreme to require (base tier).
        G : float
            Required absolute isolation gap between top two z-scores (base tier).
        R_hi : float
            Higher ratio triggering immediate WEAK (no gap needed).
        """
        if len(corr_arr) == 0:
            return {
                'label': 'STRONG', 'label_icon': '',
                'metrics': {'z_max':0,'q_ext':0,'ratio':0,'gap':0},
                'params': {'alpha':alpha,'R':R,'G':G,'R_hi':R_hi,'borderline_factor':borderline_factor},
                'rationale': 'Empty input.'
            }
        z = corr_arr / self.sqrt_n
        s = len(z)
        # Extreme quantile for max of half-normal samples: solve 2s(1-Φ(q)) = α
        target = 1 - alpha / (2*s)
        q_ext = self._norm_ppf(target)
        z_sorted = np.sort(z)[::-1]
        z_max = float(z_sorted[0])
        z_second = float(z_sorted[1]) if s > 1 else 0.0
        gap = z_max - z_second
        ratio = z_max / q_ext if q_ext > 0 else float('inf')

        # Multi-peak diagnostics
        top_k = min(5, s)
        top_k_mean = float(np.mean(z_sorted[:top_k])) if top_k > 0 else 0.0
        z_extreme_threshold = max(q_ext, 3.0)
        z_dense_threshold = max(2.8, z_extreme_threshold * 0.92)
        count_extreme = int(np.sum(z >= z_extreme_threshold))
        count_dense = int(np.sum(z >= z_dense_threshold))
        count_three_sigma = int(np.sum(z >= 3.0))

        multi_peak_extreme = (count_extreme >= 3) or (count_extreme >= 2 and top_k_mean >= 0.96 * z_extreme_threshold)
        multi_peak_dense = (count_dense >= 5) or (count_dense >= 4 and top_k_mean >= 0.9 * z_extreme_threshold)

        borderline_target = max(1.05, borderline_factor * max(R, 1.15)) if borderline_factor > 0 else float('inf')

        if ratio >= R_hi:
            label = 'WEAK'
            icon = ''
            rationale = f'Extreme max ({z_max:.2f}) is {ratio:.2f}× theoretical extreme; surpasses R_hi={R_hi:.2f}.'
        elif (ratio >= R) and (gap >= G):
            label = 'WEAK'
            icon = ''
            rationale = f'Isolated spike: ratio {ratio:.2f} ≥ R={R:.2f} with gap {gap:.2f} ≥ G={G:.2f}.'
        elif multi_peak_extreme:
            label = 'WEAK'
            icon = ''
            rationale = f'Sustained extreme activity: {count_extreme} peaks ≥ z_extreme (≥{z_extreme_threshold:.2f}).'
        elif (ratio >= R) or (count_extreme >= 2) or multi_peak_dense:
            label = 'PROBABLE-WEAK'
            icon = ''
            triggers = []
            if ratio >= R:
                triggers.append(f'ratio {ratio:.2f} ≥ R')
            if count_extreme >= 2:
                triggers.append(f'{count_extreme} peaks ≥ z_extreme')
            if multi_peak_dense and count_extreme < 2:
                triggers.append(f'dense high peaks (top{top_k} mean {top_k_mean:.2f})')
            rationale = ' and '.join(triggers) if triggers else 'Elevated peak cluster without clear isolation.'
        elif (ratio > 1.0) and ((ratio >= borderline_target) or (count_three_sigma >= 3) or (count_dense >= 3)):
            label = 'POSSIBLE-WEAK'
            icon = ''
            reasons = []
            if ratio >= borderline_target:
                reasons.append(f'ratio {ratio:.2f} ≥ borderline target {borderline_target:.2f}')
            if count_three_sigma >= 3:
                reasons.append(f'{count_three_sigma} peaks ≥ 3σ')
            elif count_dense >= 3:
                reasons.append(f'{count_dense} peaks above dense threshold (≥{z_dense_threshold:.2f})')
            rationale = '; '.join(reasons) if reasons else 'Moderate excess over theoretical extreme.'
        else:
            label = 'STRONG'
            icon = ''
            rationale = 'Peak behaviour consistent with null (max near theoretical extreme and no dense cluster).'

        return {
            'label': label,
            'label_icon': icon,
            'metrics': {
                'z_max': z_max,
                'z_second': z_second,
                'q_ext': q_ext,
                'ratio': ratio,
                'gap': gap,
                'count_extreme': count_extreme,
                'count_dense': count_dense,
                'count_three_sigma': count_three_sigma,
                'topk_mean': top_k_mean,
                'threshold_extreme': z_extreme_threshold,
                'threshold_dense': z_dense_threshold
            },
            'params': {
                'alpha': alpha,
                'R': R,
                'G': G,
                'R_hi': R_hi,
                'borderline_factor': borderline_factor,
                'borderline_target': borderline_target
            },
            'rationale': rationale
        }

    # ---------------- Multi-file Visualizations -----------------
    def create_key_comparison_plot(self, files, save_plot=True, max_cols=3):
        """Compare multiple key hypothesis outputs side-by-side.

        Parameters
        ----------
        files : list[str]
            Filenames (out_*.txt) to compare. Should include at least one strong-key (flat) and one or more weak-key (peaked) examples.
        save_plot : bool
            Whether to save the resulting figure.
        max_cols : int
            Maximum number of subplot columns before wrapping to next row.

        Layout
        ------
        Each subplot: scatter of counts vs test index with z-based color, running max, threshold lines, top peak annotated.
        Bottom row (if >1 row of scatters): aggregated bar chart of counts of exceedances per file at z in {2,2.5,3,4}.
        """
        if not files:
            raise ValueError("No files provided for comparison")

        # Load all data first
        data = []
        for fname in files:
            try:
                records, header_pairs, key_info = self.load_correlation_data(fname)
                if not records:
                    continue
                tests = [rec['line'] for rec in records]
                vals = np.array([rec['corr'] for rec in records], dtype=float)
                z_vals = vals / self.sqrt_n
                good_mask = np.array([rec['is_good'] for rec in records], dtype=bool)
                data.append((fname, tests, vals, z_vals, key_info, good_mask))
            except FileNotFoundError:
                print(f"Warning: file not found {fname}, skipping")
        if not data:
            raise ValueError("No valid data loaded for comparison")

        k = len(data)
        cols = min(max_cols, k)
        rows = int(np.ceil(k / cols))

        # Determine if we add an aggregate bar panel
        add_bar = k > 1
        fig_height = 3.0 * rows + (2.2 if add_bar else 0.6)
        fig = plt.figure(figsize=(5.5 * cols, fig_height))

        # Threshold lines (raw counts)
        z_lines = [2, 2.5, 3, 4]
        count_lines = {z: z * self.sqrt_n for z in z_lines}

        exceed_matrix = []  # rows per file, columns per z-line

        for idx, (fname, tests, vals, z_vals, key_info, good_mask) in enumerate(data, start=1):
            ax = fig.add_subplot(rows + (1 if add_bar else 0), cols, idx)
            # Color mapping by z bands
            colors = []
            for idxp, z in enumerate(z_vals):
                if good_mask[idxp]:
                    colors.append(GOOD_PAIR_COLOR)
                elif z >= 4:
                    colors.append('purple')
                elif z >= 3:
                    colors.append('red')
                elif z >= 2.5:
                    colors.append('orange')
                elif z >= 2:
                    colors.append('gold')
                else:
                    colors.append('lightgray')
            ax.scatter(tests, vals, c=colors, s=10, edgecolor='none', alpha=0.85)
            # Running max
            run_max = np.maximum.accumulate(vals)
            ax.plot(tests, run_max, color='black', linewidth=0.7, alpha=0.6)
            # Threshold lines
            for z, c in [(2,'gold'), (2.5,'orange'), (3,'red'), (4,'purple')]:
                ax.axhline(count_lines[z], color=c, linestyle='--', linewidth=0.9)
            # Annotate top peak
            peak_idx = int(np.argmax(vals))
            peak_val = vals[peak_idx]
            peak_z = z_vals[peak_idx]
            ax.annotate(f"{int(peak_val)}\nz={peak_z:.2f}", (tests[peak_idx], peak_val),
                        textcoords='offset points', xytext=(0,5), ha='center', fontsize=7)
            # Title with key info & simple classifier (consistent with single-file view)
            cls = self.classify_simple(
                vals,
                alpha=getattr(self, 'cli_alpha', 0.01),
                R=getattr(self, 'cli_R', 3.0),
                G=getattr(self, 'cli_G', 1.2),
                R_hi=getattr(self, 'cli_R_hi', 6.0),
                borderline_factor=getattr(self, 'cli_borderline_factor', 0.85)
            )
            key_str = f"({key_info.get('key1_val1','?')},{key_info.get('key1_val2','?')})-({key_info.get('key2_val1','?')},{key_info.get('key2_val2','?')})"
            ax.set_title(
                f"{fname}\n{key_str} {cls['label']} z={peak_z:.2f} r={cls['metrics']['ratio']:.2f} g={cls['metrics']['gap']:.2f}",
                fontsize=9
            )
            ax.set_xlabel('test')
            if idx % cols == 1:
                ax.set_ylabel('|c0-c1|')
            ax.grid(alpha=0.25)

            exceed_matrix.append([np.sum(z_vals >= zc) for zc in z_lines])

        # Aggregate bar chart for exceedances
        if add_bar:
            ax_bar = fig.add_subplot(rows + 1, 1, rows + 1)
            exceed_arr = np.array(exceed_matrix)
            width = 0.18
            x = np.arange(len(z_lines))
            for i, (fname, *_rest) in enumerate(data):
                ax_bar.bar(x + (i - k/2)*width + width/2, exceed_arr[i], width=width,
                           label=f"{fname}")
            ax_bar.set_xticks(x)
            ax_bar.set_xticklabels([f"z≥{z}" for z in z_lines])
            ax_bar.set_ylabel('count')
            ax_bar.set_title('Exceedance counts across keys')
            ax_bar.legend(fontsize=8, ncol=min(k,5))
            ax_bar.grid(axis='y', alpha=0.3)

        fig.suptitle('Key Strength Correlation Comparison', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0,0.03,1,0.96])

        if save_plot:
            # Save both PDF (vector) and high-res PNG
            base_name = 'key_correlation_comparison'
            
            pdf_name = f'{base_name}.pdf'
            fig.savefig(pdf_name, format='pdf', bbox_inches='tight', dpi=600)
            
            png_name = f'{base_name}.png'
            fig.savefig(png_name, format='png', bbox_inches='tight', dpi=600, 
                       facecolor='white', edgecolor='none')
            
            print(f"Saved comparison figures:")
            print(f"  PDF: {pdf_name}")
            print(f"  PNG (600 DPI): {png_name}")

        plt.show()
        return fig

    def _load_multiple(self, files):
        data = []
        for fname in files:
            try:
                records, header_pairs, key_info = self.load_correlation_data(fname)
                if not records:
                    continue
                tests = [rec['line'] for rec in records]
                vals = np.array([rec['corr'] for rec in records], dtype=float)
                z_vals = vals / self.sqrt_n
                good_mask = np.array([rec['is_good'] for rec in records], dtype=bool)
                data.append((fname, tests, vals, z_vals, key_info, good_mask))
            except FileNotFoundError:
                print(f"Warning: file not found {fname}, skipping")
        return data

    def create_peak_timeline_plot(self, files, z_cut=2.0, save_plot=True):
        """Plot only significant peak events (z >= z_cut) for multiple files on a timeline.

        Each file is a horizontal categorical band; marker color encodes z band.
        Annotates first z>=3 event and max peak per file when present.
        """
        data = self._load_multiple(files)
        if not data:
            raise ValueError("No valid data for peak timeline")

        # Sort by max z descending (weak-like first)
        data.sort(key=lambda x: np.max(x[3]), reverse=True)

        fig, ax = plt.subplots(figsize=(12, 1.1 + 0.65*len(data)))
        y_positions = np.arange(len(data))

        for i,(fname, tests, vals, z_vals, key_info, good_mask) in enumerate(data):
            mask = z_vals >= z_cut
            if not np.any(mask):
                # plot a faint baseline marker to keep row visible
                ax.plot([0],[y_positions[i]], marker='x', color='lightgray')
                continue
            indices = np.where(mask)[0]
            sel_tests = np.array(tests)[indices]
            sel_z = z_vals[indices]
            colors=[]
            for idx_mask, z in zip(indices, sel_z):
                if good_mask[idx_mask]:
                    colors.append(GOOD_PAIR_COLOR)
                elif z >= 4: colors.append('purple')
                elif z >= 3: colors.append('red')
                elif z >= 2.5: colors.append('orange')
                elif z >= 2: colors.append('gold')
                else: colors.append('lightgray')
            ax.scatter(sel_tests, np.full_like(sel_tests, y_positions[i]), c=colors, s=28, edgecolor='none', alpha=0.9)
            # annotate first z>=3
            high_mask = sel_z >= 3
            if np.any(high_mask):
                first_idx = np.argwhere(high_mask)[0][0]
                ax.annotate(f"first z>=3 @ {int(sel_tests[first_idx])}",
                            (sel_tests[first_idx], y_positions[i]), textcoords='offset points', xytext=(5,6),
                            fontsize=7, color='red')
            # annotate max
            peak_idx = np.argmax(z_vals)
            peak_test = tests[peak_idx]
            peak_z = z_vals[peak_idx]
            ax.annotate(f"max z={peak_z:.2f}", (peak_test, y_positions[i]), textcoords='offset points', xytext=(5,-10),
                        fontsize=7, color='black')
        ax.set_yticks(y_positions)
        ax.set_yticklabels([d[0] for d in data])
        ax.set_xlabel('Test index')
        ax.set_title(f'Peak Timeline (z≥{z_cut})')
        ax.grid(axis='x', alpha=0.25)
        # Color legend proxy
        from matplotlib.lines import Line2D
        legend_elems = [Line2D([0],[0], marker='o', color='w', label='2≤z<2.5', markerfacecolor='gold', markersize=6),
                        Line2D([0],[0], marker='o', color='w', label='2.5≤z<3', markerfacecolor='orange', markersize=6),
                        Line2D([0],[0], marker='o', color='w', label='3≤z<4', markerfacecolor='red', markersize=6),
                        Line2D([0],[0], marker='o', color='w', label='z≥4', markerfacecolor='purple', markersize=6)]
        ax.legend(handles=legend_elems, bbox_to_anchor=(1.02,1), loc='upper left', frameon=False, title='Bands')
        plt.tight_layout()
        if save_plot:
            fig.savefig('peak_timeline.pdf', format='pdf', bbox_inches='tight', dpi=600)
            fig.savefig('peak_timeline.png', format='png', bbox_inches='tight', dpi=600, 
                       facecolor='white', edgecolor='none')
            print('Saved peak_timeline.pdf and peak_timeline.png (600 DPI)')
        plt.show()
        return fig

    def create_tail_exceedance_plot(self, files, save_plot=True):
        """Plot empirical exceedance curves P(Z>=t) with theoretical half-normal tail and bar counts."""
        data = self._load_multiple(files)
        if not data:
            raise ValueError("No valid data for tail exceedance plot")

        # Build exceedance curves
        fig = plt.figure(figsize=(10, 6))
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(2,1, height_ratios=[2.7,1])
        ax_curve = fig.add_subplot(gs[0])
        ax_bar = fig.add_subplot(gs[1])

        z_lines = [2,2.5,3,4]
        exceed_counts = []
        labels = []

        # Theoretical tail: 2(1-Phi(t))
        t_grid = np.linspace(0,4.5,300)
        # Robust tail function: try scipy.special.erfc for vectorization; fallback to numpy frompyfunc
        try:
            from scipy.special import erfc
            def phi_tail(z):
                # Phi(z) = 1 - 0.5*erfc(z/√2); tail two-sided for |Z|
                return 2*(1 - (1 - 0.5*erfc(z/np.sqrt(2))))
        except Exception:
            import math
            erf_vec = np.frompyfunc(lambda x: math.erf(float(x)), 1, 1)
            def phi_tail(z):
                z = np.asarray(z, dtype=float)
                erf_vals = erf_vec(z/np.sqrt(2)).astype(float)
                return 2*(1 - 0.5*(1 + erf_vals))
        theory = phi_tail(t_grid)
        ax_curve.plot(t_grid, theory, color='black', linestyle='-', linewidth=1.2, label='Half-normal tail')

        for (fname, _tests, _vals, z_vals, _ki) in data:
            z_sorted = np.sort(z_vals)
            # Empirical survival: for each point z_i, S(z_i)= (# >= z_i)/N
            N = len(z_sorted)
            surv_z = z_sorted
            surv_p = (N - np.arange(N)) / N
            ax_curve.step(surv_z, surv_p, where='post', label=fname, linewidth=1.1)
            exceed_counts.append([np.sum(z_vals >= thr) for thr in z_lines])
            labels.append(fname)

        ax_curve.set_yscale('log')
        ax_curve.set_xlabel('z')
        ax_curve.set_ylabel('P(Z ≥ z)')
        ax_curve.set_title('Empirical Exceedance vs Half-Normal Tail')
        ax_curve.set_xlim(0,4.5)
        ax_curve.set_ylim(1e-4,1)
        ax_curve.grid(True, which='both', axis='both', alpha=0.3)
        ax_curve.legend(fontsize=8, ncol=2)

        # Bars
        exceed_arr = np.array(exceed_counts)
        x = np.arange(len(z_lines))
        width = 0.8 / max(1,len(labels))
        for i,label in enumerate(labels):
            ax_bar.bar(x + (i - len(labels)/2)*width + width/2, exceed_arr[i], width=width, label=label)
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels([f"z≥{z}" for z in z_lines])
        ax_bar.set_ylabel('count')
        ax_bar.set_title('Exceedance Counts')
        ax_bar.grid(axis='y', alpha=0.3)
        if len(labels) > 1:
            ax_bar.legend(fontsize=8, ncol=min(4,len(labels)))
        plt.tight_layout()
        if save_plot:
            fig.savefig('tail_exceedance.pdf', format='pdf', bbox_inches='tight', dpi=600)
            fig.savefig('tail_exceedance.png', format='png', bbox_inches='tight', dpi=600,
                       facecolor='white', edgecolor='none')
            print('Saved tail_exceedance.pdf and tail_exceedance.png (600 DPI)')
        plt.show()
        return fig

    # ------------------- Synthetic Data Utilities -------------------
    def generate_synthetic_null(self, s=5000, filename=None, seed=None):
        """Generate a synthetic strong-key (null) dataset and write an out_*.txt file.

        The synthetic model draws Z_i ~ |N(0,1)| independently, then rescales to
        counts C_i = round(Z_i * sqrt(n)). We fabricate the first four integer
        columns as zeros so that the 5th column matches expected correlation format.

        Parameters
        ----------
        s : int
            Number of tests (lines) to generate.
        filename : str or None
            Destination filename (within data_dir). If None, auto-name as
            out_SYNTH-NULL_0-0.txt (with numeric suffix if collision).
        seed : int or None
            Optional RNG seed for reproducibility.

        Returns
        -------
        str : actual filename created (basename only).
        dict : classification result for the synthetic dataset.
        """
        rng = np.random.default_rng(seed)
        # Half-normal via absolute of normal draws
        z = np.abs(rng.standard_normal(size=s))
        counts = np.round(z * self.sqrt_n).astype(int)

        # Ensure data directory exists
        self.data_dir.mkdir(parents=True, exist_ok=True)

        if filename is None:
            base = 'out_SYNTH-NULL_0-0'
            candidate = base
            counter = 1
            while (self.data_dir / f"{candidate}.txt").exists():
                candidate = f"{base}_{counter}"
                counter += 1
            filename = f"{candidate}.txt"
        else:
            if not filename.endswith('.txt'):
                filename += '.txt'

        path = self.data_dir / filename
        with open(path, 'w') as f:
            for c in counts:
                # Minimal 5-column line; 5th column is correlation count
                f.write(f"0 0 0 0 {int(c)}\n")

        cls = self.classify_simple(
            counts,
            alpha=getattr(self, 'cli_alpha', 0.01),
            R=getattr(self, 'cli_R', 3.0),
            G=getattr(self, 'cli_G', 1.2),
            R_hi=getattr(self, 'cli_R_hi', 6.0),
            borderline_factor=getattr(self, 'cli_borderline_factor', 0.85)
        )
        print(f"Generated synthetic null file: {filename} (s={s}) -> {cls['label']} ratio={cls['metrics']['ratio']:.2f} gap={cls['metrics']['gap']:.2f}")
        return filename, cls

def main():
    """Main function for visualization."""
    import argparse
    import os
    # We'll normalize data directory BEFORE possibly changing directories so that
    # user-supplied paths like "keyrecovery/Data/..." do not get duplicated when
    # we chdir into the script directory (which is already keyrecovery/).
    script_dir = Path(__file__).parent
    cwd_before = Path.cwd()
    
    parser = argparse.ArgumentParser(description='Visualize differential-linear distinguisher correlations')
    parser.add_argument('--data-dir', '-d', default='Data/R4_V1/Output10', 
                       help='Data directory containing output files')
    parser.add_argument('--file', '-f', help='Specific file to analyze')
    parser.add_argument('--makefile', '-m', default='makefile', help='Path to makefile')
    parser.add_argument('--no-save', action='store_true', help='Do not save plots')
    parser.add_argument('--compare', nargs='+', help='List of out_*.txt files to produce multi-key comparison figure')
    parser.add_argument('--timeline', nargs='+', help='List of out_*.txt files for peak timeline figure')
    parser.add_argument('--tails', nargs='+', help='List of out_*.txt files for tail exceedance figure')
    parser.add_argument('--synth-null', type=int, metavar='S', help='Generate synthetic null dataset with S tests and visualize')
    parser.add_argument('--alpha', type=float, default=0.01, help='Family-wise error rate target (default 0.01)')
    parser.add_argument('--R', type=float, default=1.30, help='Base ratio threshold (default 1.30; lower -> more sensitive)')
    parser.add_argument('--G', type=float, default=1.05, help='Isolation gap threshold (default 1.05)')
    parser.add_argument('--R-hi', dest='R_hi', type=float, default=2.00, help='High-tier ratio threshold (default 2.00)')
    parser.add_argument('--borderline-factor', type=float, default=0.72,
                        help='Fraction of R to flag POSSIBLE-WEAK (default 0.72; set 0 to disable)')
    parser.add_argument('--sensitivity', choices=['conservative','balanced','aggressive'], default=None,
                        help='Preset for (R,G,R_hi,borderline). Overrides individual flags if set.')
    parser.add_argument('--auto-tune', action='store_true',
                        help='Automatically tune R/G/R_hi/borderline based on number of tests & estimated extreme.')
    parser.add_argument('--label-top-k', type=int, default=1, help='Number of top peaks to label in scatter (default 1; 0 disables)')
    parser.add_argument('--batch-all', action='store_true', help='Generate and save plots for all out_*.txt files (no display)')
    parser.add_argument('--output-label', '-l', type=str, default=None, 
                       help='Custom label for output files (e.g., "weak_key_1" produces weak_key_1.pdf/png instead of out_X-Y_Z-W_analysis.pdf/png)')
    
    # Experiment parameter arguments for descriptive filenames
    parser.add_argument('--deg', type=int, help='DEG parameter value')
    parser.add_argument('--rounds', type=int, help='Number of rounds')
    parser.add_argument('--offset', type=int, help='Round offset')
    parser.add_argument('--output-mask', type=str, help='Output mask pattern')
    parser.add_argument('--input-pattern', type=str, help='Input difference pattern')
    parser.add_argument('--trail-name', type=str, help='Trail/experiment name')
    
    args = parser.parse_args()
    
    try:
        # Resolve absolute paths first (with respect to original cwd) then chdir.
        data_dir_abs = (cwd_before / args.data_dir).resolve() if not Path(args.data_dir).is_absolute() else Path(args.data_dir)
        makefile_abs = (cwd_before / args.makefile).resolve() if not Path(args.makefile).is_absolute() else Path(args.makefile)

        # Change working directory to script dir for any relative imports / ancillary files
        os.chdir(script_dir)

        visualizer = CorrelationVisualizer(
            data_dir=str(data_dir_abs),
            makefile_path=str(makefile_abs)
        )
        # Attach CLI threshold tuning parameters for downstream classification calls
        visualizer.cli_alpha = args.alpha
        visualizer.cli_R = args.R
        visualizer.cli_G = args.G
        visualizer.cli_R_hi = args.R_hi
        visualizer.cli_label_top_k = args.label_top_k
        # Apply sensitivity presets if requested
        if args.sensitivity:
            if args.sensitivity == 'conservative':
                args.R, args.G, args.R_hi, args.borderline_factor = 1.45, 1.20, 2.20, 0.80
            elif args.sensitivity == 'balanced':
                args.R, args.G, args.R_hi, args.borderline_factor = 1.30, 1.05, 2.00, 0.72
            elif args.sensitivity == 'aggressive':
                args.R, args.G, args.R_hi, args.borderline_factor = 1.18, 0.90, 1.80, 0.65
        visualizer.cli_borderline_factor = args.borderline_factor
        visualizer.cli_auto_tune = args.auto_tune
        
        # Store experiment parameters for descriptive filenames
        visualizer.exp_params = {
            'deg': args.deg,
            'rounds': args.rounds,
            'offset': args.offset,
            'output_mask': args.output_mask,
            'input_pattern': args.input_pattern,
            'trail_name': args.trail_name
        }
        
        # Mutually exclusive precedence: timeline > tails > compare > single
        if args.batch_all:
            files = sorted([p.name for p in Path(str(data_dir_abs)).glob('out_*.txt')])
            if not files:
                print('No out_*.txt files found for batch.')
                return 1
            print(f'Batch generating plots for {len(files)} files...')
            for f in files:
                try:
                    visualizer.create_correlation_plot(filename=f, save_plot=True, show_plot=False, print_report=False)
                except Exception as e:
                    print(f'  [WARN] Failed {f}: {e}')
            print('Batch generation complete.')
        elif args.synth_null:
            fname, _cls = visualizer.generate_synthetic_null(s=args.synth_null)
            visualizer.create_correlation_plot(filename=fname, save_plot=not args.no_save, show_plot=not args.no_save, output_label=args.output_label)
        elif args.timeline:
            visualizer.create_peak_timeline_plot(args.timeline, save_plot=not args.no_save)
        elif args.tails:
            visualizer.create_tail_exceedance_plot(args.tails, save_plot=not args.no_save)
        elif args.compare:
            visualizer.create_key_comparison_plot(args.compare, save_plot=not args.no_save)
        else:
            visualizer.create_correlation_plot(filename=args.file, save_plot=not args.no_save, show_plot=not args.no_save, output_label=args.output_label)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
