#!/usr/bin/env python3
"""
Copyright (C) 2025 Hosein Hadipour
Email: hsn.hadipour@gmail.com
Date: September 24, 2025

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

This program creates publication-ready visualizations for correlation and probability datasets

Features:
- Intelligent bin selection using multiple statistical rules
- Automatic curve fitting (Normal, Log-Normal, Beta, Gamma distributions)  
- Professional styling with customizable themes
- High-DPI PDF output optimized for academic publications
- Comprehensive statistical analysis and reporting
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from matplotlib.ticker import FuncFormatter, LogFormatter, ScalarFormatter
import argparse
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# Set high-quality defaults for academic publications
plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.family': 'serif',
    'font.serif': ['Computer Modern', 'Times', 'DejaVu Serif'],
    'text.usetex': False,  # Keep False for compatibility, use raw strings for math
    'mathtext.fontset': 'cm',  # Use Computer Modern for math
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.linewidth': 0.4,
    'axes.axisbelow': True,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.direction': 'out',
    'ytick.direction': 'out'
})


class MeasurementAnalyzer:
    """Professional analysis and visualization system for correlation or probability data"""
    
    def __init__(self, data, quantity='correlation', config=None):
        """Initialize analyzer with measurement data"""
        self.data = np.array(data)
        self.quantity = quantity
        self.config = config or {}
        self.quantity_title = 'Probability' if quantity == 'probability' else 'Correlation'
        self.quantity_label = 'Probability' if quantity == 'probability' else 'Correlation Magnitude'
        self.quantity_column = 'probability' if quantity == 'probability' else 'correlation'
        self.stats = self._compute_statistics()
        self.optimal_bins, self.bin_method = self._calculate_optimal_bins()
        self.bin_edges = self._generate_histogram_edges()
        
    def _compute_statistics(self):
        """Compute comprehensive statistical measures"""
        stats_dict = {
            'n': len(self.data),
            'mean': np.mean(self.data),
            'std': np.std(self.data, ddof=1),
            'median': np.median(self.data),
            'min': np.min(self.data),
            'max': np.max(self.data),
            'range': np.max(self.data) - np.min(self.data),
            'q25': np.percentile(self.data, 25),
            'q75': np.percentile(self.data, 75),
            'iqr': np.percentile(self.data, 75) - np.percentile(self.data, 25),
            'skewness': stats.skew(self.data),
            'kurtosis': stats.kurtosis(self.data),
            'cv': np.std(self.data, ddof=1) / np.abs(np.mean(self.data)) if np.mean(self.data) != 0 else np.inf
        }
        
        # Additional robust statistics
        stats_dict['mad'] = stats.median_abs_deviation(self.data)  # Median Absolute Deviation
        stats_dict['trimmed_mean'] = stats.trim_mean(self.data, 0.1)  # 10% trimmed mean
        
        return stats_dict
    
    def _calculate_optimal_bins(self):
        """Calculate optimal number of bins using multiple statistical rules."""
        n = max(1, self.stats['n'])

        if self.stats['iqr'] > 0:
            fd_bins = int(np.ceil(self.stats['range'] / (2 * self.stats['iqr'] * n ** (-1 / 3))))
        else:
            fd_bins = int(np.sqrt(n))

        sturges_bins = int(np.ceil(np.log2(n) + 1))

        if self.stats['std'] > 0:
            scott_bins = int(np.ceil(self.stats['range'] / (3.5 * self.stats['std'] * n ** (-1 / 3))))
        else:
            scott_bins = sturges_bins

        sqrt_bins = int(np.ceil(np.sqrt(n)))
        rice_bins = int(np.ceil(2 * n ** (1 / 3)))

        if abs(self.stats['skewness']) > 0.5 or abs(self.stats['kurtosis']) > 0.5:
            optimal = fd_bins
            method = "Freedman-Diaconis"
        elif n > 1000:
            optimal = rice_bins
            method = "Rice"
        else:
            optimal = scott_bins
            method = "Scott"

        # Adaptive minimum bins based on sample size - use statistical rules!
        if n < 50:
            default_min = max(10, n // 3)  # At least 3 samples per bin on average
        elif n < 200:
            # For 100 samples: use ~20-25 bins for good balance between detail and clarity
            default_min = max(15, int(np.sqrt(n)))  # Square root rule: sqrt(100) = 10, but use 15 minimum
        else:
            default_min = 30
        
        min_bins = int(self.config.get('min_bins', default_min))
        # Cap maximum bins at n/3 for good statistical representation (avg 3+ samples per bin)
        max_bins = int(self.config.get('max_bins', max(50, n // 3)))
        optimal = int(np.ceil(optimal))
        optimal = int(np.clip(optimal, min_bins, max_bins))
        method_annotation = method

        positive_data = self.data[self.data > 0]
        if self.quantity == 'correlation' and positive_data.size > 0:
            min_pos = np.min(positive_data)
            max_pos = np.max(positive_data)
            if min_pos > 0:
                dynamic_range = max_pos / max(min_pos, np.finfo(float).tiny)
                log2_span = np.log2(dynamic_range)
                
                # For narrow distributions, use many more bins to show detail
                if log2_span < 1.0:
                    # Very narrow: boost significantly
                    optimal = min(max_bins, max(optimal * 3, 256))
                elif log2_span < 2.0:
                    # Narrow: boost moderately
                    optimal = min(max_bins, max(optimal * 2, 128))
                elif dynamic_range >= 16:
                    optimal = min(max_bins, int(np.ceil(optimal * 1.5)))
                elif dynamic_range >= 64:
                    optimal = min(max_bins, int(np.ceil(optimal * 2.0)))

                # Align to power of two for log-scale axis
                if log2_span >= 2.0:
                    max_power = int(np.floor(np.log2(max_bins))) if max_bins > 0 else 10
                    min_power = int(np.ceil(np.log2(min_bins))) if min_bins > 0 else 7
                    power_target = int(np.clip(np.ceil(np.log2(optimal)), min_power, max_power))
                    optimal = 2 ** power_target
                    method_annotation = f"{method} (power-of-two aligned)"

        print("Bin selection analysis:")
        print(f"  Freedman-Diaconis: {fd_bins}")
        print(f"  Sturges: {sturges_bins}")
        print(f"  Scott: {scott_bins}")
        print(f"  Square-root: {sqrt_bins}")
        print(f"  Rice: {rice_bins}")
        print(f"  Selected: {optimal} ({method_annotation} rule)")

        return int(optimal), method_annotation

    def _power_of_two_range(self):
        """Return integer exponent range covering the positive data values."""
        positive = self.data[self.data > 0]
        if positive.size == 0:
            return None

        exp_min = int(np.floor(np.log2(np.min(positive))))
        exp_max = int(np.ceil(np.log2(np.max(positive))))

        if exp_min == exp_max:
            exp_min -= 1
            exp_max += 1

        return exp_min, exp_max

    def _generate_histogram_edges(self):
        """Generate histogram bin edges, aligning to powers of two when appropriate."""
        finite = self.data[np.isfinite(self.data)]
        if finite.size == 0:
            return np.linspace(0.0, 1.0, max(self.optimal_bins, 1) + 1)

        # For correlation data, use linear spacing in the actual data range
        # Don't extend to power-of-2 boundaries - that creates too much empty space
        min_val = np.min(finite)
        max_val = np.max(finite)
        
        if np.isclose(min_val, max_val):
            offset = abs(min_val) * 0.05 if min_val != 0 else 1e-6
            min_val -= offset
            max_val += offset

        return np.linspace(min_val, max_val, self.optimal_bins + 1)

    @staticmethod
    def _format_power_label(value, precision=2, include_decimal=False):
        """Format a positive value as a power-of-two string with optional decimal context."""
        if value <= 0 or not np.isfinite(value):
            return r"$0$" if not include_decimal else r"$0$" + "\n(0)"

        exponent = np.log2(value)
        if not np.isfinite(exponent):
            return r"$0$" if not include_decimal else r"$0$" + "\n(0)"

        # Use LaTeX formatting for clean rendering
        if abs(exponent - round(exponent)) < 0.05:
            exponent_text = rf"$2^{{{int(round(exponent))}}}$"
        else:
            exponent_text = rf"$2^{{{exponent:.{precision}f}}}$"

        if include_decimal:
            return f"{exponent_text}\n({value:.{precision + 1}e})"
        return exponent_text

    def _configure_correlation_axis(self, ax):
        """Configure correlation axis - use log scale only for wide dynamic ranges."""
        positive = self.data[self.data > 0]
        has_zero = np.any(self.data == 0)
        
        if positive.size == 0:
            return False

        # Calculate actual log2 span (not rounded)
        min_val = np.min(positive)
        max_val = np.max(positive)
        actual_log2_span = np.log2(max_val) - np.log2(min_val)
        
        # Only use log scale if the dynamic range is wide enough (> 2 orders of magnitude)
        # For narrow ranges, linear scale shows details better
        if actual_log2_span < 2.0:
            # Use linear scale for narrow distributions
            ax.set_xscale('linear')
            
            # Zoom to actual data range (min to max), not bin edges
            # Add moderate margin (5%) for a wider context view
            data_range = max_val - min_val
            margin = data_range * 0.05
            
            # If data includes zero, start from 0; otherwise from actual min
            x_min = 0 if has_zero else max(0, min_val - margin)
            x_max = max_val + margin
            ax.set_xlim(x_min, x_max)
            
            # Calculate bin centers from the histogram bin edges
            # bin_edges has n+1 values, bin centers are midpoints
            bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
            
            # Select key ticks: show more ticks for better resolution
            num_bins = len(bin_centers)
            tick_positions = []
            tick_labels = []
            
            # First bar center (or zero if present)
            if has_zero:
                tick_positions.append(0)
                tick_labels.append(r"$0$")
            
            # Determine how many ticks to show (aim for 12-15 ticks for better resolution)
            num_ticks = min(15, max(12, num_bins // 7))
            tick_indices = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
            
            for idx in tick_indices:
                tick_positions.append(bin_centers[idx])
                tick_labels.append(rf"$2^{{{np.log2(bin_centers[idx]):.2f}}}$")
            
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=45, ha='right')
            ax.set_xlabel(r'$|Cr|$', fontweight='bold')
            return True

        # Wide range: use log scale
        exponent_range = self._power_of_two_range()
        if exponent_range is None:
            return False
        
        exp_min, exp_max = exponent_range
        ax.set_xscale('log', base=2)
        # Zoom to actual data range with moderate margin (5%) to show wider context
        data_range = max_val - min_val
        margin_factor = 0.05
        lower_bound = min_val * (1 - margin_factor)
        upper_bound = max_val * (1 + margin_factor)
        ax.set_xlim(lower_bound, upper_bound)

        # Calculate bin centers from the histogram bin edges
        bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
        num_bins = len(bin_centers)
        
        # Select key bin centers for ticks (show more ticks for better resolution)
        major_ticks = []
        tick_labels = []
        
        # Add zero if any samples are exactly zero
        if has_zero:
            # Can't show 0 on log scale directly, so add annotation instead
            ax.annotate(r'$0$', xy=(lower_bound * 0.95, 0), xytext=(lower_bound * 0.8, ax.get_ylim()[1] * 0.1),
                       fontsize=10, color='red', weight='bold',
                       arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
        
        # Determine how many ticks to show (aim for 12-15 ticks for better resolution)
        num_ticks = min(15, max(12, num_bins // 7))
        tick_indices = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
        
        for idx in tick_indices:
            major_ticks.append(bin_centers[idx])
            log_val = np.log2(bin_centers[idx])
            tick_labels.append(rf"$2^{{{log_val:.2f}}}$")
        
        ax.set_xticks(major_ticks)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        ax.tick_params(axis='x', which='both', length=4, width=0.8)
        ax.set_xlabel(r'$|Cr|$', fontweight='bold')
        return True
    
    def fit_distributions(self):
        """Fit multiple probability distributions to the data"""
        distributions = {
            'Normal': stats.norm,
            'Log-Normal': stats.lognorm,
            'Gamma': stats.gamma,
            'Beta': stats.beta,
            'Exponential': stats.expon,
            'Weibull': stats.weibull_min
        }
        
        results = {}
        
        # Only fit positive distributions to positive data
        positive_data = self.data[self.data > 0] if np.any(self.data > 0) else self.data
        
        for name, dist in distributions.items():
            try:
                if name in ['Log-Normal', 'Gamma', 'Exponential', 'Weibull'] and np.any(self.data <= 0):
                    # Skip distributions that require positive data if we have negative values
                    continue
                    
                if name == 'Beta' and (np.min(self.data) < 0 or np.max(self.data) > 1):
                    # Beta distribution requires data in [0,1]
                    continue
                    
                # Fit distribution
                if name == 'Log-Normal':
                    params = dist.fit(positive_data)
                else:
                    params = dist.fit(self.data)
                    
                # Calculate AIC and BIC
                loglik = np.sum(dist.logpdf(positive_data if name == 'Log-Normal' else self.data, *params))
                k = len(params)
                n = len(positive_data if name == 'Log-Normal' else self.data)
                
                aic = 2*k - 2*loglik
                bic = k*np.log(n) - 2*loglik
                
                # Kolmogorov-Smirnov test
                ks_stat, ks_p = stats.kstest(positive_data if name == 'Log-Normal' else self.data, 
                                           lambda x: dist.cdf(x, *params))
                
                results[name] = {
                    'params': params,
                    'aic': aic,
                    'bic': bic,
                    'ks_statistic': ks_stat,
                    'ks_pvalue': ks_p,
                    'loglik': loglik
                }
                
            except Exception as e:
                print(f"Warning: Could not fit {name} distribution: {e}")
                continue
        
        if results:
            # Find best distribution by AIC
            best_dist = min(results.keys(), key=lambda x: results[x]['aic'])
            print(f"\nDistribution fitting results:")
            print(f"Best fit: {best_dist} (lowest AIC: {results[best_dist]['aic']:.2f})")
            
            for name, result in sorted(results.items(), key=lambda x: x[1]['aic']):
                print(f"  {name:12s}: AIC={result['aic']:8.2f}, BIC={result['bic']:8.2f}, KS p-value={result['ks_pvalue']:.4f}")
        
        return results
    def create_visualization(self, output_path=None, theme='academic'):
        """Create professional visualization for the configured measurement"""
        # Set theme
        if theme == 'academic':
            colors = {
                'primary': '#1f77b4',    # Professional blue
                'secondary': '#ff7f0e',   # Orange
                'accent': '#d62728',      # Red  
                'fit_curve': '#2ca02c',   # Green
                'stats': '#9467bd',       # Purple
                'grid': '#e0e0e0'
            }
        elif theme == 'nature':
            colors = {
                'primary': '#2E8B57',     # Sea green
                'secondary': '#FF6347',   # Tomato
                'accent': '#4169E1',      # Royal blue
                'fit_curve': '#32CD32',   # Lime green
                'stats': '#8A2BE2',       # Blue violet
                'grid': '#f0f0f0'
            }
        else:
            colors = {
                'primary': '#1f77b4',
                'secondary': '#ff7f0e',
                'accent': '#d62728',
                'fit_curve': '#2ca02c',
                'stats': '#9467bd',
                'grid': '#e0e0e0'
            }

        self.bin_edges = self._generate_histogram_edges()
        
        # Create figure with single panel (no statistics sidebar)
        fig, ax_main = plt.subplots(1, 1, figsize=(10, 6))
        
        # Main histogram plot with NO GAPS between bars
        counts, bins, patches = ax_main.hist(
            self.data,
            bins=self.bin_edges,
            alpha=0.82,
            color=colors['primary'],
            edgecolor='#0f172a',
            linewidth=0.4,
            density=True,
            rwidth=1.0,  # Full width bars - NO GAPS!
            label=f"Data (n={self.stats['n']:,})"
        )
        ax_main.set_facecolor('#F8FAFF')
        
        # Use Kernel Density Estimation for accurate non-parametric fit
        from scipy.stats import gaussian_kde
        
        try:
            # Create KDE with automatically selected bandwidth
            kde = gaussian_kde(self.data, bw_method='scott')
            
            # Generate smooth curve
            data_range = self.stats['max'] - self.stats['min']
            x_min = max(0, self.stats['min'] - 0.05 * data_range)
            x_max = self.stats['max'] + 0.05 * data_range
            x_smooth = np.linspace(x_min, x_max, 5000)
            y_smooth = kde(x_smooth)
            
            ax_main.plot(x_smooth, y_smooth, color=colors['fit_curve'], 
                       linewidth=2.5, label='KDE fit', alpha=0.9, zorder=5)
        except Exception as e:
            print(f"Warning: KDE fitting failed ({e}), trying parametric distributions...")
            
            # Fallback to parametric distributions if KDE fails
            fit_results = self.fit_distributions()
            if fit_results:
                best_dist_name = min(fit_results.keys(), key=lambda x: fit_results[x]['aic'])
                best_result = fit_results[best_dist_name]
                
                # Get distribution object
                dist_map = {
                    'Normal': stats.norm,
                    'Log-Normal': stats.lognorm, 
                    'Gamma': stats.gamma,
                    'Beta': stats.beta,
                    'Exponential': stats.expon,
                    'Weibull': stats.weibull_min
                }
                
                if best_dist_name in dist_map:
                    dist = dist_map[best_dist_name]
                    params = best_result['params']
                    
                    # Use very fine resolution for smooth curve (5000 points)
                    # Extend slightly beyond data range for better visualization
                    data_range = self.stats['max'] - self.stats['min']
                    x_min = max(0, self.stats['min'] - 0.02 * data_range)
                    x_max = self.stats['max'] + 0.02 * data_range
                    
                    x_smooth = np.linspace(x_min, x_max, 5000)
                    y_smooth = dist.pdf(x_smooth, *params)
                    
                    # Remove any NaN or Inf values that might occur at boundaries
                    valid_mask = np.isfinite(y_smooth)
                    x_smooth = x_smooth[valid_mask]
                    y_smooth = y_smooth[valid_mask]
                    
                    ax_main.plot(x_smooth, y_smooth, color=colors['fit_curve'], 
                               linewidth=2.5, label=f'{best_dist_name} fit', alpha=0.9, zorder=5)
        
        # Add statistical reference lines with power-of-2 notation
        mean_label = self._format_power_label(abs(self.stats['mean'])) if self.stats['mean'] != 0 else r'$0$'
        median_label = self._format_power_label(abs(self.stats['median'])) if self.stats['median'] != 0 else r'$0$'
        std_label = self._format_power_label(self.stats['std']) if self.stats['std'] > 0 else r'$0$'
        
        ax_main.axvline(self.stats['mean'], color=colors['accent'], 
                       linestyle='--', linewidth=2, alpha=0.8, 
                       label=f'Mean: {mean_label}, $\\sigma$: {std_label}')
        ax_main.axvline(self.stats['median'], color=colors['secondary'], 
                       linestyle=':', linewidth=2, alpha=0.8, 
                       label=f'Median: {median_label}')
        
        # Add confidence intervals (95% CI)
        ci_lower = self.stats['mean'] - 1.96 * self.stats['std']
        ci_upper = self.stats['mean'] + 1.96 * self.stats['std']

        if self.quantity == 'correlation' and np.any(self.data > 0):
            min_positive = np.min(self.data[self.data > 0])
            ci_lower = max(ci_lower, min_positive)

        ci_lower = max(ci_lower, np.min(self.bin_edges))
        ci_upper = min(ci_upper, np.max(self.bin_edges))

        if ci_lower < ci_upper and np.isfinite(ci_lower) and np.isfinite(ci_upper):
            ax_main.axvspan(ci_lower, ci_upper, alpha=0.1, color=colors['stats'],
                           label=r'95% CI ($\pm 1.96\sigma$)')
        
        # Format main plot labels
        ax_main.set_ylabel('Probability Density', fontweight='bold')
        ax_main.grid(which='both', alpha=0.25)

        if self.quantity == 'probability':
            ax_main.set_xlabel('Probability', fontweight='bold')

            def probability_formatter(x, pos):
                if x < 0 or x > 1:
                    return f"{x:.4f}"
                # Ensure small probabilities are not rounded to zero
                if x < 1e-3:
                    return f"{x:.4e}"
                return f"{x:.4f}"

            ax_main.xaxis.set_major_formatter(FuncFormatter(probability_formatter))
        else:
            axis_configured = self._configure_correlation_axis(ax_main)

            if not axis_configured:
                ax_main.set_xlabel(self.quantity_label, fontweight='bold')

                def linear_log2_formatter(x, pos):
                    if x <= 0:
                        return f"{x:.1e}" if x < 0 else "0"
                    return self._format_power_label(x, include_decimal=True)

                ax_main.xaxis.set_major_formatter(FuncFormatter(linear_log2_formatter))
        
        # Legend with better positioning
        ax_main.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, framealpha=0.9)
        
        plt.tight_layout()
        
        # Save with high quality
        if output_path:
            # Save PDF for publications
            plt.savefig(output_path.replace('.png', '.pdf'), format='pdf', 
                       bbox_inches='tight', dpi=300)
            # Also save PNG for presentations
            plt.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
            print(f"Visualization saved to: {output_path}")
            print(f"PDF version saved to: {output_path.replace('.png', '.pdf')}")
        
        return fig
    
    def export_data_summary(self, output_path):
        """Export comprehensive data analysis to CSV and text files"""
        # CSV export for further analysis
        edges = self.bin_edges if isinstance(self.bin_edges, np.ndarray) and len(self.bin_edges) > 1 else self._generate_histogram_edges()
        csv_data = {
            self.quantity_column: self.data,
            'bin_assignment': np.digitize(self.data, edges, right=False)
        }
        
        df = pd.DataFrame(csv_data)
        df.to_csv(output_path.replace('.txt', '.csv'), index=False)
        
        # Detailed text report with power-of-2 representation
        with open(output_path, 'w') as f:
            f.write(f"Cryptographic Cipher {self.quantity_title} Analysis Report\n")
            f.write("="*50 + "\n\n")
            
            # Power-of-2 calculations
            mean_log2 = np.log2(abs(self.stats['mean'])) if self.stats['mean'] != 0 else -np.inf
            std_log2 = np.log2(self.stats['std']) if self.stats['std'] > 0 else -np.inf
            median_log2 = np.log2(abs(self.stats['median'])) if self.stats['median'] != 0 else -np.inf
            
            f.write("Dataset Information:\n")
            f.write(f"  Total samples: {self.stats['n']:,}\n")
            f.write(f"  Data range: [{self.stats['min']:.6e}, {self.stats['max']:.6e}]\n")
            if self.stats['min'] != 0:
                f.write(f"  Dynamic range: {self.stats['max']/self.stats['min']:.2e}\n")
            else:
                f.write("  Dynamic range: undefined (minimum is zero)\n")
            min_log = np.log2(abs(self.stats['min'])) if self.stats['min'] not in (0, -0.0) else None
            max_log = np.log2(abs(self.stats['max'])) if self.stats['max'] != 0 else None
            min_log_str = f"{min_log:.2f}" if (min_log is not None and np.isfinite(min_log)) else "-∞"
            max_log_str = f"{max_log:.2f}" if (max_log is not None and np.isfinite(max_log)) else "-∞"
            f.write(f"  Log₂ range: [{min_log_str}, {max_log_str}]\n\n")
            
            f.write("Central Tendency (Power-of-2 Representation):\n")
            f.write(f"  Mean: 2^{mean_log2:.4f} ({self.stats['mean']:.6e})\n")
            f.write(f"  Median: 2^{median_log2:.4f} ({self.stats['median']:.6e})\n") 
            f.write(f"  Trimmed mean (10%): {self.stats['trimmed_mean']:.6e}\n\n")
            
            f.write("Dispersion:\n")
            f.write(f"  Standard deviation: 2^{std_log2:.4f} ({self.stats['std']:.6e})\n")
            f.write(f"  Median abs. deviation: {self.stats['mad']:.6e}\n")
            f.write(f"  Interquartile range: {self.stats['iqr']:.6e}\n")
            f.write(f"  Coefficient of variation: {self.stats['cv']:.4f}\n\n")
            
            f.write("Interpretation:\n")
            if self.quantity == 'probability':
                f.write(f"  Mean probability (log₂): {mean_log2:.2f} bits\n")
                f.write(f"  Peak probability: {self.stats['max']:.6e}\n")
                f.write(f"  Coverage above 50%: {np.mean(self.data >= 0.5)*100:.2f}%\n\n")
            else:
                f.write(f"  Correlation strength (log₂): {mean_log2:.2f} bits\n")
                f.write(f"  Security margin: ~{-mean_log2/2:.1f} bits\n")
                f.write(f"  Detection probability: 2^{2*mean_log2:.2f}\n\n")
            
            f.write("Shape Statistics:\n")
            f.write(f"  Skewness: {self.stats['skewness']:.4f}\n")
            f.write(f"  Kurtosis: {self.stats['kurtosis']:.4f}\n\n")
            
            f.write("Histogram Configuration:\n")
            f.write(f"  Optimal bins: {self.optimal_bins} ({self.bin_method})\n")
            if self.quantity == 'correlation' and self._power_of_two_range() is not None:
                positive_edges = edges[edges > 0]
                if len(positive_edges) > 1:
                    bin_step = np.median(np.diff(np.log2(positive_edges)))
                    f.write(f"  Bin step (log₂): {bin_step:.3f}\n")
                else:
                    f.write("  Bin step (log₂): n/a\n")
            else:
                if len(edges) > 1:
                    bin_width = np.median(np.diff(edges))
                    f.write(f"  Median bin width: {bin_width:.6e}\n")
                else:
                    f.write("  Bin width: n/a\n")
            f.write("\n")
            
            # Distribution fitting results
            fit_results = self.fit_distributions()
            if fit_results:
                f.write("Distribution Fitting Results:\n")
                for name, result in sorted(fit_results.items(), key=lambda x: x[1]['aic']):
                    f.write(f"  {name}:\n")
                    f.write(f"    AIC: {result['aic']:.2f}\n")
                    f.write(f"    BIC: {result['bic']:.2f}\n")
                    f.write(f"    KS p-value: {result['ks_pvalue']:.6f}\n")
                    f.write(f"    Parameters: {result['params']}\n\n")


def parse_orthros_csv(file_path, header_lines):
    """Parse Orthros-specific CSV format with metadata extraction"""
    metadata = {}
    
    # Extract metadata from header comments
    for line in header_lines:
        if line.startswith('# Mode'):
            metadata['mode'] = int(line.split(':')[1].strip())
        elif line.startswith('# Offset'):
            metadata['offset'] = int(line.split(':')[1].strip())
        elif line.startswith('# ΔL'):
            metadata['delta_left'] = line.split(':')[1].strip()
        elif line.startswith('# ΔR'):
            metadata['delta_right'] = line.split(':')[1].strip()
        elif line.startswith('# Output'):
            metadata['output_mask'] = line.split(':')[1].strip()
    
    # Try multiple parsing strategies
    correlations = None
    
    # Strategy 1: Try normal CSV with comma separation and headers
    try:
        df = pd.read_csv(file_path, comment='#')
        if 'correlation' in df.columns:
            correlations = df['correlation'].values
        elif 'Correlation' in df.columns:
            correlations = df['Correlation'].values
        else:
            # Look for any column with correlation-like name
            corr_cols = [col for col in df.columns if 'corr' in col.lower()]
            if corr_cols:
                correlations = df[corr_cols[0]].values
            else:
                correlations = df.iloc[:, -1].values  # Last column
    except:
        pass
    
    # Strategy 2: Try space-separated format (old format)
    if correlations is None:
        try:
            df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None, 
                            names=['ExperimentID', 'Key', 'Correlation'])
            correlations = df['Correlation'].values
        except:
            pass
    
    # Strategy 3: Load as raw numeric data
    if correlations is None:
        try:
            data = np.loadtxt(file_path, comments='#')
            if data.ndim == 1:
                correlations = data
            else:
                correlations = data[:, -1]  # Last column should be correlation
        except:
            pass
    
    if correlations is None:
        raise ValueError("Could not parse correlation data from Orthros CSV file")
    
    # Store metadata for potential use in visualization
    if hasattr(correlations, '__dict__'):
        correlations.metadata = metadata
    
    return correlations


def load_measurement_data(file_path):
    """Load correlation or probability data from various file formats"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Data file not found: {file_path}")
    
    file_ext = os.path.splitext(file_path)[1].lower()
    
    try:
        if file_ext == '.csv':
            # First, check if this is an Orthros file by reading the first few lines
            with open(file_path, 'r') as f:
                header_lines = [f.readline().strip() for _ in range(10)]
            
            is_orthros = any('Orthros' in line for line in header_lines)
            
            if is_orthros:
                # Parse Orthros-specific format
                return parse_orthros_csv(file_path, header_lines)
            else:
                # Standard CSV format
                df = pd.read_csv(file_path, comment='#')
                # Try to automatically detect measurement column
                measurement_cols = [
                    col for col in df.columns
                    if any(keyword in col.lower() for keyword in ['corr', 'prob', 'value', 'metric'])
                ]
                if measurement_cols:
                    return df[measurement_cols[0]].values
                else:
                    return df.iloc[:, -1].values  # Last column
                
        elif file_ext in ['.txt', '.dat']:
            # Try structured format first (with headers)
            try:
                df = pd.read_csv(file_path, sep=r'\s+', comment='#')
                measurement_cols = [
                    col for col in df.columns
                    if any(keyword in col.lower() for keyword in ['corr', 'prob', 'value', 'metric'])
                ]
                if measurement_cols:
                    return df[measurement_cols[0]].values
                else:
                    return df.iloc[:, -1].values  # Last column
            except:
                # Fallback: load as simple numeric data
                data = np.loadtxt(file_path, comments='#')
                if data.ndim == 1:
                    return data
                else:
                    return data[:, -1]  # Last column
                
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path)
            measurement_cols = [
                col for col in df.columns
                if any(keyword in col.lower() for keyword in ['corr', 'prob', 'value', 'metric'])
            ]
            if measurement_cols:
                return df[measurement_cols[0]].values
            else:
                return df.iloc[:, -1].values  # Last column
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
            
    except Exception as e:
        raise ValueError(f"Error loading data from {file_path}: {e}")


def infer_quantity(data, mode='auto'):
    """Infer whether the dataset represents correlation or probability values"""
    if mode in ('correlation', 'probability'):
        return mode

    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return 'correlation'

    eps = 1e-9
    min_val = np.min(finite)
    max_val = np.max(finite)

    # If values are in [0,1] range, check if they're small (likely correlations)
    # or near 0.5-1.0 (likely probabilities)
    if min_val >= -eps and max_val <= 1.0 + eps:
        if np.any(finite < -eps):
            return 'correlation'
        # Small positive values (< 0.01) are likely correlations
        if max_val < 0.01:
            return 'correlation'
        return 'probability'

    return 'correlation'


def main():
    parser = argparse.ArgumentParser(
        description='Cryptographic Cipher Visualization Toolkit (Correlation / Probability)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualizer.py data.csv
  python visualizer.py cipher_8r_keys_*.txt --theme nature
  python visualizer.py probabilities.csv --quantity probability --output results/probability.png
        """
    )
    
    parser.add_argument('input_file', help='Input data file (CSV, TXT, Excel)')
    parser.add_argument('--output', '-o', default=None, 
                       help='Output file path for visualization (default: auto-generated)')
    parser.add_argument('--theme', choices=['academic', 'nature'], default='academic',
                       help='Visual theme for plots')
    parser.add_argument('--quantity', choices=['auto', 'correlation', 'probability'], default='auto',
                       help='Type of values to visualize; auto attempts detection based on the dataset')
    parser.add_argument('--export-data', action='store_true',
                       help='Export detailed analysis data')
    parser.add_argument('--no-display', action='store_true',
                       help='Do not display plot (useful for batch processing)')
    
    args = parser.parse_args()
    
    try:
        # Load data
        print(f"Loading data from: {args.input_file}")
        measurement_data = load_measurement_data(args.input_file)
        
        # Check for and handle NaN values
        initial_count = len(measurement_data)
        nan_count = np.sum(np.isnan(measurement_data))
        if nan_count > 0:
            print(f"Warning: Found {nan_count} NaN values, removing them from analysis")
            measurement_data = measurement_data[~np.isnan(measurement_data)]
        
        final_count = len(measurement_data)
        print(f"Loaded {final_count:,} valid samples (removed {initial_count - final_count} invalid entries)")
        
        # Verify we have data
        if final_count == 0:
            raise ValueError("No valid data found")

        # Determine quantity type
        quantity = infer_quantity(measurement_data, args.quantity)
        print(f"Detected quantity: {quantity}")
        
        # Show data range for verification
        min_val = np.min(measurement_data)
        max_val = np.max(measurement_data)
        print(f"Data range: [{min_val:.6e}, {max_val:.6e}]")
        zero_count = np.sum(measurement_data == 0.0)
        negative_count = np.sum(measurement_data < 0.0)
        positive_count = np.sum(measurement_data > 0.0)
        print(f"Sample distribution: {positive_count} positive, {zero_count} zero, {negative_count} negative")
        
        # Analyze data
        analyzer = MeasurementAnalyzer(measurement_data, quantity=quantity)
        
        # Generate output filename if not provided
        if args.output is None:
            base_name = os.path.splitext(os.path.basename(args.input_file))[0]
            args.output = f"{base_name}_{quantity}_analysis.png"
        
        # Create visualization
        print("\nCreating visualization...")
        fig = analyzer.create_visualization(args.output, theme=args.theme)
        
        # Export detailed analysis if requested
        if args.export_data:
            report_path = args.output.replace('.png', '_report.txt')
            analyzer.export_data_summary(report_path)
            print(f"Detailed analysis exported to: {report_path}")
        
        # Display plot unless disabled
        if not args.no_display:
            plt.show()
            
        print("\nAnalysis completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
