#!/usr/bin/env python3
"""
Copyright (C) 2025 Hosein Hadipour and Mostafizar Rahman
Email: hsn.hadipour@gmail.com
Date: September 30, 2025

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

Dynamic Statistical Threshold Calculator for Differential Cryptanalysis
Revised model: CORR = |c0 - c1| is modeled as a half-normal variable under H0.

Previous version incorrectly centered thresholds at n/2 using a binomial
count interpretation. We now treat D = c0 - c1 where c0 + c1 = n and
under the null D ~ approx Normal(0, n). Hence CORR = |D| follows a
HalfNormal(sigma=√n) distribution.

Key formulas:
    E[CORR] = √n * √(2/π)
    Var[CORR] = n * (1 - 2/π)
    Tail: P(CORR ≥ t) = 2 * (1 - Φ(t/√n))
Multiple testing (Bonferroni over T tests): choose t such that
    2T (1 - Φ(t/√n)) ≤ α  ⇒  t = √n * Φ^{-1}(1 - α/(2T))

Parameter usage:
    CORR : Expected correlation parameter (approx 2ε). We set ε = CORR / 2.
    (Removed EPSILON / ALPHA overrides for simplicity; using fixed global α=0.05.)
"""

import math
import re
import os
from pathlib import Path
import argparse

try:
    from statistics import NormalDist as _NormalDist
    STATISTICS_AVAILABLE = True
except Exception:
    STATISTICS_AVAILABLE = False

class DynamicCorrelationThresholdCalculator:
    def __init__(self, makefile_path="makefile", data_dir=None):
        """
        Initialize dynamic threshold calculator.
        
        Args:
            makefile_path: Path to makefile containing attack parameters
            data_dir: Directory containing attack output data (optional)
        """
        self.makefile_path = makefile_path
        self.data_dir = data_dir
        
        # Extract parameters from makefile
        self.params = self.extract_attack_parameters()

        # Core sample size
        self.n = 2 ** self.params['DEG']  # Number of pairs per test
        # Half-normal parameters for CORR = |c0 - c1|
        self.sqrt_n = math.sqrt(self.n)
        self.expected_mean = self.sqrt_n * math.sqrt(2 / math.pi)
        self.expected_std = self.sqrt_n * math.sqrt(1 - 2 / math.pi)

        # Convert CORR from -log2 format to actual correlation value
        # CORR in makefile is -log2(correlation), e.g., 11 means 2^{-11}
        if 'CORR' in self.params:
            self.corr_log2 = self.params['CORR']  # Store original -log2 value
            self.corr_value = 2 ** (-self.corr_log2)  # Actual correlation value
        else:
            self.corr_log2 = None
            self.corr_value = 0.0

        # Estimate number of tests from data if available
        self.num_tests = self.estimate_num_tests()
        
    def extract_attack_parameters(self):
        """Extract attack parameters from makefile."""
        params = {}
        
        try:
            with open(self.makefile_path, 'r') as f:
                content = f.read()
            
            # Extract key parameters - use more specific patterns
            patterns = {
                'DEG': r'^\s*DEG\s*[=:]\s*(\d+)',
                'M': r'^\s*M\s*[=:]\s*(\d+(?:\.\d+)?)',
                'CORR': r'^\s*CORR\s*[=:]\s*(\d+(?:\.\d+)?)',
                'C0': r'^\s*C0\s*[=:]\s*(\d+(?:\.\d+)?)',
                'EPSILON': r'^\s*EPSILON\s*[=:]\s*(\d+(?:\.\d+)?)',
                'ALPHA': r'^\s*ALPHA\s*[=:]\s*(\d+(?:\.\d+)?)',
                'N': r'^\s*N\s*[=:]\s*(\d+)',
                'NUM_OF_ROUNDS': r'^\s*NUM_OF_ROUNDS\s*[=:]\s*(\d+)',
                'RND_OFFSET': r'^\s*RND_OFFSET\s*[=:]\s*(\d+)',
                'NAME': r'^\s*NAME\s*[=:]\s*"([^"]+)"'
            }
            
            for param, pattern in patterns.items():
                match = re.search(pattern, content, re.MULTILINE)  # Use MULTILINE flag
                if match:
                    value = match.group(1)
                    # Convert to appropriate type
                    if param in ['DEG', 'N', 'NUM_OF_ROUNDS', 'RND_OFFSET']:
                        params[param] = int(value)
                    elif param == 'NAME':
                        params[param] = value
                    else:
                        params[param] = float(value)

            dy_patterns = {
                'DY1': r'^\s*dy1\s*=\s*((?:0x[0-9a-fA-F]+\s+)*0x[0-9a-fA-F]+)',
                'DY2': r'^\s*dy2\s*=\s*((?:0x[0-9a-fA-F]+\s+)*0x[0-9a-fA-F]+)'
            }

            for dy_label, pattern in dy_patterns.items():
                match = re.search(pattern, content, re.MULTILINE)
                if match:
                    tokens = match.group(1).split()
                    nibs = [int(tok, 16) & 0xF for tok in tokens]
                    params[f'{dy_label}_NIBS'] = nibs
                    if len(nibs) >= 2:
                        params[f'{dy_label}_BYTE'] = ((nibs[0] & 0xF) << 4) | (nibs[1] & 0xF)
                        
        except FileNotFoundError:
            print(f"Warning: Makefile not found at {self.makefile_path}")
            # Use default values
            params = {
                'DEG': 10, 'M': 5.0, 'CORR': 0.203, 'C0': 3.0,
                'EPSILON': 0.01, 'ALPHA': 0.5, 'N': 1,
                'NUM_OF_ROUNDS': 5, 'RND_OFFSET': 2, 'NAME': 'unknown'
            }
        
        return params
    
    def estimate_num_tests(self):
        """Estimate number of tests from data directory or makefile."""
        
        # Try to count from actual data files
        if self.data_dir:
            data_path = Path(self.data_dir)
            if data_path.exists():
                output_files = list(data_path.glob("out_*.txt"))
                if output_files:
                    # Count lines in first file
                    try:
                        with open(output_files[0], 'r') as f:
                            return sum(1 for line in f if line.strip())
                    except:
                        pass
        
        # Estimate from parameters
        m_value = self.params.get('M', 5.0)
        estimated_tests = int(m_value ** 2 * 200)  # Rough estimation
        
        return max(estimated_tests, 1000)  # Minimum reasonable number
    
    def norm_ppf(self, p):
        """Normal inverse CDF using Python's standard library (statistics.NormalDist)."""
        if not STATISTICS_AVAILABLE:
            raise RuntimeError("statistics.NormalDist not available; requires Python 3.8+.")
        return _NormalDist().inv_cdf(p)
    
    def norm_cdf(self, x):
        """Normal CDF using Python's standard library (statistics.NormalDist)."""
        if not STATISTICS_AVAILABLE:
            raise RuntimeError("statistics.NormalDist not available; requires Python 3.8+.")
        return _NormalDist().cdf(x)
    
    def compute_thresholds(self):
        """Compute thresholds in half-normal model.

        Thresholds are expressed as t = z * √n (no centering) because mean≈0 for |D| before absolute.
        We supply liberal (z=2), moderate (z=2.5), conservative (z=3) plus
        classical significance levels translated to half-normal tails and Bonferroni.
        """

        thresholds = {}

        # Core z-levels (practical cryptanalytic)
        practical_levels = [
            (2.0, 'liberal', 'z=2 (≈4.55% two-tailed; half-normal tail ≈4.55%)'),
            (2.5, 'moderate', 'z=2.5 (≈1.24% tail)'),
            (3.0, 'conservative', 'z=3 (≈0.27% tail)')
        ]
        for z, name, desc in practical_levels:
            thresholds[name] = {
                'threshold': z * self.sqrt_n,
                'z_score': z,
                'description': desc,
                'type': 'practical'
            }

        # Classical significance levels (convert α to threshold using half-normal tail: P=2(1-Φ(z)))
        alpha_levels = [0.10, 0.05, 0.01, 0.001]
        for alpha in alpha_levels:
            # solve 2(1-Φ(z)) = alpha => Φ(z) = 1 - alpha/2
            z_critical = self.norm_ppf(1 - alpha / 2)
            thresholds[f'alpha_{alpha:.3f}'] = {
                'threshold': z_critical * self.sqrt_n,
                'z_score': z_critical,
                'description': f'Half-normal tail α={alpha}',
                'type': 'statistical'
            }

        # Multiple testing Bonferroni with fixed global α=0.05 (removed Makefile override)
        global_alpha = 0.05
        per_test_alpha = global_alpha / self.num_tests
        # Solve 2(1-Φ(z)) = per_test_alpha => z = Φ^{-1}(1 - per_test_alpha/2)
        z_bonf = self.norm_ppf(1 - per_test_alpha / 2)
        thresholds['bonferroni'] = {
            'threshold': z_bonf * self.sqrt_n,
            'z_score': z_bonf,
            'description': f'Bonferroni α=0.05 over {self.num_tests} tests',
            'type': 'multiple_testing'
        }

        # Information-theoretic (expected bias based):
        # CORR in the makefile is provided as -log2(correlation). We parsed it as
        #   self.corr_log2 (original) and self.corr_value = 2^{-self.corr_log2} (actual correlation).
        # Expected correlation corr ≈ 2ε ⇒ ε ≈ corr/2. Bias in D mean ≈ ε n, so E[|D|] under H1 ≈ |ε| n.
        corr_param = self.corr_value if self.corr_log2 is not None else 0.0
        epsilon = corr_param / 2.0
        expected_shift = abs(epsilon) * self.n
        thresholds['information_theoretic'] = {
            'threshold': expected_shift,  # This can be very large relative to √n
            'z_score': expected_shift / self.sqrt_n if self.sqrt_n else 0.0,
            'description': f'Expected |D| under H1 (ε={epsilon:.3f})',
            'type': 'information_theoretic'
        }

        return thresholds
    
    def analyze_correlation_strength(self, correlation_value):
        """Analyze statistical significance using half-normal model."""
        z_score = correlation_value / self.sqrt_n if self.sqrt_n else 0.0
        # Half-normal tail P(|Z| >= z) = 2(1-Φ(z))
        p_value = 2 * (1 - self.norm_cdf(z_score))

        if z_score >= 5:
            sig = "EXTREME"
        elif z_score >= 4:
            sig = "VERY STRONG"
        elif z_score >= 3:
            sig = "STRONG"
        elif z_score >= 2.5:
            sig = "MODERATE"
        elif z_score >= 2.0:
            sig = "WEAK"
        else:
            sig = "NOT SIGNIFICANT"

        return {
            'z_score': z_score,
            'p_value_half_normal': p_value,
            'significance': sig
        }
    
    def recommend_threshold(self, context='general'):
        """Recommend threshold key based on usage context (half-normal)."""
        thresholds = self.compute_thresholds()
        mapping = {
            'general': 'conservative',
            'screening': 'moderate',
            'exploratory': 'liberal',
            'production': 'bonferroni',
            'research': 'alpha_0.010'
        }
        key = mapping.get(context, 'conservative')
        if key not in thresholds:
            key = 'conservative'
        return key, thresholds[key]
    
    def print_attack_parameters(self):
        """Print the extracted attack parameters and half-normal stats."""
        print("EXTRACTED ATTACK PARAMETERS (Half-Normal Model):")
        print("-" * 70)
        for key, value in self.params.items():
            print(f"  {key:<18}: {value}")
        print(f"  {'DATA_COMPLEXITY':<18}: 2^{self.params['DEG']} = {self.n}")
        print(f"  {'sqrt(n)':<18}: {self.sqrt_n:.2f}")
        print(f"  {'E[CORR]':<18}: {self.expected_mean:.2f}")
        print(f"  {'SD[CORR]':<18}: {self.expected_std:.2f}")
        print(f"  {'NUM_TESTS':<18}: {self.num_tests}")
        print()
    
    def print_threshold_table(self):
        """Print a comprehensive table of all thresholds."""
        
        thresholds = self.compute_thresholds()
        
        print("=" * 100)
        print("DYNAMIC STATISTICAL THRESHOLDS FOR DIFFERENTIAL CRYPTANALYSIS")
        print("=" * 100)
        
        # Print attack parameters first
        self.print_attack_parameters()
        
        print("-" * 100)
        print(f"{'Threshold Name':<25} {'Value':<12} {'Z':<8} {'Description':<45}")
        print("-" * 100)
        
        # Sort thresholds by value
        sorted_thresholds = sorted(thresholds.items(), 
                                 key=lambda x: x[1]['threshold'])
        
        for name, data in sorted_thresholds:
            print(f"{name:<25} {data['threshold']:<12.2f} {data['z_score']:<8.2f} {data['description']:<45}")
        
        print("-" * 100)
        
        # Recommendations
        print("\nRECOMMENDATIONS FOR DIFFERENT CONTEXTS (threshold = z * √n):")
        contexts = ['general', 'screening', 'exploratory', 'production', 'research']
        for context in contexts:
            key, rec = self.recommend_threshold(context)
            print(f"  {context.upper():<12}: {rec['threshold']:.2f} (z={rec['z_score']:.2f}, {rec['description']})")
        
        print("=" * 100)
    
    def get_visualization_thresholds(self):
        """Get the main thresholds for visualization."""
        _, conservative = self.recommend_threshold('general')
        _, moderate = self.recommend_threshold('screening')
        _, liberal = self.recommend_threshold('exploratory')
        return {
            'conservative': conservative['threshold'],
            'moderate': moderate['threshold'],
            'liberal': liberal['threshold'],
            'random': self.expected_mean  # Baseline mean of half-normal
        }

def main():
    """Main function to compute and display thresholds dynamically."""
    
    parser = argparse.ArgumentParser(description='Dynamic threshold calculator for differential cryptanalysis')
    parser.add_argument('--makefile', '-m', default='makefile', help='Path to makefile')
    parser.add_argument('--data-dir', '-d', help='Data directory for empirical analysis')
    parser.add_argument('--analyze-value', '-v', type=int, help='Analyze specific correlation value')
    
    args = parser.parse_args()
    
    # Initialize calculator
    calc = DynamicCorrelationThresholdCalculator(
        makefile_path=args.makefile,
        data_dir=args.data_dir
    )
    
    # Print comprehensive threshold table
    calc.print_threshold_table()
    
    # Example analysis of correlation values
    if args.analyze_value:
        print(f"\nANALYSIS OF CORRELATION VALUE: {args.analyze_value}")
        print("-" * 50)
        analysis = calc.analyze_correlation_strength(args.analyze_value)
        for key, value in analysis.items():
            print(f"  {key}: {value}")
    else:
        print("EXAMPLE CORRELATION ANALYSIS (using half-normal):")
        print("-" * 50)
        # Provide example z multiples
        example_z = [1, 2, 2.5, 3, 4, 5]
        for z in example_z:
            value = z * calc.sqrt_n
            analysis = calc.analyze_correlation_strength(value)
            print(f"Value={value:7.2f} (z={z:3.1f}) -> P={analysis['p_value_half_normal']:.2e}, {analysis['significance']}")
    
    # Get visualization thresholds
    viz_thresholds = calc.get_visualization_thresholds()
    print(f"\nVISUALIZATION THRESHOLDS:")
    print("-" * 50)
    print(f"Half-normal mean (black):   {viz_thresholds['random']:.2f}")
    print(f"Liberal threshold (yellow):  {viz_thresholds['liberal']:.2f}")
    print(f"Moderate threshold (orange): {viz_thresholds['moderate']:.2f}")
    print(f"Conservative threshold (red): {viz_thresholds['conservative']:.2f}")

if __name__ == "__main__":
    main()
