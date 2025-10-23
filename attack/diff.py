#!/usr/bin/env python3

"""
Copyright (C) 2024 Hosein Hadipour
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

In case you use this tool please include the above copyright
information (name, contact, license)

Disclaimer: We acknowledge that the Orthros block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of Orthros against differential and differential-linear cryptanalysis.
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import yaml
import time
from gurobipy import *
import math
import os

class Diff:
    """
    This class is used to find differential trail as well as
    computing the differential effect of Orthros block cipher.

    x_roundNumber_bitNumber
    x_roundNumber_0: msb

    [x]-- |S-box| --> [y] --> |Permutation| --> [z] --> |MixColumns| --> [x_next_round]

    """

    diff_count = 0

    def __init__(self, params) -> None:
        Diff.diff_count += 1

        self.nrounds = params["nrounds"]
        self.branch_type = params["branchtype"]
        self.offset = params["offset"]
        self.time_limit = params["timelimit"]
        self.start_weight = params["startweight"]
        self.end_weight = params["endweight"]
        self.fixed_variables = params['fixedVariables']
        self.mode = params['mode']
        self.number_of_trails = params["numberoftrails"]
        self.eps = 1e-3
        self.lp_file_name = f"orthros_nr_{self.nrounds}.lp"
        self.result_file_name = f"result_nr_{self.nrounds}.txt"
        self.milp_variables = []

        # Define some constant lookup tables

        self.bit_permutation_1 = [6,46,62,126,70,52,28,14,36,125,72,83,106,95,4,35,
                                25,41,10,76,87,74,120,42,88,21,11,67,64,38,112,50,
                                85,109,24,65,99,0,49,37,8,66,114,47,127,100,56,40,
                                13,117,78,86,92,58,124,101,55,89,97,9,18,116,59,15,
                                20,45,75,2,77,27,1,60,115,107,26,69,119,3,84,51,
                                123,110,31,82,113,53,81,102,63,118,93,12,30,94,108,32,
                                5,111,29,43,91,19,79,33,73,44,98,48,22,61,68,105,
                                34,71,54,104,17,57,80,103,96,121,23,39,122,90,7,16]
        
        self.bit_permutation_2 = [20,122,74,62,119,35,15,66,9,85,32,117,21,83,127,106,
                                11,98,115,59,71,90,56,26,2,44,103,121,114,107,68,16,
                                84,1,102,33,80,52,76,36,27,94,37,55,82,12,112,64,
                                105,14,91,17,108,124,6,93,29,86,123,79,72,53,19,99,
                                50,18,81,73,67,88,4,61,111,49,24,45,57,78,100,22,
                                110,47,116,54,60,70,97,39,3,41,48,96,23,42,113,87,
                                126,13,31,40,51,25,65,125,8,101,118,28,38,89,5,104,
                                109,120,69,43,7,77,58,34,10,63,30,95,75,46,0,92]
        
        self.nibble_permutation_1 = [10, 27, 5, 1, 30, 23, 16, 13, 21, 31, 6, 14, 0, 25, 11, 18, 15, 28, 19, 24, 7, 8, 22, 3, 4, 29, 9, 2, 26, 20, 12, 17]

        self.nibble_permutation_2 = [26, 13, 7, 11, 29, 0, 17, 21, 23, 5, 18, 25, 12, 10, 28, 2, 14, 19, 24, 22, 1, 8, 4, 31, 15, 6, 27, 9, 16, 30, 20, 3]

        """
        We use S-box Analyzer to encode the DDT of the S-box
        sage: from sboxanalyzer import *
        sage: from sage.crypto.sboxes import SBox
        sage: sa = SboxAnalyzer([0x1,0x0,0x2,0x4,0x3,0x8,0x6,0xd,0x9,0xa,0xb,0xe,0xf,0xc,0x7,0x5])
        sage: cnf, milp = sa.minimized_diff_constraints()
        Simplifying the MILP/SAT constraints ...
        Time used to simplify the constraints: 0.02 seconds
        Number of constraints: 64
        Input:	a0||a1||a2||a3; a0: msb
        Output:	b0||b1||b2||b3; b0: msb
        Weight: 3.0000 p0 + 2.0000 p1
        sage: pretty_print(milp)
        ['- p0 - p1 >= -1',
        '- a0 - a1 + a2 + p0 >= -1',
        'a0 + a3 - b0 + p0 >= 0',
        'a1 - a2 + b1 + p0 >= 0',
        '- a1 + b0 + b1 + p0 >= 0',
        '- a1 - b0 - b2 + p0 >= -2',
        'a1 - a3 + b2 + p0 >= 0',
        '- a0 + b1 - b3 + p0 >= -1',
        'a3 - b2 - b3 + p0 >= -1',
        'a2 - a3 + b3 + p0 >= 0',
        'b0 + b1 + b3 - p1 >= 0',
        'a0 + a2 + a3 + b0 - b3 >= 0',
        'a0 + a1 - b0 - b1 - b3 >= -2',
        'a2 + a3 - b1 - b2 - b3 >= -2',
        '- a0 - a3 + b0 - b2 + p0 >= -2',
        'a1 + b0 - b1 - b2 + p0 >= -1',
        '- a0 + a3 - b1 + b2 + p0 >= -1',
        '- a1 - a2 + b2 - b3 + p0 >= -2',
        '- a0 - a3 - b0 + b3 + p0 >= -2',
        'a2 - b1 + b2 + b3 + p0 >= 0',
        'a0 + a1 + a2 - b0 + p1 >= 0',
        'a1 + a2 + a3 - b2 + p1 >= 0',
        'a1 - b0 - b1 - b2 + p1 >= -2',
        'a0 - a3 + b2 + b3 + p1 >= 0',
        'a0 - a1 - a2 - b0 + b1 - b3 >= -3',
        '- a0 - a1 + a3 - b0 + b1 - b3 >= -3',
        '- a1 - a2 + a3 + b1 - b2 - b3 >= -3',
        '- a0 - a2 - a3 + b1 + b2 - b3 >= -3',
        '- a1 - a2 - a3 + b0 + b1 + b3 >= -2',
        '- a0 + a2 - a3 + b0 + b1 + b3 >= -1',
        '- a0 + a1 + a3 + b0 + b1 + b3 >= 0',
        '- a0 + a1 + a2 - b1 - b2 + b3 >= -2',
        'a0 - a1 - a2 + b1 + b2 + b3 >= -1',
        'a0 - a1 - a2 + b0 + b3 + p0 >= -1',
        '- a0 - a1 - a2 + a3 + b1 + p1 >= -2',
        '- a1 - a3 + b0 - b1 + b2 + p1 >= -2',
        'a0 - a3 - b0 + b1 + b2 + p1 >= -1',
        'a1 + a2 - b0 - b1 - b3 + p1 >= -2',
        'a1 + a3 + b0 - b1 - b3 + p1 >= -1',
        '- a0 - a1 - a2 + a3 + b3 + p1 >= -2',
        'a0 + a1 - a2 - b1 + b3 + p1 >= -1',
        'a1 - a2 - b0 - b1 + b3 + p1 >= -2',
        '- a0 - a1 - a2 + b1 + b3 + p1 >= -2',
        'a3 + b0 - b1 - b2 + b3 + p1 >= -1',
        '- a0 - a1 + b0 + b2 + b3 + p1 >= -1',
        '- a0 + a3 + b1 + b2 + b3 + p1 >= 0',
        'a0 - a1 + a2 - a3 + b0 - b1 - b2 >= -3',
        'a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -3',
        '- a0 - a1 - a2 - a3 + b0 - b2 + b3 >= -4',
        'a0 - a1 - a3 + b0 - b1 - b3 + p1 >= -3',
        '- a2 - a3 - b0 - b1 - b2 - b3 + p1 >= -5',
        'a0 + a3 + b0 - b1 + b2 - b3 + p1 >= -1',
        'a0 - a2 - a3 - b0 - b1 + b3 + p1 >= -3',
        '- a0 + a1 + a2 + b1 + b2 + b3 + p1 >= 0',
        'a0 + a1 + a3 + b0 + b2 - p0 + p1 >= 0',
        'a0 + a1 + a2 - a3 + b0 - b2 - b3 + p1 >= -2',
        '- a1 + a2 - a3 - b0 + b1 - b2 - p0 + p1 >= -4',
        'a0 - a2 + a3 - b0 + b1 - b2 - p0 + p1 >= -3',
        'a0 - a1 + a2 + b0 + b1 + b2 - p0 + p1 >= -1',
        '- a0 + a1 - a2 + b0 - b2 - b3 - p0 + p1 >= -4',
        'a0 - a2 - a3 + b0 + b1 - b2 - b3 - p0 + p1 >= -4',
        '- a0 - a1 + a2 - a3 - b0 - b1 + b3 - p0 + p1 >= -5',
        'a0 - a1 + a2 + a3 - b0 - b2 + b3 - p0 + p1 >= -3',
        '- a0 + a1 - a3 + b0 + b2 >= -1']

        """
        self.sbox_inequalities = ['- p0 - p1 >= -1',
                                '- a0 - a1 + a2 + p0 >= -1',
                                'a0 + a3 - b0 + p0 >= 0',
                                'a1 - a2 + b1 + p0 >= 0',
                                '- a1 + b0 + b1 + p0 >= 0',
                                '- a1 - b0 - b2 + p0 >= -2',
                                'a1 - a3 + b2 + p0 >= 0',
                                '- a0 + b1 - b3 + p0 >= -1',
                                'a3 - b2 - b3 + p0 >= -1',
                                'a2 - a3 + b3 + p0 >= 0',
                                'b0 + b1 + b3 - p1 >= 0',
                                'a0 + a2 + a3 + b0 - b3 >= 0',
                                'a0 + a1 - b0 - b1 - b3 >= -2',
                                'a2 + a3 - b1 - b2 - b3 >= -2',
                                '- a0 - a3 + b0 - b2 + p0 >= -2',
                                'a1 + b0 - b1 - b2 + p0 >= -1',
                                '- a0 + a3 - b1 + b2 + p0 >= -1',
                                '- a1 - a2 + b2 - b3 + p0 >= -2',
                                '- a0 - a3 - b0 + b3 + p0 >= -2',
                                'a2 - b1 + b2 + b3 + p0 >= 0',
                                'a0 + a1 + a2 - b0 + p1 >= 0',
                                'a1 + a2 + a3 - b2 + p1 >= 0',
                                'a1 - b0 - b1 - b2 + p1 >= -2',
                                'a0 - a3 + b2 + b3 + p1 >= 0',
                                'a0 - a1 - a2 - b0 + b1 - b3 >= -3',
                                '- a0 - a1 + a3 - b0 + b1 - b3 >= -3',
                                '- a1 - a2 + a3 + b1 - b2 - b3 >= -3',
                                '- a0 - a2 - a3 + b1 + b2 - b3 >= -3',
                                '- a1 - a2 - a3 + b0 + b1 + b3 >= -2',
                                '- a0 + a2 - a3 + b0 + b1 + b3 >= -1',
                                '- a0 + a1 + a3 + b0 + b1 + b3 >= 0',
                                '- a0 + a1 + a2 - b1 - b2 + b3 >= -2',
                                'a0 - a1 - a2 + b1 + b2 + b3 >= -1',
                                'a0 - a1 - a2 + b0 + b3 + p0 >= -1',
                                '- a0 - a1 - a2 + a3 + b1 + p1 >= -2',
                                '- a1 - a3 + b0 - b1 + b2 + p1 >= -2',
                                'a0 - a3 - b0 + b1 + b2 + p1 >= -1',
                                'a1 + a2 - b0 - b1 - b3 + p1 >= -2',
                                'a1 + a3 + b0 - b1 - b3 + p1 >= -1',
                                '- a0 - a1 - a2 + a3 + b3 + p1 >= -2',
                                'a0 + a1 - a2 - b1 + b3 + p1 >= -1',
                                'a1 - a2 - b0 - b1 + b3 + p1 >= -2',
                                '- a0 - a1 - a2 + b1 + b3 + p1 >= -2',
                                'a3 + b0 - b1 - b2 + b3 + p1 >= -1',
                                '- a0 - a1 + b0 + b2 + b3 + p1 >= -1',
                                '- a0 + a3 + b1 + b2 + b3 + p1 >= 0',
                                'a0 - a1 + a2 - a3 + b0 - b1 - b2 >= -3',
                                'a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -3',
                                '- a0 - a1 - a2 - a3 + b0 - b2 + b3 >= -4',
                                'a0 - a1 - a3 + b0 - b1 - b3 + p1 >= -3',
                                '- a2 - a3 - b0 - b1 - b2 - b3 + p1 >= -5',
                                'a0 + a3 + b0 - b1 + b2 - b3 + p1 >= -1',
                                'a0 - a2 - a3 - b0 - b1 + b3 + p1 >= -3',
                                '- a0 + a1 + a2 + b1 + b2 + b3 + p1 >= 0',
                                'a0 + a1 + a3 + b0 + b2 - p0 + p1 >= 0',
                                'a0 + a1 + a2 - a3 + b0 - b2 - b3 + p1 >= -2',
                                '- a1 + a2 - a3 - b0 + b1 - b2 - p0 + p1 >= -4',
                                'a0 - a2 + a3 - b0 + b1 - b2 - p0 + p1 >= -3',
                                'a0 - a1 + a2 + b0 + b1 + b2 - p0 + p1 >= -1',
                                '- a0 + a1 - a2 + b0 - b2 - b3 - p0 + p1 >= -4',
                                'a0 - a2 - a3 + b0 + b1 - b2 - b3 - p0 + p1 >= -4',
                                '- a0 - a1 + a2 - a3 - b0 - b1 + b3 - p0 + p1 >= -5',
                                'a0 - a1 + a2 + a3 - b0 - b2 + b3 - p0 + p1 >= -3',
                                '- a0 + a1 - a3 + b0 + b2 >= -1']


    @staticmethod
    def ordered_set(seq):
        """
        This method eliminates duplicated elements in a given list,
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    
    @staticmethod
    def flatten_state(state):
        """
        Flatten the state variables
        """

        return [item for sublist in state for item in sublist]

    def generate_round_variables(self, prefix='x', rn=0):
        """
        Generate the input variables of rn'th round
        """

        state_var = [f"{prefix}_{rn}_{bit}" for bit in range(128)]
        self.milp_variables.extend(state_var)
        return state_var

    def generate_round_pr_variables(self, rn):
        """
        Generate the variables encoding the probability of S-boxes
        """

        pr = [[f"pr_{rn}_{nibble}_{bit}" for bit in range(2)] for nibble in range(32)]
        self.milp_variables.extend(self.flatten_state(pr))
        return pr

    def constraints_by_equality(self, a, b):
        """
        Generate constraints for equality
        a = b
        """
        constraint = f"{a} - {b} = 0\n"
        return constraint

    def constraints_by_nibble_equality(self, a, b):
        """
        Generate constraints corresponding
        to equality of two nibbles
        """

        constraints = ""
        for bit in range(4):
            constraints += f"{a[bit]} - {b[bit]} = 0\n"
        return constraints

    def constraints_by_xor(self, a, b, c):
        """
        a + b = c
        model:
        - a - b - c >= -2
          a + b - c >= 0
          a - b + c >= 0
        - a + b + c >= 0
        """

        constraints = f"- {a} - {b} - {c} >= -2\n"
        constraints += f"{a} + {b} - {c} >= 0\n"
        constraints += f"{a} - {b} + {c} >= 0\n"
        constraints += f"- {a} + {b} + {c} >= 0\n"
        return constraints

    def constraints_by_nibble_xor(self, a, b, c):
        """
        Generate constraints for XOR of nibbles
        """

        constraints = ""
        for bit in range(4):
            constraints += self.constraints_by_xor(a[bit], b[bit], c[bit])
        return constraints

    def constraints_by_xor3(self, a0, a1, a2, b):
        """
        Generate the constraints of a three-input XOR  (b = a0 xor a1 xor a2)    
        b - a2 - a1 - a0 >= -2
        - b + a2 - a1 - a0 >= -2
        - b - a2 + a1 - a0 >= -2
        b + a2 + a1 - a0 >= 0
        - b - a2 - a1 + a0 >= -2
        b + a2 - a1 + a0 >= 0
        b - a2 + a1 + a0 >= 0
        - b + a2 + a1 + a0 >= 0
        """

        constraints = f"{b} - {a2} - {a1} - {a0} >= -2\n"
        constraints += f"- {b} + {a2} - {a1} - {a0} >= -2\n"
        constraints += f"- {b} - {a2} + {a1} - {a0} >= -2\n"
        constraints += f"{b} + {a2} + {a1} - {a0} >= 0\n"
        constraints += f"- {b} - {a2} - {a1} + {a0} >= -2\n"
        constraints += f"{b} + {a2} - {a1} + {a0} >= 0\n"
        constraints += f"{b} - {a2} + {a1} + {a0} >= 0\n"
        constraints += f"- {b} + {a2} + {a1} + {a0} >= 0\n"
        return constraints
    
    def constraints_by_nibble_xor3(self, a0, a1, a2, b):
        """
        Generate the constraints of a three-input XOR  (b = a0 xor a1 xor a2), 
        where a0, a1, a2, and b are nibbles  
        """

        constraints = ""
        for bit in range(4):
            constraints += self.constraints_by_xor3(a0[bit], a1[bit], a2[bit], b[bit])
        return constraints

    def constraints_by_sbox(self, di, do, pr):
        """
        Generate constraints modeling the DDT of S-box

        :param str[4] di: input difference
        :param str[4] do: output difference
        :param str[3] pr: probability of (di --> do) such that
                          hamming_weight(pr) = -log2(pr(di --> do))
        :return constraints encoding the DDT of S-box:
        :rtype str:
        """

        constraints = ""
        for ineq in self.sbox_inequalities:
            temp = ineq
            for i in range(4):
                temp = temp.replace(f"a{i}", di[i])
            for i in range(4):
                temp = temp.replace(f"b{i}", do[i])
            for i in range(2):
                temp = temp.replace(f"p{i}", pr[i])
            constraints += temp + "\n"
        return constraints

    def generate_objective_function(self):
        """
        Generate the objective function of MILP model
        The objective is minimizing the summation of variables
        which encode the weight (or probability exponent) the
        differential trail
        """

        objective_function = "minimize\n"
        weight = []
        for r in range(self.nrounds):
            pr_vars = self.generate_round_pr_variables(rn=r)
            weight += [f"3 {pr_vars[nibble][0]} + 2 {pr_vars[nibble][1]}" for nibble in range(32)]
        weight = " + ".join(weight)
        objective_function += weight + "\n"
        return objective_function

    def generate_constraints(self):
        """
        Generate the constraints describing the propagation
        of differential trails through a reduced-round Orthros
        """

        constraints = "subject to\n"
        for rn in range(self.nrounds):
            x_in = self.generate_round_variables(rn=rn, prefix='x')
            pr = self.generate_round_pr_variables(rn=rn)
            y = self.generate_round_variables(rn=rn, prefix='y')
            z = self.generate_round_variables(rn=rn, prefix='z')
            x_out = self.generate_round_variables(rn=rn + 1, prefix='x')
            # model the S-box layer
            for nibble in range(32):
                constraints += self.constraints_by_sbox(di=x_in[4*nibble:4*nibble + 4], do=y[4*nibble:4*nibble + 4], pr=pr[nibble])
            # model the permuation layer
            if (self.offset + rn) < 4:
                for bit in range(128):
                    if self.branch_type == 0:
                        constraints += self.constraints_by_equality(z[self.bit_permutation_1[bit]], y[bit])
                    else:
                        constraints += self.constraints_by_equality(z[self.bit_permutation_2[bit]], y[bit])
            else:
                for nibble in range(32):
                    for bit in range(4):
                        if self.branch_type == 0:
                            constraints += self.constraints_by_equality(z[4*self.nibble_permutation_1[nibble] + bit], y[4*nibble + bit])
                        else:
                            constraints += self.constraints_by_equality(z[4*self.nibble_permutation_2[nibble] + bit], y[4*nibble + bit])
            # model the MixColumns layer
            for i in range(8):
                constraints += self.constraints_by_nibble_xor3(z[16*i + 4*1: 16*i + 4*1 + 4], z[16*i + 4*2: 16*i + 4*2 + 4], z[16*i + 4*3: 16*i + 4*3 + 4], x_out[16*i + 4*0: 16*i + 4*0 + 4])
                constraints += self.constraints_by_nibble_xor3(z[16*i + 4*0: 16*i + 4*0 + 4], z[16*i + 4*2: 16*i + 4*2 + 4], z[16*i + 4*3: 16*i + 4*3 + 4], x_out[16*i + 4*1: 16*i + 4*1 + 4])
                constraints += self.constraints_by_nibble_xor3(z[16*i + 4*0: 16*i + 4*0 + 4], z[16*i + 4*1: 16*i + 4*1 + 4], z[16*i + 4*3: 16*i + 4*3 + 4], x_out[16*i + 4*2: 16*i + 4*2 + 4])
                constraints += self.constraints_by_nibble_xor3(z[16*i + 4*0: 16*i + 4*0 + 4], z[16*i + 4*1: 16*i + 4*1 + 4], z[16*i + 4*2: 16*i + 4*2 + 4], x_out[16*i + 4*3: 16*i + 4*3 + 4])
        return constraints

    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """

        self.milp_variables = self.ordered_set(self.milp_variables)
        constraints = "Binary\n"
        constraints += "\n".join(self.milp_variables) + "\n"
        return constraints

    def exclude_trivial_trail(self):
        """
        Exclude all-zero solution from the solution space
        """

        input_diff = self.generate_round_variables(rn=0, prefix='x')
        input_diff = " + ".join(input_diff)
        constraint = f"{input_diff} >= 1\n"
        return constraint

    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            if len(var) == 2:
                state_vars = self.generate_round_variables(rn=var[1], prefix=var[0])
                state_values = list(bin(int(val, 16))[2:].zfill(128))
                for i in range(128):
                    lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
            elif len(var) == 3:
                state_vars = [f"{var[0]}_{var[1]}_{var[2]}"]
                lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
            else:
                pass
        return lp_contents

    def make_model(self):
        """
        Build the MILP model to find the best differential trail
        """

        lp_contents = "\\ Differential attack on {} rounds of Orthros\n".format(self.nrounds)
        lp_contents += self.generate_objective_function()
        lp_contents += self.generate_constraints()
        lp_contents += self.exclude_trivial_trail()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_binary_vars()
        lp_contents += "end"
        with open(self.lp_file_name, "w") as lp_file:
            lp_file.write(lp_contents)

    def exclude_the_previous_sol(self):
        '''
        Let x{S} be the binary variables. Suppose you have a binary
        solution x* in available from the most recent optimization.
        Let N be the subset of S such that x*[n] = 1 for all n in N
        Then, add the following constraint:
        sum{n in N} x[n] - sum{s in S-N} x[s] <= |N|-1
        '''

        all_vars = self.milp_model.getVars()
        nonzero_vars = [v for v in all_vars if v.x == 1]
        zero_vars = [v for v in all_vars if v.x == 0]
        support = len(nonzero_vars)
        first_term = sum(nonzero_vars)
        second_term = sum(zero_vars)
        lhs = first_term - second_term
        self.milp_model.addConstr(lhs <= support - 1)

    def solve(self):
        output = None
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        if self.mode == 0:
            output = self.find_characteristic()
        elif self.mode == 1:
            self.find_multiple_characteristics(self.number_of_trails)
        elif self.mode == 2:
            output = self.compute_differential_effect()
            # self.compute_differential_effect_classic_method()
        else:
            print('Enter a number in [0, 1, 2], for the mode parameter please!')
        return output

    def parse_solver_output(self):
        """
        Extract the differential characteristic from the solver output
        """

        characteristic = dict()
        for r in range(self.nrounds + 1):
            x = self.generate_round_variables(rn=r, prefix='x')
            x_value = hex(int("0b" + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x))), 2))[2:].zfill(32)
            characteristic[f"x_{r}"] = x_value
        for r in range(self.nrounds):
            round_probability = 0
            for nibble in range(32):
                round_probability += sum([3*int(self.milp_model.getVarByName(f"pr_{r}_{nibble}_0").Xn) +\
                                          2*int(self.milp_model.getVarByName(f"pr_{r}_{nibble}_1").Xn)])
            characteristic[f"pr_{r}"] = f"-{round_probability}"
        characteristic["total_weight"] = "%0.02f" % self.total_weight
        characteristic["nrounds"] = self.nrounds
        return characteristic

    @staticmethod
    def print_trail(trail):
        """
        Print out the discovered differential characteristic
        """

        header = ['x', 'pr']
        # Print everthing
        trail_values = map(str, trail.values())
        col_width = max(len(s) for s in trail_values) + 2
        header_str = "Rounds\t"
        data_str = ""
        current_row = 0
        for entry in header[0:-2]:
            header_str += entry.ljust(col_width)
        header_str += header[-2].ljust(col_width)
        header_str += header[-1].ljust(7)
        for r in range(trail["nrounds"] + 1):
            data_str += str(current_row) + '\t'
            data_str += trail.get(f"x_{r}", 'none').ljust(col_width)
            data_str += trail.get(f"pr_{r}", 'none').ljust(col_width)
            data_str += '\n'
            current_row += 1
        stroutput = header_str
        stroutput += "\n" + "-"*len(header_str) + "\n"
        stroutput += data_str
        total_weight = trail["total_weight"]
        stroutput += f"Weight: -{total_weight}" + "\n"
        return stroutput

    def find_characteristic(self):
        """
        Find the best differential trail for reduced-round Orthros
        """
        diff_trail = None
        self.milp_model.Params.OutputFlag = True
        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        #m.setParam(GRB.Param.Threads, 16)
        self.milp_model.optimize()
        # Gurobi syntax: m.Status == 2 represents the model is feasible.
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            self.total_weight = self.milp_model.objVal
            print(f"\nThe probability of the best differential characteristic: 2^-({self.total_weight})")
            print("\nDifferential trail:\n")
            diff_trail = self.parse_solver_output()
            stroutput = self.print_trail(trail=diff_trail)
            print(stroutput)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.Status.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Time used: %0.02f" % elapsed_time)
        return diff_trail

    def find_multiple_characteristics(self, number_of_trails=2):
        """
        Find multiple differential trails for reduced-round of Orthros
        """
        self.milp_model.Params.PreSolve = 1
        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        self.milp_model.Params.OutputFlag = False
        self.milp_model.Params.PoolSearchMode = 2
        # Limit number of solutions
        self.milp_model.Params.PoolSolutions = number_of_trails
        time_start = time.time()
        self.milp_model.optimize()
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            # First Method:
            for sol_number in range(number_of_trails):
                if (self.milp_model.Status == GRB.OPTIMAL):
                    self.total_weight = self.milp_model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    output = self.print_trail(trail=diff_trail)
                    print(output)
                elif (self.milp_model.Status == GRB.TIME_LIMIT or self.milp_model.Status == GRB.INTERRUPTED):
                    self.total_weight = self.milp_model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    output = self.print_trail(trail=diff_trail)
                    print(output)
                    break
                else:
                    break
                self.exclude_the_previous_sol()
                print("#"*50)
                self.milp_model.optimize()
            # Second Method:
            # number_of_trails = self.milp_model.SolCount
            # for sol_number in range(number_of_trails):
            #     self.milp_model.Params.SolutionNumber = sol_number
            #     # PoolObjVal : This attribute is used to query the objective value of the <span>$</span>k<span>$</span>-th solution stored in the pool of feasible solutions found so far for the problem
            #     self.total_weight = self.milp_model.PoolObjVal
            #     diff_trail = self.parse_solver_output()
            #     self.print_trail(trail=diff_trail)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Total time to find %s differential trails: %0.02f" % (number_of_trails, elapsed_time))

    def compute_differential_effect(self):
        """
        Compute the differential effect for a given input/output differences
        Some general information about Gurobi:
        PoolSolutions: It controls the size of the solution pool.
        Changing this parameter won't affect the number of solutions that are found -
        it simply determines how many of those are retained
        You can use the PoolSearchMode parameter to control the approach used to find solutions.
        In its default setting (0), the MIP search simply aims to find one optimal solution.
        Setting the parameter to 2 causes the MIP to do a systematic search for the n best solutions.
        With a setting of 2, it will find the n best solutions,
        where n is determined by the value of the PoolSolutions parameter
        SolCount: Number of solutions found during the most recent optimization.

        Model status:
        LOADED	1	Model is loaded, but no solution information is available.
        OPTIMAL	2	Model was solved to optimality (subject to tolerances), and an optimal solution is available.
        INFEASIBLE	3	Model was proven to be infeasible.
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        #self.milp_model.Params.PreSolve = 0 # Activating this flag causes the performance to be decreased, but the accuracy will be increased
        self.milp_model.Params.PoolSearchMode = 2
        self.milp_model.Params.PoolSolutions = 1
        self.milp_model.Params.OutputFlag = False

        self.milp_model.printStats()

        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        current_probability = 0
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.PoolObjVal
                self.milp_model.Params.PoolSolutions = 2000000000 #GRB.MAXINT
                temp_constraint = self.milp_model.addConstr(obj == self.total_weight, name='temp_constraint')
                # self.milp_model.Params.PoolGap = 0
                # self.milp_model.Params.PreSolve = 0
                # self.milp_model.printStats()
                self.milp_model.update()
                self.milp_model.optimize()
                diff_prob += math.pow(2, -self.total_weight) * self.milp_model.SolCount
                print(f"Current weight: {self.total_weight}")
                print(f"Number of trails: {self.milp_model.SolCount}")
                current_probability = math.log(diff_prob, 2)
                print(f"\tCurrent Probability: 2^({current_probability})")
                elapsed_time = time.time() - time_start
                print("Time used = %0.04f seconds\n" % elapsed_time)
                self.milp_model.remove(temp_constraint)
                self.milp_model.Params.PoolSolutions = 1
                self.milp_model.addConstr(obj >= (self.total_weight + self.eps), name='temp_cond')
                #self.milp_model.Params.PreSolve = 0
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print("The model is infeasible!")
        else:
            print("Unknown Error!")
        return current_probability

    def compute_differential_effect_classic_method(self):
        """
        Compute differential effect by enumerating all possible differential trails
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        self.milp_model.Params.OutputFlag = False
        # self.milp_model.printStats()
        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        # self.milp_model.Params.Quad = 1
        sol_dict = dict()
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.objVal
                diff_prob += math.pow(2, -self.total_weight)
                total_weight_st = 'ntrails_%0.2f' % self.total_weight
                sol_dict[total_weight_st] = sol_dict.get(total_weight_st, 0) + 1
                print('Current weight: %s' % str(self.total_weight))
                print('Number of trails: %d' % sol_dict[total_weight_st])
                print('\tCurrent Probability: 2^(' + str(math.log(diff_prob, 2)) + ')')
                time_end = time.time()
                print('Time used = %0.4f seconds\n' % (time_end - time_start))
                self.exclude_the_previous_sol()
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print('The model is infeasible!')
        else:
            print('Unknown Error!')

def loadparameters(args):
        """
        Get parameters from the argument list and inputfile.
        """
        # Load default values
        params = {"nrounds" : 3,
                  "branchtype" : 0,
                  "offset" : 4,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 2,
                  "fixedVariables" : {}}

        # Check if there is an input file specified
        if args.inputfile:
            with open(args.inputfile[0], 'r') as input_file:
                doc = yaml.load(input_file, Loader=yaml.FullLoader)
                params.update(doc)
                if "fixedVariables" in doc:
                    fixed_vars = {}
                    for variable in doc["fixedVariables"]:
                        fixed_vars = dict(list(fixed_vars.items()) +
                                        list(variable.items()))
                    params["fixedVariables"] = fixed_vars

        # Override parameters if they are set on commandline
        if args.nrounds:
            params["nrounds"] = args.nrounds[0]

        if args.branchtype:
            params["branchtype"] = args.branchtype[0]

        if args.offset:
            params["offset"] = args.offset[0]

        if args.startweight:
            params["startweight"] = args.startweight[0]

        if args.endweight:
            params["endweight"] = args.endweight[0]

        if args.mode:
            params["mode"] = args.mode[0]

        if args.timelimit:
            params["timelimit"] = args.timelimit[0]

        if args.numberoftrails:
            params["numberoftrails"] = args.numberoftrails[0]

        return params

def main():
    """
    Parse the arguments and start the request functionality with the provided
    parameters.
    """

    parser = ArgumentParser(description="This tool finds the best differential"
                                        "trail in a cryptographic primitive"
                                        "using Gurobi",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('--nrounds', nargs=1, type=int,
                        help="The number of rounds for the cipher")
    parser.add_argument('--branchtype', nargs=1, type=int,
                        choices=[0, 1], help=
                        "0 = branch type 1\n"
                        "1 = branch type 2")
    parser.add_argument('--offset', nargs=1, type=int,
                        help="The starting round for distinguisher.")
    parser.add_argument('--startweight', nargs=1, type=int,
                        help="Starting weight for the trail search.")
    parser.add_argument('--endweight', nargs=1, type=int,
                        help="Stop search after reaching endweight.")
    parser.add_argument('--mode', nargs=1, type=int,
                        choices=[0, 1], help=
                        "0 = search characteristic for fixed round\n"
                        "1 = determine the probability of the differential\n")
    parser.add_argument('--timelimit', nargs=1, type=int,
                        help="Set a timelimit for the search in seconds.")
    parser.add_argument('--inputfile', nargs=1, help="Use an yaml input file to"
                                                     "read the parameters.")
    parser.add_argument('--numberoftrails', nargs=1, type=int,
                        help="Number of trails.")

    # Parse command line arguments and construct parameter list.
    args = parser.parse_args()
    params = loadparameters(args)
    orthros = Diff(params)
    orthros.make_model()
    orthros.solve()

if __name__ == "__main__":
    main()