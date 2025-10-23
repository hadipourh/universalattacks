#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

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
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import yaml
import time
from gurobipy import read
from gurobipy import GRB
import os
from drawprfmilp import DrawDL
import subprocess
try:
    subprocess.run(['gurobi_cl', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    gurobi_available = True
except FileNotFoundError:
    gurobi_available = False

if gurobi_available:
    print("Gurobi is available")
    from diff import Diff

class DL:
    '''
    Convert the problem of finding differential-linear distinguisher to a MILP problem
    '''

    count = 0
    def __init__(self, param):
        self.RB = param['RB']
        self.RU = param['RU']
        self.RM = param['RM']
        self.RL = param['RL']
        self.WU = param['WU']
        self.WM = param['WM']
        self.WL = param['WL']
        self.offset = param['offset']
        self.compute_differential_effect_using_gurobi = param["differential_effect"]        
        self.number_of_parallel_threads = param['np']
        self.time_limit = param['timelimit']
        self.fixed_variables = param['fixedVariables']
        self.output_file_name = [param["output1"], param["output2"]]

        self.RT = self.RB + self.RU + self.RM + self.RL
        self.block_size = 128
        self.sbox_size = 4
        self.bit_permutation = [[6,46,62,126,70,52,28,14,36,125,72,83,106,95,4,35, 25,41,10,76,87,74,120,42,88,21,11,67,64,38,112,50, 85,109,24,65,99,0,49,37,8,66,114,47,127,100,56,40, 13,117,78,86,92,58,124,101,55,89,97,9,18,116,59,15, 20,45,75,2,77,27,1,60,115,107,26,69,119,3,84,51, 123,110,31,82,113,53,81,102,63,118,93,12,30,94,108,32, 5,111,29,43,91,19,79,33,73,44,98,48,22,61,68,105, 34,71,54,104,17,57,80,103,96,121,23,39,122,90,7,16],
                                [20,122,74,62,119,35,15,66,9,85,32,117,21,83,127,106, 11,98,115,59,71,90,56,26,2,44,103,121,114,107,68,16, 84,1,102,33,80,52,76,36,27,94,37,55,82,12,112,64, 105,14,91,17,108,124,6,93,29,86,123,79,72,53,19,99, 50,18,81,73,67,88,4,61,111,49,24,45,57,78,100,22, 110,47,116,54,60,70,97,39,3,41,48,96,23,42,113,87, 126,13,31,40,51,25,65,125,8,101,118,28,38,89,5,104, 109,120,69,43,7,77,58,34,10,63,30,95,75,46,0,92]]
    
        self.nibble_permutation = [[10, 27, 5, 1, 30, 23, 16, 13, 21, 31, 6, 14, 0, 25, 11, 18, 15, 28, 19, 24, 7, 8, 22, 3, 4, 29, 9, 2, 26, 20, 12, 17],
                                   [26, 13, 7, 11, 29, 0, 17, 21, 23, 5, 18, 25, 12, 10, 28, 2, 14, 19, 24, 22, 1, 8, 4, 31, 15, 6, 27, 9, 16, 30, 20, 3]]
        
        self.ksch_bit_permutation = [[0,53,87,73,22,95,99,48,61,36,108,1,24,67,119,93,54,103,69,112,16,111,94,122,31,66,33,83,47,3,65,62,123,9,101,19,5,58,89,37,38,51,28,106,82,76,121,4,70,7,42,92,104,80,45,75,114,17,2,97,46,107,63,18,109,15,127,43,13,59,29,125,77,11,50,30,12,90,118,64,20,35,57,10,124,56,68,91,116,21,84,98,52,81,126,34,105,27,120,74,6,85,40,72,113,41,23,49,79,55,102,8,117,39,88,26,25,110,14,32,115,100,86,71,78,44,96,60],
                                     [76,30,53,35,31,46,2,79,11,125,110,87,39,91,14,101,97,118,36,48,29,80,57,115,49,18,74,85,61,82,105,126,70,12,47,111,51,17,66,1,60,96,116,71,81,114,104,15,42,124,100,4,113,44,75,89,23,0,84,107,32,26,88,8,69,121,38,94,37,86,54,21,62,123,41,10,16,95,117,65,45,50,72,20,109,58,7,67,108,28,3,55,92,103,24,5,77,9,27,102,122,6,106,22,99,34,90,56,43,83,120,64,78,59,119,93,40,98,52,68,112,33,63,25,19,73,127,13]]
        
        self.inv_ksch_bit_permutation = [[0, 11, 58, 29, 47, 36, 100, 49, 111, 33, 83, 73, 76, 68, 118, 65, 20, 57, 63, 35, 80, 89, 4, 106, 12, 116, 115, 97, 42, 70, 75, 24, 119, 26, 95, 81, 9, 39, 40, 113, 102, 105, 50, 67, 125, 54, 60, 28, 7, 107, 74, 41, 92, 1, 16, 109, 85, 82, 37, 69, 127, 8, 31, 62, 79, 30, 25, 13, 86, 18, 48, 123, 103, 3, 99, 55, 45, 72, 124, 108, 53, 93, 44, 27, 90, 101, 122, 2, 114, 38, 77, 87, 51, 15, 22, 5, 126, 59, 91, 6, 121, 34, 110, 17, 52, 96, 43, 61, 10, 64, 117, 21, 19, 104, 56, 120, 88, 112, 78, 14, 98, 46, 23, 32, 84, 71, 94, 66],
                                         [57, 39, 6, 90, 51, 95, 101, 86, 63, 97, 75, 8, 33, 127, 14, 47, 76, 37, 25, 124, 83, 71, 103, 56, 94, 123, 61, 98, 89, 20, 1, 4, 60, 121, 105, 3, 18, 68, 66, 12, 116, 74, 48, 108, 53, 80, 5, 34, 19, 24, 81, 36, 118, 2, 70, 91, 107, 22, 85, 113, 40, 28, 72, 122, 111, 79, 38, 87, 119, 64, 32, 43, 82, 125, 26, 54, 0, 96, 112, 7, 21, 44, 29, 109, 58, 27, 69, 11, 62, 55, 106, 13, 92, 115, 67, 77, 41, 16, 117, 104, 50, 15, 99, 93, 46, 30, 102, 59, 88, 84, 10, 35, 120, 52, 45, 23, 42, 78, 17, 114, 110, 65, 100, 73, 49, 9, 31, 126]]

        self.inv_kperm_at_round = [[[0 for _ in range(self.block_size)] for _ in range(2)] for _ in range(self.offset + 1)]
        for branch in range(2):
            for bit in range(self.block_size):
                self.inv_kperm_at_round[0][branch][bit] = self.inv_ksch_bit_permutation[branch][bit]
        for r in range(1, self.offset + 1):
            for branch in range(2):
                for bit in range(self.block_size):
                    self.inv_kperm_at_round[r][branch][bit] = self.inv_ksch_bit_permutation[branch][self.inv_kperm_at_round[r - 1][branch][bit]]

        self.translate_joint_vars_to_single_vars = {(0, 0): 0, (0, 1): 1, (1, 0): -1}
        self.used_variables = []
        self.lp_constraints_for_joint_variables = []
        self.model_filename = f"dl-{self.RT}r.lp"
    
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 64
        # Input:	a0||a1||a2||a3; a0: msb
        # Output:	b0||b1||b2||b3; b0: msb
        # Weight: 3.0000 p0 + 2.0000 p1
        self.ddt_constraints = ['- p0 - p1 >= -1',
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
        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 54
        # Input:	a0||a1||a2||a3; a0: msb
        # Output:	b0||b1||b2||b3; b0: msb
        self.star_ddt_constraints = ['a0 + a1 + a3 - b0 - b1 >= -1',
                                    '- a0 + a1 + a3 + b0 - b2 >= -1',
                                    'a1 + a2 + a3 + b0 - b2 >= 0',
                                    'a0 + a1 + a2 - b0 + b2 >= 0',
                                    'a0 + a1 - a3 - b0 + b2 >= -1',
                                    'a0 + a2 + a3 + b0 - b3 >= 0',
                                    'a0 + a1 - b0 - b1 - b3 >= -2',
                                    'a1 + a2 + a3 - b2 - b3 >= -1',
                                    'a0 + a1 + a2 - b0 + b3 >= 0',
                                    'a0 - a3 + b0 + b2 + b3 >= 0',
                                    '- a0 + b0 + b1 + b2 + b3 >= 0',
                                    '- a0 - a1 - a2 + a3 + b0 + b1 >= -2',
                                    '- a1 + a2 - a3 - b0 + b1 - b2 >= -3',
                                    'a0 - a2 + a3 - b0 + b1 - b2 >= -2',
                                    'a0 - a1 - a2 - b0 + b1 - b3 >= -3',
                                    '- a0 - a1 + a3 - b0 + b1 - b3 >= -3',
                                    '- a0 + a1 - a2 + b0 - b2 - b3 >= -3',
                                    '- a1 + a2 + a3 - b1 - b2 - b3 >= -3',
                                    'a0 + a2 + b0 - b1 - b2 - b3 >= -2',
                                    'a1 + a2 - b0 - b1 + b2 - b3 >= -2',
                                    '- a0 + a1 + b0 - b1 + b2 - b3 >= -2',
                                    '- a0 - a2 - a3 + b1 + b2 - b3 >= -3',
                                    '- a0 - a1 - a2 - a3 + b1 + b3 >= -3',
                                    'a0 + a2 - a3 - b0 + b1 + b3 >= -1',
                                    '- a1 - a2 - a3 + b0 + b1 + b3 >= -2',
                                    '- a0 + a2 - a3 + b0 + b1 + b3 >= -1',
                                    'a0 + a2 + a3 - b0 - b2 + b3 >= -1',
                                    '- a0 + a1 + a2 - b1 - b2 + b3 >= -2',
                                    'a0 + a1 + a2 - b1 + b2 + b3 >= 0',
                                    '- a1 + a2 - a3 - b1 + b2 + b3 >= -2',
                                    'a0 - a1 - a2 + b1 + b2 + b3 >= -1',
                                    'a0 + a1 - a2 + a3 + b0 + b1 + b2 >= 0',
                                    'a0 - a2 - a3 + b0 + b1 - b2 - b3 >= -3',
                                    'a0 - a1 + a3 + b0 + b1 - b2 - b3 >= -2',
                                    'a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -3',
                                    'a0 - a1 - a2 + b0 - b1 + b2 - b3 >= -3',
                                    '- a0 - a1 - a3 + b0 - b1 + b2 - b3 >= -4',
                                    '- a0 + a1 - a2 - a3 - b0 - b1 + b3 >= -4',
                                    'a0 - a1 + a2 - a3 + b0 - b1 + b3 >= -2',
                                    '- a0 - a1 + a2 + a3 + b0 - b1 + b3 >= -2',
                                    '- a0 - a1 - a2 + a3 - b0 - b2 + b3 >= -4',
                                    '- a0 + a2 - a3 - b0 - b1 - b2 + b3 >= -4',
                                    'a0 + a1 - a2 + b0 - b1 - b2 + b3 >= -2',
                                    '- a0 + a1 - a2 + a3 - b0 + b2 + b3 >= -2',
                                    '- a0 - a1 - a2 + a3 - b1 + b2 + b3 >= -3',
                                    '- a0 - a1 + a2 + a3 + b1 + b2 + b3 >= -1',
                                    '- a0 - a1 - a2 - a3 + b0 - b2 + b3 >= -4',
                                    '- a1 - a2 - a3 - b0 - b1 - b2 - b3 >= -6',
                                    'a0 - a1 - a2 - a3 - b0 - b1 - b2 >= -5',
                                    'a0 - a2 + a3 + b0 - b1 - b2 + b3 >= -2',
                                    'a1 + a2 - a3 + b1 + b2 + b3 >= 0',
                                    'a0 - a1 + a2 + b0 + b1 + b2 >= 0',
                                    'a1 + a3 - b1 - b2 - b3 >= -2',
                                    '- a0 + a1 - a3 + b0 + b2 >= -1']
        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 71
        # Input:	a0||a1||a2||a3; a0: msb
        # Output:	b0||b1||b2||b3; b0: msb
        # Weight: 4.0000 p0 + 2.0000 p1
        self.square_lat_constraints = ['- p0 - p1 >= -1',
                                    'a0 - b0 + b2 + p0 >= 0',
                                    'a0 + a1 + a2 - b0 - b1 >= -1',
                                    '- a0 - a1 - a2 + b0 + b1 >= -2',
                                    '- a0 - a1 - a2 + b0 + p0 >= -2',
                                    'a0 + a2 - b1 - b2 + p0 >= -1',
                                    'a1 + a2 - b1 - b2 + p0 >= -1',
                                    '- a0 + a3 + b1 - b2 + p0 >= -1',
                                    'a1 - b0 + b1 - b2 + p0 >= -1',
                                    '- a0 + a1 + b0 + b2 + p0 >= 0',
                                    '- a1 - a2 + b1 + b2 + p0 >= -1',
                                    '- a1 + b0 + b1 + b2 + p0 >= 0',
                                    '- a1 + a3 - b0 - b3 + p0 >= -2',
                                    'a2 + a3 + b1 - b3 + p0 >= 0',
                                    'a3 + b0 + b2 - b3 + p0 >= 0',
                                    'a2 - a3 - b0 + b3 + p0 >= -1',
                                    '- a0 + a2 - b1 + b3 + p0 >= -1',
                                    'a0 - a2 + b1 + b3 + p0 >= 0',
                                    '- a0 + a2 - b2 + b3 + p0 >= -1',
                                    'a3 + b0 - b2 + b3 + p0 >= 0',
                                    '- a3 + b0 + b2 + b3 + p0 >= 0',
                                    '- a0 - a1 - a2 - a3 + p1 >= -3',
                                    'a0 + a1 - a3 + b0 - b1 + p0 >= -1',
                                    'a1 - a2 - a3 - b0 + b2 + p0 >= -2',
                                    'a2 - a3 - b0 - b1 + b2 + p0 >= -2',
                                    'a0 + a3 + b0 - b1 - b3 + p0 >= -1',
                                    '- a0 - a1 + a2 + b1 - b3 + p0 >= -2',
                                    'a1 - a2 - a3 - b2 - b3 + p0 >= -3',
                                    'a0 + a1 + b0 - b2 - b3 + p0 >= -1',
                                    '- a1 - a3 + b0 - b2 - b3 + p0 >= -3',
                                    'a1 - a2 + a3 - b0 + b3 + p0 >= -1',
                                    'a1 + a2 + b0 + b2 + b3 - p1 >= 0',
                                    '- a0 - a1 - a2 - b0 - b1 + p1 >= -4',
                                    '- a0 - a1 - a3 - b0 - b2 + p1 >= -4',
                                    'a0 + a1 + a3 - b0 - b2 + p1 >= -1',
                                    '- a0 - a1 - a3 + b0 + b2 + p1 >= -2',
                                    'a0 - a1 + a2 + a3 + b3 + p1 >= 0',
                                    '- a0 - a2 - b0 - b1 + b3 + p1 >= -3',
                                    '- a0 - a3 - b0 - b2 + b3 + p1 >= -3',
                                    'a0 + a3 - b0 - b2 + b3 + p1 >= -1',
                                    '- a2 - a3 - b1 - b2 + b3 + p1 >= -3',
                                    'a0 + a1 + a2 + a3 - p0 + p1 >= 0',
                                    'a0 + a1 - a2 + a3 + b0 + b1 + b2 >= 0',
                                    'a0 + a1 - a2 - a3 - b1 - b2 - b3 >= -4',
                                    'a0 - a1 - a2 - a3 + b1 - b2 - b3 >= -4',
                                    'a0 - a1 - a2 + a3 + b0 + b2 + b3 >= -1',
                                    '- a1 - a2 - a3 - b0 - b1 - b2 + p0 >= -5',
                                    '- a0 + a1 - a2 - b0 - b1 - b3 + p0 >= -4',
                                    'a0 + a1 + a3 + b0 - b1 + b2 + p1 >= 0',
                                    '- a0 - a1 + a2 + a3 + b1 + b2 + p1 >= -1',
                                    'a0 + a1 + a2 + b0 + b1 - b3 + p1 >= 0',
                                    '- a1 - a2 - a3 - b1 + b2 - b3 + p1 >= -4',
                                    'a1 + a2 + a3 - b1 + b2 - b3 + p1 >= -1',
                                    'a0 + a2 - a3 + b0 + b1 + b3 + p1 >= 0',
                                    '- a0 + a2 + a3 + b1 + b2 + b3 + p1 >= 0',
                                    '- a0 + a1 - a2 + a3 + b0 - b1 - b3 + p1 >= -3',
                                    'a0 - a1 + a2 - a3 - b0 + b1 - b3 + p1 >= -3',
                                    '- a0 + a1 - a2 + a3 - b0 + b1 - b3 + p1 >= -3',
                                    '- a0 + a1 + a2 - a3 + b0 - b2 - b3 + p1 >= -3',
                                    'a0 - a1 - a2 + a3 + b0 - b2 - b3 + p1 >= -3',
                                    '- a0 - a1 + a2 + a3 - b1 - b2 - b3 + p1 >= -4',
                                    '- a0 + a1 + a2 - a3 - b0 + b2 - b3 + p1 >= -3',
                                    'a0 + a1 - a2 - a3 + b1 + b2 - b3 + p1 >= -2',
                                    '- a0 + a1 + a2 + a3 + b1 - b2 - b3 - p0 + p1 >= -3',
                                    'a0 - a1 + a2 - a3 + b0 - b1 - b3 >= -3',
                                    'a0 - a1 - a2 + a3 - b0 + b2 - b3 >= -3',
                                    '- a0 - a2 + b0 + b1 + b3 >= -1',
                                    '- a0 - a3 + b0 + b2 + b3 >= -1',
                                    '- a2 - a3 + b1 + b2 + b3 >= -1',
                                    'a2 + a3 - b1 - b2 + b3 >= -1',
                                    'a0 + a2 - b0 - b1 + b3 >= -1']
        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 44
        # Input:	a0||a1||a2||a3; a0: msb
        # Output:	b0||b1||b2||b3; b0: msb
        self.star_lat_constraints = ['a0 + a1 + a2 + a3 - b0 >= 0',
                                    'a0 + a1 + a2 + a3 - b3 >= 0',
                                    'a0 + a1 - a2 + a3 + b0 + b1 + b2 >= 0',
                                    '- a0 - a1 - a2 + a3 - b0 - b1 - b3 >= -5',
                                    '- a0 - a1 - a2 - b0 - b1 - b2 - b3 >= -6',
                                    'a0 + a1 - a3 - b0 - b1 - b2 - b3 >= -4',
                                    '- a0 - a1 + a3 - b0 - b1 - b2 - b3 >= -5',
                                    '- a0 + a1 + a2 - a3 - b0 - b2 + b3 >= -3',
                                    'a1 - a2 + a3 - b0 - b1 - b2 + b3 >= -3',
                                    '- a0 + a1 - a3 - b0 + b1 - b2 + b3 >= -3',
                                    'a0 + a2 + a3 + b0 + b1 - b2 + b3 >= 0',
                                    '- a0 + a1 - a2 - b0 - b1 + b2 + b3 >= -3',
                                    'a0 - a1 + a2 + a3 + b1 + b2 + b3 >= 0',
                                    '- a0 + a1 + a2 - a3 + b0 - b1 - b2 - b3 >= -4',
                                    'a0 - a1 + a2 - a3 - b0 + b1 + b2 - b3 >= -3',
                                    'a0 - a1 - a2 + a3 + b0 - b1 - b2 - b3 >= -4',
                                    '- a0 - a1 + a2 + a3 + b1 + b2 - b3 >= -2',
                                    'a1 + a2 - a3 - b0 - b1 + b2 - b3 >= -3',
                                    'a0 - a1 + a2 - a3 + b0 - b1 - b3 >= -3',
                                    '- a0 - a1 + a2 - a3 - b0 + b1 - b2 >= -4',
                                    'a0 + a1 + a2 + b0 + b1 - b2 - b3 >= -1',
                                    'a0 - a1 - a2 + a3 + b0 + b2 + b3 >= -1',
                                    'a0 - a1 - a2 - b0 - b1 + b2 - b3 >= -4',
                                    'a0 - a1 - a2 - a3 + b1 - b2 - b3 >= -4',
                                    'a0 + a1 - a2 - a3 + b0 - b1 - b2 >= -3',
                                    '- a1 - a3 - b0 - b1 - b2 + b3 >= -4',
                                    'a0 - a1 - a2 + a3 - b0 + b2 - b3 >= -3',
                                    'a0 + a1 - a2 - a3 - b0 + b1 + b2 >= -2',
                                    'a0 - a2 + a3 - b0 + b1 - b2 + b3 >= -2',
                                    'a1 + a3 + b0 - b1 + b2 - b3 >= -1',
                                    '- a0 - a2 + b0 + b1 + b3 >= -1',
                                    'a1 + a2 + a3 + b1 - b2 - b3 >= -1',
                                    '- a0 - a1 - a2 - a3 + b1 + b2 >= -3',
                                    'a1 + a3 - b0 + b1 - b2 - b3 >= -2',
                                    '- a2 - a3 + b1 + b2 + b3 >= -1',
                                    '- a0 - a3 + b0 + b2 + b3 >= -1',
                                    '- a0 + b0 + b1 + b2 + b3 >= 0',
                                    '- a0 - a1 + b0 + b1 + b2 >= -1',
                                    '- a3 + b0 + b1 + b2 + b3 >= 0',
                                    '- a0 - a1 - a2 - a3 + b0 >= -3',
                                    '- a0 - a1 - a2 + b0 + b1 >= -2',
                                    'a2 + a3 - b1 - b2 + b3 >= -1',
                                    'a0 + a2 - b0 - b1 + b3 >= -1',
                                    'a0 + a1 + a2 + a3 - b1 >= 0']
        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 20
        # Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
        # Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
        self.deterministic_differential_forward_constraints = ['- a0 + b6 >= 0',
                                                        '- a1 + b6 >= 0',
                                                        '- a2 + b6 >= 0',
                                                        '- a3 + b6 >= 0',
                                                        '- a4 + b6 >= 0',
                                                        '- a5 + b6 >= 0',
                                                        '- a6 + b6 >= 0',
                                                        '- a7 + b6 >= 0',
                                                        '- b0 + b1 + b2 + b3 + b5 + b7 >= 0',
                                                        'b0 + b1 + b3 - b4 + b5 + b7 >= 0',
                                                        'b1 + b3 + b4 + b5 - b6 + b7 >= 0',
                                                        'a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + b1 - b2 + b3 + b5 + b7 >= 0',
                                                        '- a6 - a7 >= -1',
                                                        '- a4 - a5 >= -1',
                                                        '- a2 - a3 >= -1',
                                                        '- a0 - a1 >= -1',
                                                        '- b7 >= 0',
                                                        '- b5 >= 0',
                                                        '- b3 >= 0',
                                                        '- b1 >= 0']
        
        self.deterministic_differential_backward_constraints = ['- a0 + b6 >= 0',
                                                                '- a1 + b6 >= 0',
                                                                '- a2 + b6 >= 0',
                                                                '- a3 + b6 >= 0',
                                                                '- a4 + b6 >= 0',
                                                                '- a5 + b6 >= 0',
                                                                '- a6 + b6 >= 0',
                                                                '- a7 + b6 >= 0',
                                                                '- b0 + b1 + b2 + b3 + b5 + b7 >= 0',
                                                                'b0 + b1 + b3 - b4 + b5 + b7 >= 0',
                                                                'b1 + b3 + b4 + b5 - b6 + b7 >= 0',
                                                                'a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + b1 - b2 + b3 + b5 + b7 >= 0',
                                                                '- a6 - a7 >= -1',
                                                                '- a4 - a5 >= -1',
                                                                '- a2 - a3 >= -1',
                                                                '- a0 - a1 >= -1',
                                                                '- b7 >= 0',
                                                                '- b5 >= 0',
                                                                '- b3 >= 0',
                                                                '- b1 >= 0']

        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer
        # Number of constraints: 23
        # Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
        # Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
        self.deterministic_linear_backward_constraints = ['- a6 + b4 >= 0',
                                                        '- a7 + b4 >= 0',
                                                        'a5 - b5 >= 0',
                                                        '- b4 - b5 >= -1',
                                                        '- a0 + b6 >= 0',
                                                        '- a2 + b6 >= 0',
                                                        '- a4 + b6 >= 0',
                                                        '- a5 + b6 >= 0',
                                                        '- b4 + b6 >= 0',
                                                        '- a1 + a3 + b4 >= 0',
                                                        'a1 - b2 + b4 >= 0',
                                                        '- a3 + b4 + b5 >= 0',
                                                        '- b0 + b1 + b2 + b3 + b7 >= 0',
                                                        'b0 + b1 + b3 - b6 + b7 >= 0',
                                                        '- a1 - a3 - a5 + a6 + a7 + b5 >= -2',
                                                        'a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 - b2 >= 0',
                                                        '- a6 - a7 >= -1',
                                                        '- a4 - a5 >= -1',
                                                        '- a2 - a3 >= -1',
                                                        '- a0 - a1 >= -1',
                                                        '- b7 >= 0',
                                                        '- b3 >= 0',
                                                        '- b1 >= 0']
        # The following constraints encode the overlap between the two propagations
        # F: {-1, 0, 1}^2 --> {0, 1}
        # -1 -1   1
        # -1  0   0
        # -1  1   1
        #  0 -1   0
        #  0  0   0
        #  0  1   0
        #  1 -1   1
        #  1  0   0
        #  1  1   0
        
        # Generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer        
        # sage: binary_vectors =[(1, 0, 1, 0, 1), 
                                # (1, 0, 0, 0, 0), 
                                # (1, 0, 0, 1, 1),
                                # (0, 0, 1, 0, 0),
                                # (0, 0, 0, 0, 0),
                                # (0, 0, 0, 1, 0),
                                # (0, 1, 1, 0, 1),
                                # (0, 1, 0, 0, 0),
                                # (0, 1, 0, 1, 0)]
        # sage: cnf, milp = sb.encode_set_of_binary_vectors(binary_vectors, input_variables=['a0', 'a1', 'b0', 'b1', 'c0'])

        # Generateing and simplifying the MILP/SAT constraints ...
        # Time used to simplify the constraints: 0.00 seconds
        # Number of constraints: 8
        # Variables: a0||a1||b0||b1||c0; msb: a0
        self.overlap_detector_constraints = ['- a1 - b1 - c0 >= -2',
                                            'b0 + b1 - c0 >= 0',
                                            'a0 + a1 - c0 >= 0',
                                            '- a0 - b1 + c0 >= -1',
                                            '- a1 - b0 + c0 >= -1',
                                            '- a0 - b0 + c0 >= -1',
                                            '- b0 - b1 >= -1',
                                            '- a0 - a1 >= -1']
    
    def create_nibble_variables(self, s, r):
        '''
        Generate the nibble variables
        '''

        array = [f"{s}_{r}_{i}" for i in range(32)]
        self.used_variables.extend(array)
        return array
        
    def create_state_variables(self, s, r):
        '''
        Generate the state variables
        '''

        array = [f"{s}_{r}_{i}" for i in range(self.block_size)]
        self.used_variables.extend(array)
        return array
    
    def create_probability_variables(self, s, r):
        '''
        Generate the probability variables
        '''

        array = [[f"{s}_{r}_{nibble}_{0}", f"{s}_{r}_{nibble}_{1}"] for nibble in range(self.block_size//self.sbox_size)]
        self.used_variables.extend([item for sublist in array for item in sublist])
        return array

    def create_state_tuple_variables(self, s, r):
        '''
        Generate the state variables for deterministic propagation
        '''

        array = [[f"{s}_{r}_{i}_0", f"{s}_{r}_{i}_1"] for i in range(self.block_size)]
        self.used_variables.extend([item for sublist in array for item in sublist])
        for item in array:
            self.lp_constraints_for_joint_variables.append(f"{item[0]} + {item[1]} <= 1\n")
        return array
    
    def xor2(self, a, b, c):
        """
        Generate the constraints of a two-input XOR  (c = a xor b)
        - a - b - c >= -2
          a + b - c >= 0
          a - b + c >= 0
        - a + b + c >= 0
        """

        lp_contents = ""
        lp_contents += f"-1 {a} - {b} - {c} >= -2\n"
        lp_contents += f"{a} + {b} - {c} >= 0\n"
        lp_contents += f"{a} - {b} + {c} >= 0\n"
        lp_contents += f"-1 {a} + {b} + {c} >= 0\n"
        return lp_contents

    def xor3(self, a0, a1, a2, b):
        '''
        Generate the constraints of a three-input XOR  (b = a0 xor a1 xor a2)    
        b - a2 - a1 - a0 >= -2
        - b + a2 - a1 - a0 >= -2
        - b - a2 + a1 - a0 >= -2
        b + a2 + a1 - a0 >= 0
        - b - a2 - a1 + a0 >= -2
        b + a2 - a1 + a0 >= 0
        b - a2 + a1 + a0 >= 0
        - b + a2 + a1 + a0 >= 0
        The above inequalities are derived with QuineMcCluskey algorithm
        '''

        lp_contents = ""
        lp_contents += f"{b} - {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} + {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} - {a2} + {a1} - {a0} >= -2\n"
        lp_contents += f"{b} + {a2} + {a1} - {a0} >= 0\n"
        lp_contents += f"-1 {b} - {a2} - {a1} + {a0} >= -2\n"
        lp_contents += f"{b} + {a2} - {a1} + {a0} >= 0\n"
        lp_contents += f"{b} - {a2} + {a1} + {a0} >= 0\n"
        lp_contents += f"-1 {b} + {a2} + {a1} + {a0} >= 0\n"
        return lp_contents
    
    def xor3_deterministic(self, a0, a1, b0, b1, c0, c1, d0, d1):
        """
        Generate the constraints for deterministic propagations through a three-input XOR:
        d0||d1 = a0||a1 XOR b0||b1 XOR c0||c1
        Constraints are generated by Sbox Analyzer: https://github.com/hadipourh/sboxanalyzer    
        sage: import itertools
        sage: domain = {-1, 0, 1}
        sage: cartesian_product = list(itertools.product(domain, repeat=3))
        sage: valid_binary_vectors = []
        sage: translate = {-0: [0, 0], 1: [0, 1], -1: [1, 0]}
        ....: for item in cartesian_product:
        ....:     if -1 in item:
        ....:         output = -1
        ....:     else:
        ....:         output = sum(item) % 2
        ....:     valid_binary_vectors.append(tuple(translate[item[0]] + translate[item[1]] + translate[item[2]] + translate[output]))
        sage: cnf, milp = sb.encode_set_of_binary_vectors(valid_binary_vectors, input_variables=['a0', 'a1', 'b0', 'b1', 'c0', 'c1', 'd0', 'd1'])
        Generateing and simplifying the MILP/SAT constraints ...
        Time used to simplify the constraints: 0.00 seconds
        Number of constraints: 16
        Variables: a0||a1||b0||b1||c0||c1||d0||d1; msb: a0
        sage: pretty_print(milp)
        '- d0 - d1 >= -1',
        'a0 + a1 + b0 + b1 - c1 + d1 >= 0',
        'a0 + a1 - b1 + c0 + c1 + d1 >= 0',
        '- a1 + b0 + b1 + c0 + c1 + d1 >= 0',
        'a0 + b0 + c0 - d0 >= 0',
        '- a1 - b1 - c1 + d1 >= -2',
        'a1 + b1 + c1 - d1 >= 0',
        '- a1 - b1 + c1 - d1 >= -2',
        '- a1 + b1 - c1 - d1 >= -2',
        'a1 - b1 - c1 - d1 >= -2',
        '- b0 - b1 >= -1',
        '- c0 - c1 >= -1',
        '- a0 - a1 >= -1',
        '- c0 + d0 >= 0',
        '- b0 + d0 >= 0',
        '- a0 + d0 >= 0'
        """

        lp_contents = ""
        lp_contents += f"- {d0} - {d1} >= -1\n"
        lp_contents += f"{a0} + {a1} + {b0} + {b1} - {c1} + {d1} >= 0\n"
        lp_contents += f"{a0} + {a1} - {b1} + {c0} + {c1} + {d1} >= 0\n"
        lp_contents += f"- {a1} + {b0} + {b1} + {c0} + {c1} + {d1} >= 0\n"
        lp_contents += f"{a0} + {b0} + {c0} - {d0} >= 0\n"
        lp_contents += f"- {a1} - {b1} - {c1} + {d1} >= -2\n"
        lp_contents += f"{a1} + {b1} + {c1} - {d1} >= 0\n"
        lp_contents += f"- {a1} - {b1} + {c1} - {d1} >= -2\n"
        lp_contents += f"- {a1} + {b1} - {c1} - {d1} >= -2\n"
        lp_contents += f"{a1} - {b1} - {c1} - {d1} >= -2\n"
        lp_contents += f"- {b0} - {b1} >= -1\n"
        lp_contents += f"- {c0} - {c1} >= -1\n"
        lp_contents += f"- {a0} - {a1} >= -1\n"
        lp_contents += f"- {c0} + {d0} >= 0\n"
        lp_contents += f"- {b0} + {d0} >= 0\n"
        lp_contents += f"- {a0} + {d0} >= 0\n"
        return lp_contents
    
    def equality(self, x, y):
        '''
        Generate the MILP constraints modeling the equality of two bits
        '''

        lp_contents = f"{x} - {y} = 0\n"
        return lp_contents

    def permutation_layer_branch(self, x=None, y=None, branch=0, bitwise=False, double=False):
        '''
        Model the permutation layer for branch 0
        '''

        lp_contents = ""
        if bitwise and double:
            for i in range(32):
                for j in range(4):
                    lp_contents += self.equality(y[self.bit_permutation[branch][4*i + j]][0], x[4*i + j][0])
                    lp_contents += self.equality(y[self.bit_permutation[branch][4*i + j]][1], x[4*i + j][1])
        elif bitwise and not double:
            for i in range(32):
                for j in range(4):
                    lp_contents += self.equality(y[self.bit_permutation[branch][4*i + j]], x[4*i + j])
        elif not bitwise and double:
            for i in range(32):
                for j in range(4):
                    lp_contents += self.equality(y[4*self.nibble_permutation[branch][i] + j][0], x[4*i + j][0])
                    lp_contents += self.equality(y[4*self.nibble_permutation[branch][i] + j][1], x[4*i + j][1])
        else:
            for i in range(32):
                for j in range(4):
                    lp_contents += self.equality(y[4*self.nibble_permutation[branch][i] + j], x[4*i + j])
        return lp_contents

    def mix_columns(self, x, y, deterministic=False):
        '''
        Model the MixColumns in differential part
        '''

        lp_contents = ""
        if deterministic:
            for i in range(8):
                for j in range(4):
                    lp_contents += self.xor3_deterministic(
                        x[16*i + 4*1 + j][0],
                        x[16*i + 4*1 + j][1],
                        x[16*i + 4*2 + j][0],
                        x[16*i + 4*2 + j][1],
                        x[16*i + 4*3 + j][0],
                        x[16*i + 4*3 + j][1],
                        y[16*i + 4*0 + j][0],
                        y[16*i + 4*0 + j][1]
                    )
                    lp_contents += self.xor3_deterministic(
                        x[16*i + 4*0 + j][0],
                        x[16*i + 4*0 + j][1],                        
                        x[16*i + 4*2 + j][0],
                        x[16*i + 4*2 + j][1],
                        x[16*i + 4*3 + j][0],
                        x[16*i + 4*3 + j][1],
                        y[16*i + 4*1 + j][0],
                        y[16*i + 4*1 + j][1]
                    )
                    lp_contents += self.xor3_deterministic(
                        x[16*i + 4*0 + j][0],
                        x[16*i + 4*0 + j][1],
                        x[16*i + 4*1 + j][0],
                        x[16*i + 4*1 + j][1],
                        x[16*i + 4*3 + j][0],
                        x[16*i + 4*3 + j][1],
                        y[16*i + 4*2 + j][0],
                        y[16*i + 4*2 + j][1]
                    )
                    lp_contents += self.xor3_deterministic(
                        x[16*i + 4*0 + j][0],
                        x[16*i + 4*0 + j][1],
                        x[16*i + 4*1 + j][0],
                        x[16*i + 4*1 + j][1],
                        x[16*i + 4*2 + j][0],
                        x[16*i + 4*2 + j][1],
                        y[16*i + 4*3 + j][0],
                        y[16*i + 4*3 + j][1]
                    )
        else:
            for i in range(8):
                for j in range(4):
                    lp_contents += self.xor3(
                        x[16*i + 4*1 + j],
                        x[16*i + 4*2 + j],
                        x[16*i + 4*3 + j],
                        y[16*i + 4*0 + j]
                    )
                    lp_contents += self.xor3(
                        x[16*i + 4*0 + j],
                        x[16*i + 4*2 + j],
                        x[16*i + 4*3 + j],
                        y[16*i + 4*1 + j]
                    )
                    lp_contents += self.xor3(
                        x[16*i + 4*0 + j],
                        x[16*i + 4*1 + j],
                        x[16*i + 4*3 + j],
                        y[16*i + 4*2 + j]
                    )
                    lp_contents += self.xor3(
                        x[16*i + 4*0 + j],
                        x[16*i + 4*1 + j],
                        x[16*i + 4*2 + j],
                        y[16*i + 4*3 + j]
                    )
        return lp_contents
    
    def sbox_probabilistic(self, table=None, di=None, do=None, pr=None, activeness=None):
        """
        Generate constraints modeling the Sbox through the probabilistic part

        :param str[4] di: input difference
        :param str[4] do: output difference
        :param str[3] pr: probability of (di --> do) such that                    
        :return constraints encoding the DDT, or LAT of the Sbox
        :rtype str:
        """

        if table == 'ddt':
            inequalities = self.ddt_constraints
        elif table == 'lat':
            inequalities = self.square_lat_constraints
        elif table == 'starddt':
            inequalities = self.star_ddt_constraints
        elif table == 'starlat':
            inequalities = self.star_lat_constraints
        else: 
            raise ValueError("Invalid table type!")

        constraints = ""
        for ineq in inequalities:
            temp = ineq
            for i in range(4):
                temp = temp.replace(f"a{i}", di[i])
            for i in range(4):
                temp = temp.replace(f"b{i}", do[i])
            if table not in ['starddt', 'starlat']:
                for i in range(2):
                    temp = temp.replace(f"p{i}", pr[i])
            constraints += temp + "\n"
        if activeness != None:
            constraints += f"{activeness}" + " - " + " - ".join(di) + " <= 0\n"
            for i in range(4):
                constraints += f"{activeness} - {di[i]} >= 0\n"
        return constraints
    
    def sbox_deterministic(self, table=None, x=None, y=None, activeness=None):
        """
        Generate constraints modeling the Sbox through the deterministic part
        """

        if table == 'ddt':
            inequalities = self.deterministic_differential_forward_constraints
        elif table == 'lat':
            inequalities = self.deterministic_linear_backward_constraints
        elif table == 'backwardddt':
            inequalities = self.deterministic_differential_backward_constraints
        else:
            raise ValueError("Invalid table type!")
        
        constraints = ""
        for ineq in inequalities:
            temp = ineq
            for i in range(4):
                temp = temp.replace(f"a{2*i}", x[i][0])
                temp = temp.replace(f"a{2*i + 1}", x[i][1])
                temp = temp.replace(f"b{2*i}", y[i][0])
                temp = temp.replace(f"b{2*i + 1}", y[i][1])
            constraints += temp + "\n"
        if activeness != None:
            constraints += " + ".join([bitvars for vartuple in x for bitvars in vartuple]) + " - " + f"{activeness}" + " >= 0\n"
            for i in range(4):
                for j in range(2):
                    constraints += f"{activeness} - {x[i][j]} >= 0\n"
        return constraints
    
    def overlap_detector(self, x, y, z):
        '''
        Model the overlap detector
        each item in x has two bits
        each item in y has two bits
        each item in z has one bit
        '''

        lp_constraints = ""
        for i in range(self.block_size):
            inequalities = self.overlap_detector_constraints 
            for ineq in inequalities:
                temp = ineq
                temp = temp.replace('a0', x[i][0])
                temp = temp.replace('a1', x[i][1])
                temp = temp.replace('b0', y[i][0])
                temp = temp.replace('b1', y[i][1])
                temp = temp.replace('c0', z[i])
                lp_constraints += temp + '\n'
        return lp_constraints

    def model_eb(self):
        '''
        Generate the model for EB
        '''
        
        lp_contents = ""
        for r in range(self.RB):
            # model branch 0
            xb0 = self.create_state_tuple_variables('x0', r)
            xb0_next = self.create_state_tuple_variables('x0', r + 1)
            yb0 = self.create_state_tuple_variables('y0', r)
            zb0 = self.create_state_tuple_variables('z0', r)
            active0 = self.create_nibble_variables('active0', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='backwardddt', 
                                                       x=yb0[4*nibble:4*nibble + 4], 
                                                       y=xb0[4*nibble:4*nibble + 4],
                                                       activeness=active0[nibble])
            # model the permutation layer
            if self.offset + r < 4:
                lp_contents += self.permutation_layer_branch(yb0, zb0, branch=0, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(yb0, zb0, branch=0, bitwise=False, double=True)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(xb0_next, zb0, deterministic=True)

            # model branch 1
            xb1 = self.create_state_tuple_variables('x1', r)
            xb1_next = self.create_state_tuple_variables('x1', r + 1)
            yb1 = self.create_state_tuple_variables('y1', r)
            zb1 = self.create_state_tuple_variables('z1', r)
            active1 = self.create_nibble_variables('active1', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='backwardddt', 
                                                       x=yb1[4*nibble:4*nibble + 4], 
                                                       y=xb1[4*nibble:4*nibble + 4],
                                                       activeness=active1[nibble])
            # model the permutation layer
            if self.offset + r < 4:
                lp_contents += self.permutation_layer_branch(yb1, zb1, branch=1, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(yb1, zb1, branch=1, bitwise=False, double=True)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(xb1_next, zb1, deterministic=True)

        # the activeness pattern at the input of key recovery should be the same
        ab0 = self.create_nibble_variables('active0', 0)
        ab1 = self.create_nibble_variables('active1', 0)
        for i in range(32):
            lp_contents += self.equality(ab0[i], ab1[i])
        
        # some artificial constraints to control the number of active bits
        lp_contents += f"{ab0[27]} + {ab0[29]} <= 1\n"
        temp = " + ".join(ab0)
        lp_contents += f"{temp} <= 2\n"
        
        # control the number of common active master key bits
        key_overlap = self.create_state_variables('key_overlap', 0)
        k0 = self.create_state_variables('k0', 0)
        k1 = self.create_state_variables('k1', 0)
        xb0 = self.create_state_tuple_variables('x0', 0)
        xb1 = self.create_state_tuple_variables('x1', 0)
        for i in range(self.block_size):
            lp_contents += f"{xb0[i][0]} + {xb0[i][1]} - {k0[self.inv_kperm_at_round[self.offset][0][i]]} = 0\n"
            lp_contents += f"{xb1[i][0]} + {xb1[i][1]} - {k1[self.inv_kperm_at_round[self.offset][1][i]]} = 0\n"
            lp_contents += self.xor2(k0[i], k1[i], key_overlap[i])
        # lp_contents += "{}".format(" + ".join(key_overlap) + " = 0\n")
        return lp_contents
    
    def model_eu(self):
        """
        Generate the model for EU
        """

        # connect the output of EB to the input of EU
        lp_contents = ""
        xb0 = self.create_state_tuple_variables('x0', self.RB)
        xb1 = self.create_state_tuple_variables('x1', self.RB)
        xu0 = self.create_state_variables('xu0', 0)
        xu1 = self.create_state_variables('xu1', 0)
        for bit in range(self.block_size):
            lp_contents += self.equality(xb0[bit][1], xu0[bit])
            lp_contents += self.equality(xb1[bit][1], xu1[bit])
            lp_contents += self.equality(xb0[bit][0], 0)
            lp_contents += self.equality(xb1[bit][0], 0)
        for r in range(self.RU):
            xu0 = self.create_state_variables('xu0', r)
            yu0 = self.create_state_variables('yu0', r)
            zu0 = self.create_state_variables('zu0', r)
            xu0_next = self.create_state_variables('xu0', r + 1)
            pu0 = self.create_probability_variables('pu0', r)
            
            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_probabilistic(table='ddt',
                                                         di=xu0[4*nibble:4*nibble + 4],
                                                         do=yu0[4*nibble:4*nibble + 4],
                                                         pr=pu0[nibble],
                                                         activeness=None)
            
            # model the permutation layer
            if self.offset + self.RB + r < 4:
                lp_contents += self.permutation_layer_branch(yu0, zu0, branch=0, bitwise=True, double=False)
            else:
                lp_contents += self.permutation_layer_branch(yu0, zu0, branch=0, bitwise=False, double=False)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zu0, xu0_next, deterministic=False)

            # model the branch 1
            xu1 = self.create_state_variables('xu1', r)
            yu1 = self.create_state_variables('yu1', r)
            zu1 = self.create_state_variables('zu1', r)
            xu1_next = self.create_state_variables('xu1', r + 1)
            pu1 = self.create_probability_variables('pu1', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_probabilistic(table='ddt',
                                                         di=xu1[4*nibble:4*nibble + 4],
                                                         do=yu1[4*nibble:4*nibble + 4],
                                                         pr=pu1[nibble],
                                                         activeness=None)
            
            # model the permutation layer
            if self.offset + self.RB + r < 4:
                lp_contents += self.permutation_layer_branch(yu1, zu1, branch=1, bitwise=True, double=False)
            else:
                lp_contents += self.permutation_layer_branch(yu1, zu1, branch=1, bitwise=False, double=False)

            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zu1, xu1_next, deterministic=False)

        # exclude all-zero input difference
        xu0 = self.create_state_variables('xu0', self.RU)
        xu1 = self.create_state_variables('xu1', self.RU)
        lp_contents += "{}".format(" + ".join(xu0) + " >= 1\n")
        lp_contents += "{}".format(" + ".join(xu1) + " >= 1\n")

        return lp_contents
    
    def model_em(self):
        """
        Generate the model for EM
        """

        lp_contents = ""
        # link the output of EU to the input of EM
        xu0 = self.create_state_variables('xu0', self.RU)
        xu1 = self.create_state_variables('xu1', self.RU)
        xmu0 = self.create_state_tuple_variables('xmu0', 0)
        xmu1 = self.create_state_tuple_variables('xmu1', 0)
        for bit in range(self.block_size):
            lp_contents += self.equality(xu0[bit], xmu0[bit][1])
            lp_contents += self.equality(xu1[bit], xmu1[bit][1])
            lp_contents += self.equality(xmu0[bit][0], 0)
            lp_contents += self.equality(xmu1[bit][0], 0)
        for r in range(self.RM):
            # model the branch 0
            xmu0 = self.create_state_tuple_variables('xmu0', r)
            ymu0 = self.create_state_tuple_variables('ymu0', r)
            zmu0 = self.create_state_tuple_variables('zmu0', r)
            xmu0_next = self.create_state_tuple_variables('xmu0', r + 1)
            active_mu0 = self.create_nibble_variables('active_mu0', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='ddt',
                                                       x=xmu0[4*nibble:4*nibble + 4],
                                                       y=ymu0[4*nibble:4*nibble + 4],
                                                       activeness=active_mu0[nibble])
            
            # model the permutation layer
            if self.offset + self.RB + self.RU + r < 4:
                lp_contents += self.permutation_layer_branch(ymu0, zmu0, branch=0, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(ymu0, zmu0, branch=0, bitwise=False, double=True)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zmu0, xmu0_next, deterministic=True)

        for r in range(self.RM):
            # model the branch 1
            xmu1 = self.create_state_tuple_variables('xmu1', r)
            ymu1 = self.create_state_tuple_variables('ymu1', r)
            zmu1 = self.create_state_tuple_variables('zmu1', r)
            xmu1_next = self.create_state_tuple_variables('xmu1', r + 1)
            active_mu1 = self.create_nibble_variables('active_mu1', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='ddt',
                                                       x=xmu1[4*nibble:4*nibble + 4],
                                                       y=ymu1[4*nibble:4*nibble + 4],
                                                       activeness=active_mu1[nibble])
                
            # model the permutation layer
            if self.offset + self.RB + self.RU + r < 4:
                lp_contents += self.permutation_layer_branch(ymu1, zmu1, branch=1, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(ymu1, zmu1, branch=1, bitwise=False, double=True)

            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zmu1, xmu1_next, deterministic=True)

        for r in range(self.RM):
            # model the branch 0
            xml0 = self.create_state_tuple_variables('xml0', r)
            yml0 = self.create_state_tuple_variables('yml0', r)
            zml0 = self.create_state_tuple_variables('zml0', r)
            xml0_next = self.create_state_tuple_variables('xml0', r + 1)
            active_ml0 = self.create_nibble_variables('active_ml0', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='lat',
                                                       x=yml0[4*nibble:4*nibble + 4],
                                                       y=xml0[4*nibble:4*nibble + 4],
                                                       activeness=active_ml0[nibble])
            
            # model the permutation layer
            if self.offset + self.RB + self.RU + r < 4:
                lp_contents += self.permutation_layer_branch(yml0, zml0, branch=0, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(yml0, zml0, branch=0, bitwise=False, double=True)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(xml0_next, zml0, deterministic=True)
            
        for r in range(self.RM):
            # model the branch 1
            xml1 = self.create_state_tuple_variables('xml1', r)
            yml1 = self.create_state_tuple_variables('yml1', r)
            zml1 = self.create_state_tuple_variables('zml1', r)
            xml1_next = self.create_state_tuple_variables('xml1', r + 1)
            active_ml1 = self.create_nibble_variables('active_ml1', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_deterministic(table='lat',
                                                       x=yml1[4*nibble:4*nibble + 4],
                                                       y=xml1[4*nibble:4*nibble + 4],
                                                       activeness=active_ml1[nibble])
                
            # model the permutation layer
            if self.offset + self.RB + self.RU + r < 4:
                lp_contents += self.permutation_layer_branch(yml1, zml1, branch=1, bitwise=True, double=True)
            else:
                lp_contents += self.permutation_layer_branch(yml1, zml1, branch=1, bitwise=False, double=True)

            # model the mix-coulumns layer
            lp_contents += self.mix_columns(xml1_next, zml1, deterministic=True)
        
        for r in range(self.RM):
            xmu0 = self.create_state_tuple_variables('xmu0', r)
            xml0 = self.create_state_tuple_variables('xml0', r)
            overlap_m0 = self.create_state_variables('overlap_m0', r)
            lp_contents += self.overlap_detector(xmu0, xml0, overlap_m0)
            xmu1 = self.create_state_tuple_variables('xmu1', r)
            xml1 = self.create_state_tuple_variables('xml1', r)
            overlap_m1 = self.create_state_variables('overlap_m1', r)
            lp_contents += self.overlap_detector(xmu1, xml1, overlap_m1)
        return lp_contents
    
    def model_el(self):
        """
        Generate the model for EL
        """

        lp_contents = ""
        # link the output of EM to the input of EL
        xml0 = self.create_state_tuple_variables('xml0', self.RM)
        xml1 = self.create_state_tuple_variables('xml1', self.RM)
        xl0 = self.create_state_variables('xl0', 0)
        xl1 = self.create_state_variables('xl1', 0)
        for bit in range(self.block_size):
            lp_contents += self.equality(xml0[bit][1], xl0[bit])
            lp_contents += self.equality(xml1[bit][1], xl1[bit])
            lp_contents += self.equality(xml0[bit][0], 0)
            lp_contents += self.equality(xml1[bit][0], 0)
        for r in range(self.RL):
            xl0 = self.create_state_variables('xl0', r)
            yl0 = self.create_state_variables('yl0', r)
            zl0 = self.create_state_variables('zl0', r)
            xl0_next = self.create_state_variables('xl0', r + 1)
            pl0 = self.create_probability_variables('pl0', r)
            
            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_probabilistic(table='lat',
                                                       di=xl0[4*nibble:4*nibble + 4],
                                                       do=yl0[4*nibble:4*nibble + 4],
                                                       pr=pl0[nibble],
                                                       activeness=None)
            
            # model the permutation layer
            if self.offset + self.RB + self.RU + self.RM + r < 4:
                lp_contents += self.permutation_layer_branch(yl0, zl0, branch=0, bitwise=True, double=False)
            else:
                lp_contents += self.permutation_layer_branch(yl0, zl0, branch=0, bitwise=False, double=False)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zl0, xl0_next, deterministic=False)

            # model the branch 1
            xl1 = self.create_state_variables('xl1', r)
            yl1 = self.create_state_variables('yl1', r)
            zl1 = self.create_state_variables('zl1', r)
            xl1_next = self.create_state_variables('xl1', r + 1)
            pl1 = self.create_probability_variables('pl1', r)

            # model the sbox layer
            for nibble in range(32):
                lp_contents += self.sbox_probabilistic(table='ddt',
                                                       di=xl1[4*nibble:4*nibble + 4],
                                                       do=yl1[4*nibble:4*nibble + 4],
                                                       pr=pl1[nibble],
                                                       activeness=None)
            
            # model the permutation layer
            if self.offset + self.RB + self.RU + self.RM + r < 4:
                lp_contents += self.permutation_layer_branch(yl1, zl1, branch=1, bitwise=True, double=False)
            else:
                lp_contents += self.permutation_layer_branch(yl1, zl1, branch=1, bitwise=False, double=False)
            
            # model the mix-coulumns layer
            lp_contents += self.mix_columns(zl1, xl1_next, deterministic=False)
           
        # exclude all-zero output mask
        xl0 = self.create_state_variables('xl0', self.RL)
        xl1 = self.create_state_variables('xl1', self.RL)
        lp_contents += "{}".format(" + ".join(xl0) + " >= 1\n")
        lp_contents += "{}".format(" + ".join(xl1) + " >= 1\n")
        # output mask of branch 0 and 1 should be the same
        for i in range(self.block_size):
            lp_contents += self.equality(xl0[i], xl1[i])
        return lp_contents

    def create_objective_function(self):
        """
        Create the objective function
        """
        
        pu0 = []
        pu1 = []
        for r in range(self.RU):
            pu0vars = self.create_probability_variables('pu0', r)
            pu1vars = self.create_probability_variables('pu1', r)
            for nibble in range(32):
                pu0.append(f"{3*self.WU} {pu0vars[nibble][0]} + {2*self.WU} {pu0vars[nibble][1]}")
                pu1.append(f"{3*self.WU} {pu1vars[nibble][0]} + {2*self.WU} {pu1vars[nibble][1]}")
        pu = " + ".join(pu0 + pu1)
        
        overlaps_0 = []
        overlaps_1 = []
        for r in range(self.RM):
            temp_0 = self.create_state_variables('overlap_m0', r)
            temp_1 = self.create_state_variables('overlap_m1', r)
            for i in range(self.block_size):
                overlaps_0.append(f"{self.WM} {temp_0[i]}")
                overlaps_1.append(f"{self.WM} {temp_1[i]}")
        overlap = " + ".join(overlaps_0 + overlaps_1)

        pl0 = []
        pl1 = []
        for r in range(self.RL):
            pl0vars = self.create_probability_variables('pl0', r)
            pl1vars = self.create_probability_variables('pl1', r)
            for nibble in range(32):
                pl0.append(f"{4*self.WL} {pl0vars[nibble][0]} + {2*self.WL} {pl0vars[nibble][1]}")
                pl1.append(f"{4*self.WL} {pl1vars[nibble][0]} + {2*self.WL} {pl1vars[nibble][1]}")
        cl = " + ".join(pl0 + pl1)
    
        correlation = ""
        if self.RU > 0 and (self.RM + self.RL) == 0:
            correlation += f"{pu}"
        elif self.RU == 0 and self.RM > 0 and self.RL == 0:
            correlation += f"{overlap}"
        elif self.RU == 0 and self.RM == 0 and self.RL > 0:
            correlation += f"{cl}"
        elif self.RU > 0 and self.RM > 0 and self.RL == 0:
            correlation += f"{pu} + {overlap}"
        elif self.RU == 0 and self.RM > 0 and self.RL > 0:
            correlation += f"{overlap} + {cl}"
        elif self.RU > 0 and self.RM == 0 and self.RL > 0:
            correlation += f"{pu} + {cl}"
        elif self.RU > 0 and self.RM > 0 and self.RL > 0:
            correlation += f"{pu} + {overlap} + {cl}"
        else:
            raise ValueError("Invalid combination of rounds!")   
        objective_function = f"minimize\n{correlation}\n"
        return objective_function

    def constraint_joint_varaibles_to_exclude_1_1(self):
        """
        Generate constraints to exclude the (x0||x1) = (1||1)
        """
        lp_constraints = "".join(self.lp_constraints_for_joint_variables)
        return lp_constraints
    
    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():            
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            if len(var) == 2:
                state_vars = [f"{var[0]}_{var[1]}_{i}" for i in range(64)]
                for i in range(64):
                    if val[i] != "?":
                        lp_contents += f"{state_vars[i]} = {val[i]}\n"
            if len(var) == 3:
                state_var = f"{var[0]}_{var[1]}_{var[2]}"                
                if val != "?":
                    lp_contents += f"{state_var} = {val}\n"
        return lp_contents

    def declare_variables_type(self):
        '''
        Specifying variables' type in the LP file
        '''
        
        lp_contents = 'binary\n'
        self.used_variables = list(set(self.used_variables))
        for var in self.used_variables:
            lp_contents += var + '\n'            
        lp_contents += "end\n"
        return lp_contents

    def make_model(self):
        '''
        Make the model
        '''
        
        lp_contents = ""
        print('Generating the MILP model ...')
        lp_contents = self.create_objective_function()
        lp_contents += "subject to\n"
        lp_contents += self.model_eb()
        lp_contents += self.model_eu()
        lp_contents += self.model_em()
        lp_contents += self.model_el()
        lp_contents += self.constraint_joint_varaibles_to_exclude_1_1()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_variables_type()
        if os.path.exists(self.model_filename):
            os.remove(self.model_filename)
        with open(self.model_filename, 'w') as fileobj:
            fileobj.write(lp_contents)
        print(f"MILP model was written into {self.model_filename}\n")  

    def find_distinguisher(self):
        '''
        Find the best differential-linear trail under the given constraints, 
        e.g., satisfying an activeness pattern
        '''
        
        status = False
        if self.time_limit != -1:
            self.model.Params.TIME_LIMIT = self.time_limit
        time_start = time.time()
        self.model.Params.Threads = self.number_of_parallel_threads
        #self.model.Params.PreSolve = 0
        self.model.Params.OutputFlag = True
        self.model.optimize()
        self.attack_summary = "Attack summary:\n"
        if (self.model.Status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.INTERRUPTED, GRB.SOLUTION_LIMIT]):
            status = True
            self.attack_summary_0, self.upper_trail_0, self.lower_trail_0 = self.parse_solution(bn=0)
            self.attack_summary_1, self.upper_trail_1, self.lower_trail_1 = self.parse_solution(bn=1)
            if self.compute_differential_effect_using_gurobi:
                diff_effect_0 = self.compute_differential_effect(0)
                diff_effect_1 = self.compute_differential_effect(1)
                self.attack_summary_0 += "Differential effect over EU: {:.2f}".format(diff_effect_0) + "\n"
                self.attack_summary_1 += "Differential effect over EU: {:.2f}".format(diff_effect_1) + "\n"
            self.attack_summary += self.attack_summary_0 + "\n" + "#"*100 + "\n" + self.attack_summary_1
            print(self.attack_summary)
            draw_0 = DrawDL(self, output_file_name=self.output_file_name[0], bn=0)
            draw_1 = DrawDL(self, output_file_name=self.output_file_name[1], bn=1)
            draw_0.generate_distinguisher_shape()
            draw_1.generate_distinguisher_shape()
            if self.compute_differential_effect_using_gurobi:
                print("-Log2(DE0)            ~= \t{:.2f}".format(-1*diff_effect_0))
                print("-Log2(DE1)            ~= \t{:.2f}".format(-1*diff_effect_1))             
        elif self.model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
            status = self.model.Status
        else:
            print("Unknown error!")
        time_end = time.time()
        print("Time used = {:0.02f}".format(time_end - time_start))
        return status

    def compute_differential_effect(self, bn):
        """
        Compute the differential effect over EU
        """

        # Load default values
        params = {"nrounds" : self.RU,
                  "branchtype" : bn,
                  "offset" : self.RB,
                  "mode" : 2,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 2,
                  "fixedVariables" : {}}
        input_diff = hex(int("".join(list(map(str, self.result[f"xu{bn}"][0]))), 2))[2:].zfill(32)
        params["fixedVariables"][f"x_0"] = input_diff
        outputdiff = hex(int("".join(list(map(str, self.result[f"xu{bn}"][self.RU]))), 2))[2:].zfill(32)
        params["fixedVariables"][f"x_{self.RU}"] = outputdiff
        UDiff = Diff(params)
        UDiff.make_model()
        output = UDiff.solve()
        return output

    def solve(self, mip_focus=None):        
        self.model = read(self.model_filename)
        if mip_focus != None:
            self.model.Params.MIPFocus = mip_focus
        status = False
        status = self.find_distinguisher()
        os.remove(self.model_filename)
        return status    

    def parse_solution(self, bn=0):
        """
        Parse the solution and print the distinguisher's specifications
        """
        
        upper_trail = {"x": [[0 for _ in range(self.block_size)] for _ in range(self.RU + self.RM + 1)],
                       "y": [[0 for _ in range(self.block_size)] for _ in range(self.RU + self.RM)],
                       "z": [[0 for _ in range(self.block_size)] for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):
            for bit in range(self.block_size):
                upper_trail["x"][r][bit] = str(int(self.model.getVarByName(f"xu{bn}_{r}_{bit}").Xn))
                upper_trail["y"][r][bit] = str(int(self.model.getVarByName(f"yu{bn}_{r}_{bit}").Xn))
                upper_trail["z"][r][bit] = str(int(self.model.getVarByName(f"zu{bn}_{r}_{bit}").Xn))
        if self.RM > 0:
            for r in range(self.RU, self.RU + self.RM + 1):
                for bit in range(self.block_size):
                    upper_trail["x"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"xmu{bn}_{r - self.RU}_{bit}_0").Xn), int(self.model.getVarByName(f"xmu{bn}_{r - self.RU}_{bit}_1").Xn))]
                    if r < self.RU + self.RM:
                        upper_trail["y"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"ymu{bn}_{r - self.RU}_{bit}_0").Xn), int(self.model.getVarByName(f"ymu{bn}_{r - self.RU}_{bit}_1").Xn))]
                        upper_trail["z"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"zmu{bn}_{r - self.RU}_{bit}_0").Xn), int(self.model.getVarByName(f"zmu{bn}_{r - self.RU}_{bit}_1").Xn))]
        lower_trail = {"x": [[0 for _ in range(self.block_size)] for _ in range(self.RM + self.RL + 1)],
                       "y": [[0 for _ in range(self.block_size)] for _ in range(self.RM + self.RL)],
                       "z": [[0 for _ in range(self.block_size)] for _ in range(self.RM + self.RL)]}
        if self.RM > 0:
            for r in range(self.RM + 1):
                for bit in range(self.block_size):
                    lower_trail["x"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"xml{bn}_{r}_{bit}_0").Xn), int(self.model.getVarByName(f"xml{bn}_{r}_{bit}_1").Xn))]
                    if r < self.RM:
                        lower_trail["y"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"yml{bn}_{r}_{bit}_0").Xn), int(self.model.getVarByName(f"yml{bn}_{r}_{bit}_1").Xn))]
                        lower_trail["z"][r][bit] = self.translate_joint_vars_to_single_vars[(int(self.model.getVarByName(f"zml{bn}_{r}_{bit}_0").Xn), int(self.model.getVarByName(f"zml{bn}_{r}_{bit}_1").Xn))]
        if self.RL > 0:
            for r in range(self.RM, self.RM + self.RL + 1):
                for bit in range(self.block_size):
                    lower_trail["x"][r][bit] = str(int(self.model.getVarByName(f"xl{bn}_{r - self.RM}_{bit}").Xn))
                    if r < self.RM + self.RL:
                        lower_trail["y"][r][bit] = str(int(self.model.getVarByName(f"yl{bn}_{r - self.RM}_{bit}").Xn))
                        lower_trail["z"][r][bit] = str(int(self.model.getVarByName(f"zl{bn}_{r - self.RM}_{bit}").Xn))

        convert_to_fixed_vector = lambda x: [item if item != -1 else 0 for item in x]
        input_diff = f"char inputdiff[] = \"" + hex(int("".join(list(map(str, convert_to_fixed_vector(upper_trail["x"][0])))), 2))[2:].zfill(32) + "\";\n"
        input_diff_middle = f"char inputdiff[] = \"" + hex(int("".join(list(map(str, convert_to_fixed_vector(upper_trail["x"][self.RU])))), 2))[2:].zfill(32) + "\";\n"
        output_mask_middle = f"char outputmask[] = \"" + hex(int("".join(list(map(str, convert_to_fixed_vector(lower_trail["x"][self.RM])))), 2))[2:].zfill(32) + "\";\n"
        output_mask = f"char outputmask[] = \"" + hex(int("".join(list(map(str, convert_to_fixed_vector(lower_trail["x"][self.RM + self.RL])))), 2))[2:].zfill(32) + "\";\n"
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: Offset: {self.offset}, RB: {self.RB}, RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{input_diff_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask. middle: \n{output_mask_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}"
        attack_summary += "#"*50 + "\n"
        attack_summary += "PU{}:  {}\n".format(bn, sum([3*int(self.model.getVarByName(f"pu{bn}_{r}_{nibble}_0").Xn) + 2*int(self.model.getVarByName(f"pu{bn}_{r}_{nibble}_1").Xn) for r in range(self.RU) for nibble in range(32)]))
        attack_summary += "CM{}:  {}\n".format(bn, sum([int(self.model.getVarByName(f"overlap_m{bn}_{r}_{bit}").Xn) for r in range(self.RM) for bit in range(128)]))
        attack_summary += "CL{}:  {}\n".format(bn, sum([4*int(self.model.getVarByName(f"pl{bn}_{r}_{nibble}_0").Xn) + 2*int(self.model.getVarByName(f"pl{bn}_{r}_{nibble}_1").Xn) for r in range(self.RL) for nibble in range(32)]))
        attack_summary += "#"*50 + "\n"
        # print the upper trail
        attack_summary += "Upper trail:\n"
        for r in range(self.RU + self.RM + 1):
            attack_summary += f"Round {r}:\n"
            attack_summary += "x{:02d} = ".format(r) + "".join(list(map(str, upper_trail["x"][r]))).replace("-1", "*") + "\n"
            if r < self.RU + self.RM:
                attack_summary += "y{:02d} = ".format(r) + "".join(list(map(str, upper_trail["y"][r]))).replace("-1", "*") + "\n"
                attack_summary += "z{:02d} = ".format(r) + "".join(list(map(str, upper_trail["z"][r]))).replace("-1", "*") + "\n"
            attack_summary += "#"*50 + "\n\n"
        # print the lower trail
        attack_summary += "Lower trail:\n"
        for r in range(self.RM + self.RL + 1):
            attack_summary += "Round {:02d}:\n".format(r)
            attack_summary += "x{:02d} = ".format(r) + "".join(list(map(str, lower_trail["x"][r]))).replace("-1", "*") + "\n"
            if r < self.RM + self.RL:
                attack_summary += "y{:02d} = ".format(r) + "".join(list(map(str, lower_trail["y"][r]))).replace("-1", "*") + "\n"
                attack_summary += "z{:02d} = ".format(r) + "".join(list(map(str, lower_trail["z"][r]))).replace("-1", "*") + "\n"
            attack_summary += "#"*50 + "\n\n"
        attack_summary += "#"*50 + "\n"
        if self.RB > 0:
            attack_summary += "Involved key bits in branch {}: ".format(bn) +\
                            ", ".join(list(map(str, [i for i in range(128) if int(self.model.getVarByName(f"k{bn}_0_{i}").Xn) == 1]))) + "\n"
        return attack_summary, upper_trail, lower_trail
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#  _   _                    ___         _                __                   
# | | | | ___   ___  _ __  |_ _| _ __  | |_  ___  _ __  / _|  __ _   ___  ___ 
# | | | |/ __| / _ \| '__|  | | | '_ \ | __|/ _ \| '__|| |_  / _` | / __|/ _ \
# | |_| |\__ \|  __/| |     | | | | | || |_|  __/| |   |  _|| (_| || (__|  __/
#  \___/ |___/ \___||_|    |___||_| |_| \__|\___||_|   |_|   \__,_| \___|\___|
                                                                            
def loadparameters(args):
    '''
    Extract parameters from the argument list and input file
    '''

    # Load default values
    params = {"offset": 0,
            "RB": 0,
            "RU": 0,
            "RM": 2,
            "RL": 0,
            "WU": 1,
            "WM": 1,
            "np": 8,
            "tl": -1,
            "solver": "ortools",
            "output1": "output1.tex",
            "output2": "output2.tex",
            "fixedVariables": {}}

    # Override parameters if they are set on command line
    if args.RB is not None:
        params["RB"] = args.RB
    if args.RU is not None:
        params["RU"] = args.RU
    if args.RM is not None:
        params["RM"] = args.RM
    if args.RL is not None:
        params["RL"] = args.RL
    if args.WU is not None:
        params["WU"] = args.WU
    if args.WM is not None:
        params["WM"] = args.WM
    if args.WL is not None:
        params["WL"] = args.WL
    if args.np is not None:
        params["np"] = args.np
    if args.timelimit is not None:
        params["timelimit"] = args.timelimit
    if args.output1 is not None:
        params["output1"] = args.output1
    if args.output2 is not None:
        params["output2"] = args.output2
    if args.differential_effect is not None:
        params["differential_effect"] = args.differential_effect
    if args.offset is not None:
        params["offset"] = args.offset

    return params

def main():
    '''
    Parse the arguments and start the request functionality with the provided
    parameters.
    '''
    
    parser = ArgumentParser(description="This tool finds a nearly differential-linear distinguishers using Gurobi",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-offset", type=int, default=2, help="Number of rounds skipped from the top part of the design")    
    parser.add_argument("-RB", type=int, default=1, help="Offset for the starting round of the distinguisher")
    parser.add_argument("-RU", type=int, default=1, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=4, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=0, help="Number of rounds for EL")

    parser.add_argument("-WU", type=int, default=6, help="Scale of the probability transition over EU")
    parser.add_argument("-WM", type=int, default=2, help="Scale of the correlation of DL distinguishers over EM ")
    parser.add_argument("-WL", type=int, default=1, help="Scale of the squared correaltion of linear approximation over EL")

    parser.add_argument("-de", "--differential_effect", type=int, default=0, help="Compute differential effect over EU using Gurobi", choices=[0, 1])


    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=20, help="Time limit in seconds")  
    parser.add_argument("-o1", "--output1", default="output1.tex", type=str, help="Output file name 1")
    parser.add_argument("-o2", "--output2", default="output2.tex", type=str, help="Output file name 2")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DL(params)
    dld.make_model()
    dld.solve()

if __name__ == "__main__":
    main()
