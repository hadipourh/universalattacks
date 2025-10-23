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

import logging
from pathlib import Path
from random import randint
logging.basicConfig(filename="minizinc-python.log", level=logging.DEBUG)
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter
from drawprf import DrawDL
import subprocess

# Check if gurobi_cl command is recognized
import subprocess
try:
    subprocess.run(['gurobi_cl', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    gurobi_available = True
except FileNotFoundError:
    gurobi_available = False

if gurobi_available:
    print("Gurobi is available")
    from diff import Diff

# Check if ./difflin is present in the current directory
if not Path("./difflin").is_file():
    raise FileNotFoundError("The file './difflin' is not present in the current directory. Please ensure it exists.")

class DiffLin:
    DL_counter = 0

    def __init__(self, param) -> None:
        DiffLin.DL_counter += 1
        self.id = DiffLin.DL_counter
        self.name = "DiffLin" + str(self.id)
        self.type = "DiffLin"

        self.RB = param["RB"]
        self.RU = param["RU"]
        self.RM = param["RM"]
        self.RL = param["RL"]
        self.RMU = param["RMU"]
        self.RML = param["RML"]
        self.WU = param["WU"]
        self.WM = param["WM"]
        self.WL = param["WL"]
        self.degree = param["degree"]
        self.offset = param["offset"]
        self.compute_differential_effect_using_gurobi = param["differential_effect"]

        self.RD = self.RU + self.RM + self.RL
        self.cp_solver_name = param["solver"]
        self.cp_solver = minizinc.Solver.lookup(self.cp_solver_name)
        self.time_limit = param["timelimit"]
        self.num_of_threads = param["np"]
        self.mzn_file_name = None
        self.output_file_name = [param["output1"], param["output2"]]
        self.mzn_file_name = "attackprf.mzn"

    #############################################################################################################################################
    #############################################################################################################################################
    #  ____                           _        __                        ____   _       _    _                       _       _
    # / ___|   ___   __ _  _ __  ___ | |__    / _|  ___   _ __    __ _  |  _ \ (_) ___ | |_ (_) _ __    __ _  _   _ (_) ___ | |__    ___  _ __
    # \___ \  / _ \ / _` || '__|/ __|| '_ \  | |_  / _ \ | '__|  / _` | | | | || |/ __|| __|| || '_ \  / _` || | | || |/ __|| '_ \  / _ \| '__|
    #  ___) ||  __/| (_| || |  | (__ | | | | |  _|| (_) || |    | (_| | | |_| || |\__ \| |_ | || | | || (_| || |_| || |\__ \| | | ||  __/| |
    # |____/  \___| \__,_||_|   \___||_| |_| |_|   \___/ |_|     \__,_| |____/ |_||___/ \__||_||_| |_| \__, | \__,_||_||___/|_| |_| \___||_|
    #                                                                                                  |___/
    # Search for a distinguisher using MiniZinc

    def search(self):
        """
        Search for a distinguisher
        """

        if self.time_limit != -1:
            time_limit = datetime.timedelta(seconds=self.time_limit)
        else:
            time_limit = None

        start_time = time.time()
        #############################################################################################################################################
        print(f"Searching for a distinguisher for {self.RD} rounds of Orthros ...")
        self.cp_model = minizinc.Model()
        self.cp_model.add_file(self.mzn_file_name)
        self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
        self.cp_inst["RU"] = self.RU
        self.cp_inst["RM"] = self.RM
        self.cp_inst["RL"] = self.RL
        self.cp_inst["RMU"] = self.RMU
        self.cp_inst["RML"] = self.RML
        self.cp_inst["WU"] = self.WU
        self.cp_inst["WM"] = self.WM
        self.cp_inst["WL"] = self.WL
        self.cp_inst["RB"] = self.RB
        self.cp_inst["offset"] = self.offset
        self.input_diff_middle = ["",""]
        self.output_mask_middle = ["",""]
        self.involved_keys_in_branch = [[], []]
        self.result = self.cp_inst.solve(timeout=time_limit,
                                         processes=self.num_of_threads,
                                         verbose=False,
                                        #  debug_output=Path("./debug_output.txt", intermediate_solutions=True),
                                         random_seed=randint(0, 100))
                                        #  optimisation_level=2)
        #############################################################################################################################################
        elapsed_time = time.time() - start_time
        self.attack_summary = "Time used to find a distinguisher: {:0.02f} seconds".format(elapsed_time) + "\n"
        self.attack_summary += f"Solver status: {self.result.status}" + "\n"
        if minizinc.Status.has_solution(self.result.status) or self.result.status == minizinc.Status.ERROR:
            self.attack_summary_0, self.upper_trail_0, self.lower_trail_0 = self.parse_solution(bn=0)
            self.attack_summary_1, self.upper_trail_1, self.lower_trail_1 = self.parse_solution(bn=1)
            if self.compute_differential_effect_using_gurobi:
                self.diff_effect_0 = self.compute_differential_effect(0)
                self.diff_effect_1 = self.compute_differential_effect(1)
                self.attack_summary_0 += "Differential effect over EU: {:.2f}".format(self.diff_effect_0) + "\n"
                self.attack_summary_1 += "Differential effect over EU: {:.2f}".format(self.diff_effect_1) + "\n"
            self.attack_summary += self.attack_summary_0 + "\n" + "#"*100 + "\n" + self.attack_summary_1

            #############################################################################################################################################
            print("Computing the correlation for the middle part using the experimental approach ...")
            command = [
                "./difflin",
                str(self.RM),
                str(self.offset + self.RB + self.RU),
                self.input_diff_middle[0],
                self.input_diff_middle[1],
                self.output_mask_middle[0],
                str(self.degree),
                "100",
                "0"
            ]
            dictionary_of_correlations_branch_0 = dict()
            dictionary_of_correlations_branch_1 = dict()
            
            # Debug: output mask values (disabled)
            # print(f"Debug - output_mask_middle[0]: {self.output_mask_middle[0]}")
            # print(f"Debug - output_mask_middle[1]: {self.output_mask_middle[1]}")
            
            try:
                for i in range(1, 16):
                    command[5] = "".join([hex(i)[2:] if c != "0" else "0" for c in self.output_mask_middle[0]])
                    # Debug print disabled: print(f"Debug - Branch 0, i={i}, command[5]={command[5]}")
                    result = subprocess.run(command, capture_output=True, text=True, check=True)
                    output_lines = result.stdout.strip().split("\n")
                    correlation_log2 = None
                    sign = None
                    for line in output_lines:
                        if "Log2(|Correlation|)" in line:
                            correlation_log2 = float(line.split(":")[1].strip())
                        elif "Sign" in line:
                            sign = int(line.split(":")[1].strip())
                    if correlation_log2 is not None and sign is not None:
                        total_correlation_0 = -1*self.result["PU0"] + correlation_log2 + -1*self.result["CL0"]
                        r_sign = "+" if sign == 1 else "-"
                        r_value = f"{r_sign}2^({correlation_log2:.2f})"
                        total_value = f"{r_sign}2^({total_correlation_0:.2f})"
                        dictionary_of_correlations_branch_0[command[5]] = (r_value, total_value)
            except subprocess.CalledProcessError as e:
                print("Error occurred while running difflin:")
                print(e.stderr)
            command[8] = "1"
            try:
                for i in range(1, 16):
                    command[5] = "".join([hex(i)[2:] if c != "0" else "0" for c in self.output_mask_middle[1]])
                    # Debug print disabled: print(f"Debug - Branch 1, i={i}, command[5]={command[5]}")
                    result = subprocess.run(command, capture_output=True, text=True, check=True)
                    output_lines = result.stdout.strip().split("\n")
                    correlation_log2 = None
                    sign = None
                    for line in output_lines:
                        if "Log2(|Correlation|)" in line:
                            correlation_log2 = float(line.split(":")[1].strip())
                        elif "Sign" in line:
                            sign = int(line.split(":")[1].strip())
                    if correlation_log2 is not None and sign is not None:
                        total_correlation_1 = -1*self.result["PU1"] + correlation_log2 + -1*self.result["CL1"]
                        r_sign = "+" if sign == 1 else "-"
                        r_value = f"{r_sign}2^({correlation_log2:.2f})"
                        total_value = f"{r_sign}2^({total_correlation_1:.2f})"
                        dictionary_of_correlations_branch_1[command[5]] = (r_value, total_value)
            except subprocess.CalledProcessError as e:
                print("Error occurred while running difflin:")
                print(e.stderr)
            
            # Add correlation summary showing the differential-linear formula breakdown
            self.attack_summary += "\n" + "="*100 + "\n"
            self.attack_summary += "Correlation Summary (Differential-Linear Formula: p × r × q²):\n"
            self.attack_summary += "="*100 + "\n"
            
            # Branch 0 summary
            self.attack_summary += f"\nBranch 0:\n"
            self.attack_summary += f"  Upper differential (p):  2^(-{self.result['PU0']:.2f})\n"
            self.attack_summary += f"  Lower linear (q):        2^(-{self.result['CL0']:.2f})\n"
            self.attack_summary += f"  q² = 2^(-{2*self.result['CL0']:.2f})\n"
            
            # Branch 1 summary
            self.attack_summary += f"\nBranch 1:\n"
            self.attack_summary += f"  Upper differential (p):  2^(-{self.result['PU1']:.2f})\n"
            self.attack_summary += f"  Lower linear (q):        2^(-{self.result['CL1']:.2f})\n"
            self.attack_summary += f"  q² = 2^(-{2*self.result['CL1']:.2f})\n"
            self.attack_summary += "\n" + "="*100 + "\n\n"
            
            self.attack_summary += f"Dictionary of correlations ({self.RD} rounds, Branch 0):\n"
            self.attack_summary += f"Format: output_mask: p × r × q² = corr_tot\n\n"
            for key, (r_value, total_value) in dictionary_of_correlations_branch_0.items():
                self.attack_summary += f"{key}: 2^(-{self.result['PU0']:.2f}) × {r_value} × 2^(-{2*self.result['CL0']:.2f}) = {total_value}\n"

            self.attack_summary += f"\nDictionary of correlations ({self.RD} rounds, Branch 1):\n"
            self.attack_summary += f"Format: output_mask: p × r × q² = corr_tot\n\n"
            for key, (r_value, total_value) in dictionary_of_correlations_branch_1.items():
                self.attack_summary += f"{key}: 2^(-{self.result['PU1']:.2f}) × {r_value} × 2^(-{2*self.result['CL1']:.2f}) = {total_value}\n"
            self.involved_keys_in_branch[0] = list(map(str, [i for i in range(128) if self.result[f"k{0}"][i] == 1]))
            self.involved_keys_in_branch[1] = list(map(str, [i for i in range(128) if self.result[f"k{1}"][i] == 1]))
            self.attack_summary += "\nInvolved key bits in branch 0: " + "{" + ", ".join(self.involved_keys_in_branch[0]) + "}" + "\n"
            self.attack_summary += "\nInvolved key bits in branch 1: " + "{" + ", ".join(self.involved_keys_in_branch[1]) + "}" + "\n"
            intersection = set(self.involved_keys_in_branch[0]) & set(self.involved_keys_in_branch[1])
            self.attack_summary += "\nIntersection of key bits in branch 0 and branch 1: " + "{" + ", ".join(list(intersection)) + "}" + "\n"
            print(self.attack_summary)
            draw_0 = DrawDL(self, output_file_name=self.output_file_name[0], bn=0)
            draw_1 = DrawDL(self, output_file_name=self.output_file_name[1], bn=1)             
            draw_0.generate_distinguisher_shape()
            draw_1.generate_distinguisher_shape()

        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable")
        elif self.result.status == minizinc.Status.UNKNOWN:
            print("Unknown error!")
        else:
            print("Solving process was interrupted")

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

    #############################################################################################################################################
    #############################################################################################################################################
    #  ____                           _    _             ____          _         _    _
    # |  _ \  __ _  _ __  ___   ___  | |_ | |__    ___  / ___|   ___  | | _   _ | |_ (_)  ___   _ __
    # | |_) |/ _` || '__|/ __| / _ \ | __|| '_ \  / _ \ \___ \  / _ \ | || | | || __|| | / _ \ | '_ \
    # |  __/| (_| || |   \__ \|  __/ | |_ | | | ||  __/  ___) || (_) || || |_| || |_ | || (_) || | | |
    # |_|    \__,_||_|   |___/ \___|  \__||_| |_| \___| |____/  \___/ |_| \__,_| \__||_| \___/ |_| |_|
    # Parse the solution and print the distinguisher's specifications

    def parse_solution(self, bn=0):
        """
        Parse the solution and print the distinguisher's specifications
        """

        upper_trail = {"x": [0 for _ in range(self.RU + self.RM + 1)],
                       "y": [0 for _ in range(self.RU + self.RM)],
                       "z": [0 for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):
            upper_trail["x"][r] = self.result[f"xu{bn}"][r]
            upper_trail["y"][r] = self.result[f"yu{bn}"][r]
            upper_trail["z"][r] = self.result[f"zu{bn}"][r]
        for r in range(self.RU, self.RU + self.RM + 1):
            upper_trail["x"][r] = self.result[f"xmu{bn}"][r - self.RU]
            if r < self.RU + self.RM:
                upper_trail["y"][r] = self.result[f"ymu{bn}"][r - self.RU]
                upper_trail["z"][r] = self.result[f"zmu{bn}"][r - self.RU]
        lower_trail = {"x": [0 for _ in range(self.RM + self.RL + 1)],
                       "y": [0 for _ in range(self.RM + self.RL)],
                       "z": [0 for _ in range(self.RM + self.RL)]}
        for r in range(self.RM + 1):
            lower_trail["x"][r] = self.result[f"xml{bn}"][r]
            if r < self.RM:
                lower_trail["y"][r] = self.result[f"yml{bn}"][r]
                lower_trail["z"][r] = self.result[f"zml{bn}"][r]
        for r in range(self.RM, self.RM + self.RL + 1):
            lower_trail["x"][r] = self.result[f"xl{bn}"][r - self.RM]
            if r < self.RM + self.RL:
                lower_trail["y"][r] = self.result[f"yl{bn}"][r - self.RM]
                lower_trail["z"][r] = self.result[f"zl{bn}"][r - self.RM]
        input_diff = ""
        input_diff += f"char inputdiff[] = \"" + hex(int("".join(list(map(str, self.result[f"xu{bn}"][0]))), 2))[2:].zfill(32) + "\";\n"
        temp = self.result[f"xmu{bn}"][0]
        temp = [0 if x == -1 else x for x in temp]
        self.input_diff_middle[bn] = hex(int("".join(list(map(str, temp))), 2))[2:].zfill(32)
        temp = self.result[f"xml{bn}"][self.RM]
        temp = [0 if x == -1 else x for x in temp]
        self.output_mask_middle[bn] = hex(int("".join(list(map(str, temp))), 2))[2:].zfill(32)
        output_mask = ""
        output_mask += f"char outputmask[] = \"" + hex(int("".join(list(map(str, self.result[f"xl{bn}"][self.RL]))), 2))[2:].zfill(32) + "\";\n"

        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: Offset: {self.offset}, RB: {self.RB}, RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}, WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{self.input_diff_middle[bn]}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask. middle: \n{self.output_mask_middle[bn]}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}"
        attack_summary += "#"*50 + "\n"
        attack_summary += "PU{}:  {}\n".format(bn, self.result[f"PU{bn}"])
        attack_summary += "CM{}:   {}\n".format(bn, self.result[f"CM{bn}"])
        attack_summary += "CL{}:   {}\n".format(bn, self.result[f"CL{bn}"])        
        attack_summary += "Number of effective S-boxes in the middle:       {}\n".format(self.result[f"NASM{bn}"])
        attack_summary += "Number of effective bit-positions in the middle: {}\n".format(self.result[f"CM{bn}"])
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
        attack_summary += "Involved key bits in branch {}: ".format(bn) + ", ".join(list(map(str, [i for i in range(128) if self.result[f"k{bn}"][i] == 1]))) + "\n"
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
    params = {"RB": 0,
            "RU": 0,
            "RM": 2,
            "RL": 0,
            "RMU": 0,
            "RML": 0,
            "WU": 1,
            "WM": 1,
            "WL": 1,
            "degree": 15,
            "np" : 8,
            "tl"  : -1,
            "solver"  : "ortools",
            "output1"  : "output1.tex",
            "output2"  : "output2.tex",
            "offset" : 0}

    # Override parameters if they are set on command line
    if args.RB is not None:
        params["RB"] = args.RB
    if args.RU is not None:
        params["RU"] = args.RU
    if args.RM is not None:
        params["RM"] = args.RM
    if args.RL is not None:
        params["RL"] = args.RL
    if args.RMU is not None:
        params["RMU"] = args.RMU
    if args.RML is not None:
        params["RML"] = args.RML
    if args.WU is not None:
        params["WU"] = args.WU
    if args.WM is not None:
        params["WM"] = args.WM
    if args.WL is not None:
        params["WL"] = args.WL
    if args.degree is not None:
        params["degree"] = args.degree
    if args.np is not None:
        params["np"] = args.np
    if args.timelimit is not None:
        params["timelimit"] = args.timelimit
    if args.solver is not None:
        params["solver"] = args.solver
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

    parser = ArgumentParser(description="This tool finds a nearly optimum boomerang"
                                        "distinguisher for SKINNY family of block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-offset", type=int, default=2, help="Number of rounds skipped from the top part of the design")
    parser.add_argument("-RB", type=int, default=1, help="Offset for the starting round of the distinguisher")
    parser.add_argument("-RU", type=int, default=1, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=4, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=0, help="Number of rounds for EL")
    
    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=0, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-WU", type=int, default=5, help="Scale of the probability transition over EU")
    parser.add_argument("-WM", type=int, default=1, help="Scale of the correlation of DL distinguishers over EM ")
    parser.add_argument("-WL", type=int, default=5, help="Scale of the squared correlation of linear approximation over EL")
    parser.add_argument("-d", "--degree", type=int, default=26, help="Number of queries to experimentally compute the correlation: 2^(deg)")
    parser.add_argument("-de", "--differential_effect", type=int, default=0, help="Compute differential effect over EU using Gurobi", choices=[0, 1])

    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=-1, help="Time limit in seconds")
    # Fetch available solvers from MiniZinc
    available_solvers = [solver_name for solver_name in minizinc.default_driver.available_solvers().keys()]
    parser.add_argument("-sl", "--solver", default="cp-sat", type=str,
                        choices=available_solvers,
                        help="Choose a CP solver")
    parser.add_argument("-o1", "--output1", default="output1.tex", type=str, help="Output file name 1")
    parser.add_argument("-o2", "--output2", default="output2.tex", type=str, help="Output file name 2")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DiffLin(params)
    dld.search()

if __name__ == "__main__":
    main()
