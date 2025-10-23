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
import subprocess
from draw import DrawDL
# Check if gurobi_cl command is recognized
try:
    subprocess.run(['gurobi_cl', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    gurobi_available = True
except FileNotFoundError:
    gurobi_available = False

if gurobi_available:
    print("Gurobi is available")
    from diff import Diff
class DiffLin:
    DL_counter = 0

    def __init__(self, param) -> None:
        DiffLin.DL_counter += 1
        self.id = DiffLin.DL_counter
        self.name = "DiffLin" + str(self.id)
        self.type = "DiffLin"
        self.RU = param["RU"]
        self.RM = param["RM"]
        self.RL = param["RL"]
        self.RMU = param["RMU"]
        self.RML = param["RML"]
        self.WU = param["WU"]
        self.WM = param["WM"]
        self.WL = param["WL"]
        self.offset = param["offset"]
        self.branchtype = param["branchtype"]
        self.RD = self.RU + self.RM + self.RL
        self.cp_solver_name = param["solver"]
        self.compute_differential_effect_using_gurobi = param["differential_effect"]
        self.cp_solver = minizinc.Solver.lookup(self.cp_solver_name)
        self.time_limit = param["timelimit"]
        self.num_of_threads = param["np"]
        self.mzn_file_name = None
        self.output_file_name = param["output"]
        self.mzn_file_name = "attackbranch.mzn"
    
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
        print(f"Searching a distinguisher for {self.RD} rounds of Orthros ...")
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
        self.cp_inst["offset"] = self.offset
        self.cp_inst["branchtype"] = self.branchtype
        self.result = self.cp_inst.solve(timeout=time_limit, 
                                         processes=self.num_of_threads, 
                                         verbose=False, 
                                        #  debug_output=Path("./debug_output.txt",intermediate_solutions=True),
                                         random_seed=randint(0, 100))
                                        #  optimisation_level=2)
        #############################################################################################################################################
        elapsed_time = time.time() - start_time
        self.attack_summary = "Time used to find a distinguisher: {:0.02f} seconds".format(elapsed_time) + "\n"
        self.attack_summary += f"Solver status: {self.result.status}" + "\n"
        if minizinc.Status.has_solution(self.result.status) or self.result.status == minizinc.Status.ERROR:
            attack_summary, self.upper_trail, self.lower_trail = self.parse_solution()
            self.attack_summary += attack_summary
            if self.compute_differential_effect_using_gurobi:
                diff_effect = self.compute_differential_effect()
                self.attack_summary += "Differential effect over EU: 2^({:.2f})".format(diff_effect) + "\n"
            print(self.attack_summary)
            draw = DrawDL(self, output_file_name=self.output_file_name)
            draw.generate_distinguisher_shape()
            print("-Log2(P)              ~= \t{:02d}".format(self.result["PU"]))
            print("-Log2(r)              ~= \t{:02d}".format(self.result["CM"]))
            print("-Log2(Q^2)            ~= \t{:02d}".format(self.result["CL"]))
            if self.compute_differential_effect_using_gurobi:
                print("-Log2(DiffEffect)     ~= \t{:.2f}".format(-1*diff_effect))  
        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable") 
        elif self.result.status == minizinc.Status.UNKNOWN:
            print("Unknown error!")
        else:
            print("Solving process was interrupted")
    
    def compute_differential_effect(self):
        """
        Compute the differential effect over EU
        """

        # Load default values
        params = {"nrounds" : self.RU,
                  "branchtype" : self.branchtype,
                  "offset" : self.offset,
                  "mode" : 2,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 2,
                  "fixedVariables" : {}}
        input_diff = hex(int("".join(list(map(str, self.result["xu"][0]))), 2))[2:].zfill(32)
        params["fixedVariables"][f"x_0"] = input_diff
        outputdiff = hex(int("".join(list(map(str, self.result["xu"][self.RU]))), 2))[2:].zfill(32)
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

    def parse_solution(self):
        """
        Parse the solution and print the distinguisher's specifications
        """
        
        upper_trail = {"x": [0 for _ in range(self.RU + self.RM + 1)],
                       "y": [0 for _ in range(self.RU + self.RM)],
                       "z": [0 for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):
            upper_trail["x"][r] = self.result["xu"][r]
            upper_trail["y"][r] = self.result["yu"][r]
            upper_trail["z"][r] = self.result["zu"][r]
        for r in range(self.RU, self.RU + self.RM + 1):
            upper_trail["x"][r] = self.result["xmu"][r - self.RU]
            if r < self.RU + self.RM:
                upper_trail["y"][r] = self.result["ymu"][r - self.RU]
                upper_trail["z"][r] = self.result["zmu"][r - self.RU]
        lower_trail = {"x": [0 for _ in range(self.RM + self.RL + 1)],
                       "y": [0 for _ in range(self.RM + self.RL)],
                       "z": [0 for _ in range(self.RM + self.RL)]}
        for r in range(self.RM + 1):
            lower_trail["x"][r] = self.result["xml"][r]
            if r < self.RM:
                lower_trail["y"][r] = self.result["yml"][r]
                lower_trail["z"][r] = self.result["zml"][r]
        for r in range(self.RM, self.RM + self.RL + 1):
            lower_trail["x"][r] = self.result["xl"][r - self.RM]
            if r < self.RM + self.RL:
                lower_trail["y"][r] = self.result["yl"][r - self.RM]
                lower_trail["z"][r] = self.result["zl"][r - self.RM]
        input_diff = ""
        input_diff += f"char inputdiff[] = \"" + hex(int("".join(list(map(str, self.result["xu"][0]))), 2))[2:].zfill(32) + "\";\n"
        input_diff_middle = ""
        temp = self.result["xmu"][0]
        temp = [0 if x == -1 else x for x in temp]
        input_diff_middle += f"char inputdiff[] = \"" + hex(int("".join(list(map(str, temp))), 2))[2:].zfill(32) + "\";\n"
        output_mask_middle = ""
        temp = self.result["xml"][self.RM]
        temp = [0 if x == -1 else x for x in temp]
        output_mask_middle += f"char outputmask[] = \"" + hex(int("".join(list(map(str, temp))), 2))[2:].zfill(32) + "\";\n"
        output_mask = ""
        output_mask += f"char outputmask[] = \"" + hex(int("".join(list(map(str, self.result["xl"][self.RL]))), 2))[2:].zfill(32) + "\";\n"
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}, WU: {self.WU}, WM: {self.WM}, WL: {self.WL}, offset: {self.offset}, branchtype: {self.branchtype}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{input_diff_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask. middle: \n{output_mask_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"PU:  {self.result['PU']}\n"
        attack_summary += f"CM:  {self.result['CM']}\n"
        attack_summary += f"Q^2: {self.result['CL']}\n"
        attack_summary += f"Number of effective S-boxes in the middle:       {self.result['NASM']}\n"
        attack_summary += f"Number of effective bit-positions in the middle: {self.result['CM']}\n"
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
    params = {"RU": 0,
            "RM": 4,
            "RL": 0,
            "RMU": 0,
            "RML": 0,
            "WU": 1,
            "WM": 1,
            "WL": 1,
            "offset": 4,
            "branchtype": 1,
            "np" : 8,
            "tl"  : -1,
            "solver"  : "ortools",
            "output"  : "output.tex",
            "differential_effect" : 0}

    # Override parameters if they are set on command line
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
    if args.offset is not None:
        params["offset"] = args.offset
    if args.branchtype is not None:
        params["branchtype"] = args.branchtype
    if args.np is not None:
        params["np"] = args.np
    if args.timelimit is not None:
        params["timelimit"] = args.timelimit
    if args.solver is not None:
        params["solver"] = args.solver
    if args.output is not None:
        params["output"] = args.output
    if args.differential_effect is not None:
        params["differential_effect"] = args.differential_effect

    return params

def main():
    '''
    Parse the arguments and start the request functionality with the provided
    parameters.
    '''
    
    parser = ArgumentParser(description="This tool finds a nearly optimum boomerang"
                                        "distinguisher for SKINNY family of block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-RU", type=int, default=1, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=6, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=0, help="Number of rounds for EL")
    parser.add_argument("-offset", type=int, default=4, help="Offset for the starting round of the distinguisher")
    parser.add_argument("-branchtype", type=int, default=0, help="Branch type of Orthros (0: Branch1, 1: Branch2)", choices=[0, 1])


    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=0, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-WU", type=int, default=1, help="Scale of the probability transition over EU")
    parser.add_argument("-WM", type=int, default=1, help="Scale of the correlation of DL distinguishers over EM ")
    parser.add_argument("-WL", type=int, default=1, help="Scale of the squared correaltion of linear approximation over EL")

    parser.add_argument("-de", "--differential_effect", type=int, default=0, help="Compute differential effect over EU using Gurobi", choices=[0, 1])

    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=36000, help="Time limit in seconds")
    # Fetch available solvers from MiniZinc
    available_solvers = [solver_name for solver_name in minizinc.default_driver.available_solvers().keys()]
    parser.add_argument("-sl", "--solver", default="cp-sat", type=str,
                        choices=available_solvers,
                        help="Choose a CP solver")      
    parser.add_argument("-o", "--output", default="output.tex", type=str, help="Output file name")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DiffLin(params)
    dld.search()

if __name__ == "__main__":
    main()