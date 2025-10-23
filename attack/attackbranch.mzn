/*
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
*/

% Strategy:
% Divide E into 4 parts: EU, EMU0, EMU1
%                            EML0, EML1, EL
% Use DDT (LAT) to encode the propagation of differential (resp. linear) trails through EU (resp. EL).
% Use *-DDT (*-LAT) to encode the propagation of differential (resp. linear) trails through EMU0 (resp. EML1).
% Enode the propagation of differential (resp. linear) trails through EMU1 (resp. EML0) with probability one.
% Minimize the probability of transtions through EU, EL as well as the number of common active S-boxes trhough
% the middle part, i.e., EMU0, EMU1, EML0, EML1.
% ###############################################################################################################
% We use S-box Analyzer to encode the differential and linear behavior of the S-boxes. 
% sage: from sboxanalyzer import *
% sage: sa = SboxAnalyzer([0x1,0x0,0x2,0x4,0x3,0x8,0x6,0xd,0x9,0xa,0xb,0xe,0xf,0xc,0x7,0x5])


% sage: cnf0, milp0 = sa.minimized_diff_constraints()
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.02 seconds
% Number of constraints: 64
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% Weight: 3.0000 p0 + 2.0000 p1

% sage: cnf0star, milp0star = sa.minimized_diff_constraints(subtable='star')
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.00 seconds
% Number of constraints: 54
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb

% sage: cnf1, milp1 = sainv.minimized_linear_constraints()
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.01 seconds
% Number of constraints: 74
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% Weight: 4.0000 p0 + 2.0000 p1

% sage: cnf1star, milp1star = sainv.minimized_linear_constraints(subtable='star')
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.00 seconds
% Number of constraints: 44
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb

% sage: detdiff = sa.encode_deterministic_differential_behavior()
% sage: cpdetdiff = sa.generate_cp_constraints(detdiff); print(cpdetdiff)
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0)
% else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1)
% endif

% sage: sainv = SboxAnalyzer(sa.inverse())
% sage: detlin = sainv.encode_deterministic_linear_behavior(); detlin
% {(0, 0, 0, 0): [0, 0, 0, 0], (1, 1, 1, 0): [-1, -1, 1, -1]}
% sage: cpdetlin = sainv.generate_cp_constraints(detlin); print(cpdetlin)
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0)
% elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1)
% else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1)
% endif
% ###############################################################################################################

int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: WU;
int: WM;
int: WL;
int: offset;
int: branchtype;

constraint assert(RU >= 0, "Invalid value for RU: " ++
        "RU must be greater than or equal to 0");
constraint assert(RMU >= 0, "Invalid value for RMU: " ++
        "RMU must be greater than or equal to 0");
constraint assert(RM >= 0, "Invalid value for RM: " ++
        "RM must be greater than or equal to 0");
constraint assert(RML >= 0, "Invalid value for RML: " ++
        "RML must be greater than or equal to 0");
constraint assert(RL >= 0, "Invalid value for RL: " ++
        "RL must be greater than or equal to 0");
constraint assert(RMU <= RM, "Invalid value for RMU: " ++
        "RMU must be less than RM");
constraint assert(RML <= RM, "Invalid value for RML: " ++ 
        "RML must be less than RM");
constraint assert(offset >= 0, "Invalid value for offset: " ++
        "offset must be greater than or equal to 0");
constraint assert(branchtype in 0..1, "Invalid value for branchtype: " ++
        "branchtype must be 0 or 1");
% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____          __  _                _                   _                    _____       _      _            
% |  _ \   ___  / _|(_) _ __    ___  | |     ___    ___  | | __ _   _  _ __   |_   _|__ _ | |__  | |  ___  ___ 
% | | | | / _ \| |_ | || '_ \  / _ \ | |    / _ \  / _ \ | |/ /| | | || '_ \    | | / _` || '_ \ | | / _ \/ __|
% | |_| ||  __/|  _|| || | | ||  __/ | |___| (_) || (_) ||   < | |_| || |_) |   | || (_| || |_) || ||  __/\__ \
% |____/  \___||_|  |_||_| |_| \___| |_____|\___/  \___/ |_|\_\ \__,_|| .__/    |_| \__,_||_.__/ |_| \___||___/
%                                                                     |_|                                      
% define lookup tables

array[0..1, 0..127] of int: bit_permutation = array2d(0..1, 0..127, 
    [6,46,62,126,70,52,28,14,36,125,72,83,106,95,4,35, 25,41,10,76,87,74,120,42,88,21,11,67,64,38,112,50, 85,109,24,65,99,0,49,37,8,66,114,47,127,100,56,40, 13,117,78,86,92,58,124,101,55,89,97,9,18,116,59,15, 20,45,75,2,77,27,1,60,115,107,26,69,119,3,84,51, 123,110,31,82,113,53,81,102,63,118,93,12,30,94,108,32, 5,111,29,43,91,19,79,33,73,44,98,48,22,61,68,105, 34,71,54,104,17,57,80,103,96,121,23,39,122,90,7,16,
     20,122,74,62,119,35,15,66,9,85,32,117,21,83,127,106, 11,98,115,59,71,90,56,26,2,44,103,121,114,107,68,16, 84,1,102,33,80,52,76,36,27,94,37,55,82,12,112,64, 105,14,91,17,108,124,6,93,29,86,123,79,72,53,19,99, 50,18,81,73,67,88,4,61,111,49,24,45,57,78,100,22, 110,47,116,54,60,70,97,39,3,41,48,96,23,42,113,87, 126,13,31,40,51,25,65,125,8,101,118,28,38,89,5,104, 109,120,69,43,7,77,58,34,10,63,30,95,75,46,0,92]);

array[0..1, 0..31] of int: nibble_permutation = array2d(0..1, 0..31, 
    [10, 27, 5, 1, 30, 23, 16, 13, 21, 31, 6, 14, 0, 25, 11, 18, 15, 28, 19, 24, 7, 8, 22, 3, 4, 29, 9, 2, 26, 20, 12, 17,
     26, 13, 7, 11, 29, 0, 17, 21, 23, 5, 18, 25, 12, 10, 28, 2, 14, 19, 24, 22, 1, 8, 4, 31, 15, 6, 27, 9, 16, 30, 20, 3]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..127] of var 0..1: xu;
array[0..(RU - 1), 0..127] of var 0..1: yu;
array[0..(RU - 1), 0..127] of var 0..1: zu;
array[0..(RU - 1), 0..31, 0..1] of var 0..1: pu;

% S-box layer
constraint forall(round in 0..(RU - 1), i in 0..31)
(
    sbox_diff(xu[round, 4*i], xu[round, 4*i + 1], xu[round, 4*i + 2], xu[round, 4*i + 3], 
              yu[round, 4*i], yu[round, 4*i + 1], yu[round, 4*i + 2], yu[round, 4*i + 3],
              pu[round, i, 0], pu[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RU - 1), i in 0..31, j in 0..3)
(
    if (offset + round) < 4 then 
    (
        zu[round, bit_permutation[branchtype, 4*i + j]] = yu[round, 4*i + j]
    ) else 
    (
        zu[round, 4*nibble_permutation[branchtype, i] + j] = yu[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RU - 1), i in 0..7, j in 0..3)
(
    xor3(zu[round, 16*i + 4*1 + j], zu[round, 16*i + 4*2 + j], zu[round, 16*i + 4*3 + j], xu[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zu[round, 16*i + 4*0 + j], zu[round, 16*i + 4*2 + j], zu[round, 16*i + 4*3 + j], xu[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zu[round, 16*i + 4*0 + j], zu[round, 16*i + 4*1 + j], zu[round, 16*i + 4*3 + j], xu[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zu[round, 16*i + 4*0 + j], zu[round, 16*i + 4*1 + j], zu[round, 16*i + 4*2 + j], xu[round  +1, 16*i + 4*3 + j])
);

% Exclude the all-zero input/output difference 
constraint sum(i in 0..127) (xu[0, i]) != 0;

% Compute the probability of the differential transition
var 0..128: PU;
constraint PU = sum(r in 0..(RU - 1), i in 0..31)(3*pu[r, i, 0] + 2*pu[r, i, 1]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU                                                                                                                                                             

array[0..RM, 0..127] of var -1..1: xmu;
array[0..(RM - 1), 0..127] of var -1..1: ymu;
array[0..(RM - 1), 0..127] of var -1..1: zmu;
array[0..(RM - 1), 0..31] of var 0..1: dmu;

% link the input of EM to the output of EU
constraint forall(i in 0..127)(xmu[0, i] = xu[RU, i]);

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu[round, 4*i], xmu[round, 4*i + 1], xmu[round, 4*i + 2], xmu[round, 4*i + 3], 
                       ymu[round, 4*i], ymu[round, 4*i + 1], ymu[round, 4*i + 2], ymu[round, 4*i + 3],
                       dmu[round, i])
    ) else 
    (
        sbox_diff_m(xmu[round, 4*i], xmu[round, 4*i + 1], xmu[round, 4*i + 2], xmu[round, 4*i + 3], 
                    ymu[round, 4*i], ymu[round, 4*i + 1], ymu[round, 4*i + 2], ymu[round, 4*i + 3],
                    dmu[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RU + round) < 4 then 
    (
        zmu[round, bit_permutation[branchtype, 4*i + j]] = ymu[round, 4*i + j]
    ) else 
    (
        zmu[round, 4*nibble_permutation[branchtype, i] + j] = ymu[round, 4*i + j]
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < RMU then 
    (
        xor3(zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*2 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*0 + j])
        /\
        xor3(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*2 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*1 + j])
        /\
        xor3(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*2 + j])
        /\
        xor3(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*2 + j], xmu[round + 1, 16*i + 4*3 + j])
    ) else 
    (
        xor3_pr1(zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*2 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*0 + j])
        /\
        xor3_pr1(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*2 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*1 + j])
        /\
        xor3_pr1(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*3 + j], xmu[round + 1, 16*i + 4*2 + j])
        /\
        xor3_pr1(zmu[round, 16*i + 4*0 + j], zmu[round, 16*i + 4*1 + j], zmu[round, 16*i + 4*2 + j], xmu[round + 1, 16*i + 4*3 + j])
    ) endif
    
);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _     
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || |    
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || |    
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |___ 
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_||_____|
% Constraints for EML                                                                                                     

array[0..RM, 0..127] of var -1..1: xml;
array[0..(RM - 1), 0..127] of var -1..1: yml;
array[0..(RM - 1), 0..127] of var -1..1: zml;
array[0..(RM - 1), 0..31] of var 0..1: dml;

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml[round, 4*i], yml[round, 4*i + 1], yml[round, 4*i + 2], yml[round, 4*i + 3], 
                   xml[round, 4*i], xml[round, 4*i + 1], xml[round, 4*i + 2], xml[round, 4*i + 3],
                   dml[round, i])
    ) else
    (
        sbox_lin_star(yml[round, 4*i], yml[round, 4*i + 1], yml[round, 4*i + 2], yml[round, 4*i + 3], 
                      xml[round, 4*i], xml[round, 4*i + 1], xml[round, 4*i + 2], xml[round, 4*i + 3],
                      dml[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RU + round) < 4 then 
    (
        zml[round, bit_permutation[branchtype, 4*i + j]] = yml[round, 4*i + j]
    ) else 
    (
        zml[round, 4*nibble_permutation[branchtype, i] + j] = yml[round, 4*i + j]
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < (RM - RML) then 
    (
        xor3_pr1(xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*2 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*0 + j])
        /\
        xor3_pr1(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*2 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*1 + j])
        /\
        xor3_pr1(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*2 + j])
        /\
        xor3_pr1(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*2 + j], zml[round, 16*i + 4*3 + j])
    ) else
    (
        xor3(xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*2 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*0 + j])
        /\
        xor3(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*2 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*1 + j])
        /\
        xor3(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*3 + j], zml[round, 16*i + 4*2 + j])
        /\
        xor3(xml[round + 1, 16*i + 4*0 + j], xml[round + 1, 16*i + 4*1 + j], xml[round + 1, 16*i + 4*2 + j], zml[round, 16*i + 4*3 + j])
    ) endif
);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _     
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| |    
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |    
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |___ 
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_____|
% constraints for EL

array[0..RL, 0..127] of var 0..1: xl;
array[0..(RL - 1), 0..127] of var 0..1: yl;
array[0..(RL - 1), 0..127] of var 0..1: zl;
array[0..(RL - 1), 0..31, 0..1] of var 0..1: cl;

% link the input of EL to the output of EML
constraint forall(i in 0..127)(xl[0, i] = xml[RM, i]);

% S-box layer
constraint forall(round in 0..(RL - 1), i in 0..31)
(
    sbox_lin(yl[round, 4*i], yl[round, 4*i + 1], yl[round, 4*i + 2], yl[round, 4*i + 3],
             xl[round, 4*i], xl[round, 4*i + 1], xl[round, 4*i + 2], xl[round, 4*i + 3],
             cl[round, i, 0], cl[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RL - 1), i in 0..31, j in 0..3)
(
    if (offset + RU + RM + round) < 4 then 
    (
        zl[round, bit_permutation[branchtype, 4*i + j]] = yl[round, 4*i + j]
    ) else 
    (
        zl[round, 4*nibble_permutation[branchtype, i] + j] = yl[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RL - 1), i in 0..7, j in 0..3)
(
    xor3(zl[round, 16*i + 4*1 + j], zl[round, 16*i + 4*2 + j], zl[round, 16*i + 4*3 + j], xl[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zl[round, 16*i + 4*0 + j], zl[round, 16*i + 4*2 + j], zl[round, 16*i + 4*3 + j], xl[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zl[round, 16*i + 4*0 + j], zl[round, 16*i + 4*1 + j], zl[round, 16*i + 4*3 + j], xl[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zl[round, 16*i + 4*0 + j], zl[round, 16*i + 4*1 + j], zl[round, 16*i + 4*2 + j], xl[round + 1, 16*i + 4*3 + j])
);

% Exclude the all-zero output mask
constraint sum(i in 0..127) (xl[RL, i]) != 0;

% Compute the squared correlation of the linear trail
var 0..128: CL;
constraint CL = sum(r in 0..(RL - 1), i in 0..31)(4*cl[r, i, 0] + 2*cl[r, i, 1]);

% #############################################################################################################################################
% #############################################################################################################################################

var 0..RM*(128): NASM;
var 0..RM*(128): CM;
constraint CM = 
                1*sum(r in 0..(RM - 1), i in 0..127)(
                bool2int(ymu[r, i] = -1 /\ yml[r, i] = -1)
                +
                bool2int(ymu[r, i] = -1 /\ yml[r, i] = 1)
                +                 
                bool2int(ymu[r, i] = 1 /\ yml[r, i] = -1));

constraint NASM = sum(r in 0..(RM - 1), i in 0..31)(bool2int(dmu[r, i] + dml[r, i] == 2));

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____                           _       ____   _                _                      
% / ___|   ___   __ _  _ __  ___ | |__   / ___| | |_  _ __  __ _ | |_  ___   __ _  _   _ 
% \___ \  / _ \ / _` || '__|/ __|| '_ \  \___ \ | __|| '__|/ _` || __|/ _ \ / _` || | | |
%  ___) ||  __/| (_| || |  | (__ | | | |  ___) || |_ | |  | (_| || |_|  __/| (_| || |_| |
% |____/  \___| \__,_||_|   \___||_| |_| |____/  \__||_|   \__,_| \__|\___| \__, | \__, |
%                                                                           |___/  |___/ 
% Search strategy
ann: search_ann;
% search_ann = seq_search([
%     int_search([yl[RL - 1, 0, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 1, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 2, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..63], indomain_min, complete)]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ___   _      _              _    _               _____                     _    _               
%  / _ \ | |__  (_)  ___   ___ | |_ (_)__   __ ___  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __  
% | | | || '_ \ | | / _ \ / __|| __|| |\ \ / // _ \ | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ 
% | |_| || |_) || ||  __/| (__ | |_ | | \ V /|  __/ |  _| | |_| || | | || (__ | |_ | || (_) || | | |
%  \___/ |_.__/_/ | \___| \___| \__||_|  \_/  \___| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_|
%             |__/                                                                                  
% Objective function

solve minimize (WU*PU + WM*CM + WL*CL);
% solve minimize (WU*PU + WM*NASM + WL*CL);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%     _                 _  _  _                       _____                     _    _                    
%    / \   _   _ __  __(_)| |(_)  __ _  _ __  _   _  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___ 
%   / _ \ | | | |\ \/ /| || || | / _` || '__|| | | | | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|
%  / ___ \| |_| | >  < | || || || (_| || |   | |_| | |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \
% /_/   \_\\__,_|/_/\_\|_||_||_| \__,_||_|    \__, | |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/
%                                             |___/       
% auxiliary functions

predicate xor3(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: b) = 
    b - a2 - a1 - a0 >= -2
    /\        
    - b + a2 - a1 - a0 >= -2
    /\
    - b - a2 + a1 - a0 >= -2
    /\
    b + a2 + a1 - a0 >= 0
    /\
    - b - a2 - a1 + a0 >= -2
    /\
    b + a2 - a1 + a0 >= 0
    /\
    b - a2 + a1 + a0 >= 0
    /\
    - b + a2 + a1 + a0 >= 0;

predicate xor3_pr1(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: y) = 
    if exists([x0 == -1, x1 == -1, x2 == -1]) then y = -1
    else y = (sum([x0, x1, x2]) mod 2) endif;

predicate sbox_diff(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, 
                   var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, 
                   var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a0 - a1 + a2 + p0 >= -1
    /\
    a0 + a3 - b0 + p0 >= 0
    /\
    a1 - a2 + b1 + p0 >= 0
    /\
    - a1 + b0 + b1 + p0 >= 0
    /\
    - a1 - b0 - b2 + p0 >= -2
    /\
    a1 - a3 + b2 + p0 >= 0
    /\
    - a0 + b1 - b3 + p0 >= -1
    /\
    a3 - b2 - b3 + p0 >= -1
    /\
    a2 - a3 + b3 + p0 >= 0
    /\
    b0 + b1 + b3 - p1 >= 0
    /\
    a0 + a2 + a3 + b0 - b3 >= 0
    /\
    a0 + a1 - b0 - b1 - b3 >= -2
    /\
    a2 + a3 - b1 - b2 - b3 >= -2
    /\
    - a0 - a3 + b0 - b2 + p0 >= -2
    /\
    a1 + b0 - b1 - b2 + p0 >= -1
    /\
    - a0 + a3 - b1 + b2 + p0 >= -1
    /\
    - a1 - a2 + b2 - b3 + p0 >= -2
    /\
    - a0 - a3 - b0 + b3 + p0 >= -2
    /\
    a2 - b1 + b2 + b3 + p0 >= 0
    /\
    a0 + a1 + a2 - b0 + p1 >= 0
    /\
    a1 + a2 + a3 - b2 + p1 >= 0
    /\
    a1 - b0 - b1 - b2 + p1 >= -2
    /\
    a0 - a3 + b2 + b3 + p1 >= 0
    /\
    a0 - a1 - a2 - b0 + b1 - b3 >= -3
    /\
    - a0 - a1 + a3 - b0 + b1 - b3 >= -3
    /\
    - a1 - a2 + a3 + b1 - b2 - b3 >= -3
    /\
    - a0 - a2 - a3 + b1 + b2 - b3 >= -3
    /\
    - a1 - a2 - a3 + b0 + b1 + b3 >= -2
    /\
    - a0 + a2 - a3 + b0 + b1 + b3 >= -1
    /\
    - a0 + a1 + a3 + b0 + b1 + b3 >= 0
    /\
    - a0 + a1 + a2 - b1 - b2 + b3 >= -2
    /\
    a0 - a1 - a2 + b1 + b2 + b3 >= -1
    /\
    a0 - a1 - a2 + b0 + b3 + p0 >= -1
    /\
    - a0 - a1 - a2 + a3 + b1 + p1 >= -2
    /\
    - a1 - a3 + b0 - b1 + b2 + p1 >= -2
    /\
    a0 - a3 - b0 + b1 + b2 + p1 >= -1
    /\
    a1 + a2 - b0 - b1 - b3 + p1 >= -2
    /\
    a1 + a3 + b0 - b1 - b3 + p1 >= -1
    /\
    - a0 - a1 - a2 + a3 + b3 + p1 >= -2
    /\
    a0 + a1 - a2 - b1 + b3 + p1 >= -1
    /\
    a1 - a2 - b0 - b1 + b3 + p1 >= -2
    /\
    - a0 - a1 - a2 + b1 + b3 + p1 >= -2
    /\
    a3 + b0 - b1 - b2 + b3 + p1 >= -1
    /\
    - a0 - a1 + b0 + b2 + b3 + p1 >= -1
    /\
    - a0 + a3 + b1 + b2 + b3 + p1 >= 0
    /\
    a0 - a1 + a2 - a3 + b0 - b1 - b2 >= -3
    /\
    a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 + b0 - b2 + b3 >= -4
    /\
    a0 - a1 - a3 + b0 - b1 - b3 + p1 >= -3
    /\
    - a2 - a3 - b0 - b1 - b2 - b3 + p1 >= -5
    /\
    a0 + a3 + b0 - b1 + b2 - b3 + p1 >= -1
    /\
    a0 - a2 - a3 - b0 - b1 + b3 + p1 >= -3
    /\
    - a0 + a1 + a2 + b1 + b2 + b3 + p1 >= 0
    /\
    a0 + a1 + a3 + b0 + b2 - p0 + p1 >= 0
    /\
    a0 + a1 + a2 - a3 + b0 - b2 - b3 + p1 >= -2
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 - p0 + p1 >= -4
    /\
    a0 - a2 + a3 - b0 + b1 - b2 - p0 + p1 >= -3
    /\
    a0 - a1 + a2 + b0 + b1 + b2 - p0 + p1 >= -1
    /\
    - a0 + a1 - a2 + b0 - b2 - b3 - p0 + p1 >= -4
    /\
    a0 - a2 - a3 + b0 + b1 - b2 - b3 - p0 + p1 >= -4
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 + b3 - p0 + p1 >= -5
    /\
    a0 - a1 + a2 + a3 - b0 - b2 + b3 - p0 + p1 >= -3
    /\
    - a0 + a1 - a3 + b0 + b2 >= -1;

predicate sbox_diff_star(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, 
                        var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, 
                        var 0..1: d) =
    a0 + a1 + a3 - b0 - b1 >= -1
    /\
    - a0 + a1 + a3 + b0 - b2 >= -1
    /\
    a1 + a2 + a3 + b0 - b2 >= 0
    /\
    a0 + a1 + a2 - b0 + b2 >= 0
    /\
    a0 + a1 - a3 - b0 + b2 >= -1
    /\
    a0 + a2 + a3 + b0 - b3 >= 0
    /\
    a0 + a1 - b0 - b1 - b3 >= -2
    /\
    a1 + a2 + a3 - b2 - b3 >= -1
    /\
    a0 + a1 + a2 - b0 + b3 >= 0
    /\
    a0 - a3 + b0 + b2 + b3 >= 0
    /\
    - a0 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 - a2 + a3 + b0 + b1 >= -2
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    a0 - a2 + a3 - b0 + b1 - b2 >= -2
    /\
    a0 - a1 - a2 - b0 + b1 - b3 >= -3
    /\
    - a0 - a1 + a3 - b0 + b1 - b3 >= -3
    /\
    - a0 + a1 - a2 + b0 - b2 - b3 >= -3
    /\
    - a1 + a2 + a3 - b1 - b2 - b3 >= -3
    /\
    a0 + a2 + b0 - b1 - b2 - b3 >= -2
    /\
    a1 + a2 - b0 - b1 + b2 - b3 >= -2
    /\
    - a0 + a1 + b0 - b1 + b2 - b3 >= -2
    /\
    - a0 - a2 - a3 + b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 + b1 + b3 >= -3
    /\
    a0 + a2 - a3 - b0 + b1 + b3 >= -1
    /\
    - a1 - a2 - a3 + b0 + b1 + b3 >= -2
    /\
    - a0 + a2 - a3 + b0 + b1 + b3 >= -1
    /\
    a0 + a2 + a3 - b0 - b2 + b3 >= -1
    /\
    - a0 + a1 + a2 - b1 - b2 + b3 >= -2
    /\
    a0 + a1 + a2 - b1 + b2 + b3 >= 0
    /\
    - a1 + a2 - a3 - b1 + b2 + b3 >= -2
    /\
    a0 - a1 - a2 + b1 + b2 + b3 >= -1
    /\
    a0 + a1 - a2 + a3 + b0 + b1 + b2 >= 0
    /\
    a0 - a2 - a3 + b0 + b1 - b2 - b3 >= -3
    /\
    a0 - a1 + a3 + b0 + b1 - b2 - b3 >= -2
    /\
    a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -3
    /\
    a0 - a1 - a2 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a3 + b0 - b1 + b2 - b3 >= -4
    /\
    - a0 + a1 - a2 - a3 - b0 - b1 + b3 >= -4
    /\
    a0 - a1 + a2 - a3 + b0 - b1 + b3 >= -2
    /\
    - a0 - a1 + a2 + a3 + b0 - b1 + b3 >= -2
    /\
    - a0 - a1 - a2 + a3 - b0 - b2 + b3 >= -4
    /\
    - a0 + a2 - a3 - b0 - b1 - b2 + b3 >= -4
    /\
    a0 + a1 - a2 + b0 - b1 - b2 + b3 >= -2
    /\
    - a0 + a1 - a2 + a3 - b0 + b2 + b3 >= -2
    /\
    - a0 - a1 - a2 + a3 - b1 + b2 + b3 >= -3
    /\
    - a0 - a1 + a2 + a3 + b1 + b2 + b3 >= -1
    /\
    - a0 - a1 - a2 - a3 + b0 - b2 + b3 >= -4
    /\
    - a1 - a2 - a3 - b0 - b1 - b2 - b3 >= -6
    /\
    a0 - a1 - a2 - a3 - b0 - b1 - b2 >= -5
    /\
    a0 - a2 + a3 + b0 - b1 - b2 + b3 >= -2
    /\
    a1 + a2 - a3 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 + a2 + b0 + b1 + b2 >= 0
    /\
    a1 + a3 - b1 - b2 - b3 >= -2
    /\
    - a0 + a1 - a3 + b0 + b2 >= -1
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1]) then d = 1
    else d = 0 endif;

predicate sbox_lin(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, 
                   var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, 
                   var 0..1: p0, var 0..1: p1) = 
    - b0 + p0 + p1 >= 0
    /\
    - b1 + p0 + p1 >= 0
    /\
    - b2 + p0 + p1 >= 0
    /\
    - b3 + p0 + p1 >= 0
    /\
    - a0 + a2 + b0 - p1 >= -1
    /\
    a0 + a1 + a2 - b0 - b1 >= -1
    /\
    - a0 + b0 + b1 + b2 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a2 - b0 + b1 - p1 >= -1
    /\
    a0 - b0 - b1 - b2 - p1 >= -3
    /\
    - a1 + a3 - b0 + b2 - p1 >= -2
    /\
    - a0 - a2 + b1 + b2 - p1 >= -2
    /\
    a2 + a3 + b0 - b3 - p1 >= -1
    /\
    - a1 - a3 + b1 - b3 - p1 >= -3
    /\
    a3 - b0 + b2 - b3 - p1 >= -2
    /\
    a0 + a1 + a3 + b3 - p1 >= 0
    /\
    a1 - a2 - b0 + b3 - p1 >= -2
    /\
    a2 - a3 - b1 + b3 - p1 >= -2
    /\
    - a2 + a3 + b1 + b3 - p1 >= -1
    /\
    a1 - a3 + b2 + b3 - p1 >= -1
    /\
    a0 + b1 + b2 + b3 - p1 >= 0
    /\
    - a1 - a2 - b0 - b1 - b2 - b3 >= -5
    /\
    a1 + a2 - b0 - b1 - b2 - b3 >= -3
    /\
    - a0 + a1 - a2 - a3 + b1 + b3 >= -2
    /\
    a0 - a1 + a2 - a3 + b1 + b3 >= -1
    /\
    - a0 - a1 + a3 - b0 - b2 - p0 >= -4
    /\
    a0 + a1 + a3 + b0 + b2 - p0 >= 0
    /\
    - a0 - a2 - b0 - b1 - b3 - p0 >= -5
    /\
    - a1 - a2 + a3 - b2 - b3 - p0 >= -4
    /\
    - a0 - a2 + a3 + b0 + b3 - p0 >= -2
    /\
    a0 + a2 + a3 + b0 + b3 - p0 >= 0
    /\
    a1 + a2 + a3 + b2 + b3 - p0 >= 0
    /\
    a1 - a3 + b0 + b1 + b3 + p0 >= 0
    /\
    - a1 - a2 - a3 + b0 - b1 - p1 >= -4
    /\
    a0 - a2 - a3 + b0 + b1 - p1 >= -2
    /\
    a2 + a3 - b0 + b1 - b2 - p1 >= -2
    /\
    a0 + a1 - a3 - b1 - b3 - p1 >= -3
    /\
    a1 - a2 + b0 - b2 - b3 - p1 >= -3
    /\
    a1 - b0 + b1 - b2 - b3 - p1 >= -3
    /\
    a0 - a1 - a2 + b2 - b3 - p1 >= -3
    /\
    - a0 + a2 - b1 + b2 - b3 - p1 >= -3
    /\
    a1 - b0 - b1 + b2 - b3 - p1 >= -3
    /\
    a0 - a1 - b1 - b2 + b3 - p1 >= -3
    /\
    - a0 + a1 - b1 - b2 + b3 - p1 >= -3
    /\
    - a0 - a2 + a3 - b0 - b3 + p1 >= -3
    /\
    a0 - a1 - a2 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 - b0 - b1 + b3 >= -5
    /\
    - a0 - a1 - a2 - a3 + b0 + b1 - p0 >= -4
    /\
    a0 + a1 + a2 - a3 + b0 + b1 - p0 >= -1
    /\
    a0 + a1 - a2 - a3 - b1 - b2 - p0 >= -4
    /\
    - a0 - a1 + a2 - a3 - b1 - b2 - p0 >= -5
    /\
    - a0 - a1 + a2 - a3 + b1 + b2 - p0 >= -3
    /\
    a0 - a1 + a2 - a3 - b1 - b3 - p0 >= -4
    /\
    - a0 - a2 + a3 + b0 - b1 - b3 - p1 >= -4
    /\
    - a0 - a1 - b0 + b1 - b2 + b3 - p1 >= -4
    /\
    a0 + a1 - a2 - a3 + b1 + b2 + p1 >= -1
    /\
    a1 + a2 + b0 + b1 - b2 - b3 + p1 >= -1
    /\
    a1 + a2 - b0 - b1 + b2 + b3 + p1 >= -1
    /\
    - a0 + a1 - a2 - a3 + b0 - b1 - b3 - p0 >= -5
    /\
    a0 - a2 - a3 - b0 + b1 + b2 - b3 - p0 >= -4
    /\
    - a0 + a1 - a3 + b0 - b1 + b2 - b3 + p1 >= -3
    /\
    a0 - a2 - a3 + b0 - b1 - b2 + b3 + p1 >= -3
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 + b2 - b3 - p0 >= -4
    /\
    a0 - a1 - a2 - a3 - b0 + b1 - b2 + b3 - p0 >= -5
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b2 + b3 - p0 >= -4
    /\
    a0 - a1 - a2 - a3 - b0 - b1 + b2 + b3 - p0 >= -5
    /\
    a0 - a1 - a3 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a2 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a1 + a3 - b0 - b2 >= -1
    /\
    - a0 - a1 + a3 + b0 + b2 >= -1
    /\
    - a1 - a2 + a3 + b2 + b3 >= -1
    /\
    a1 + a2 + a3 - b2 - b3 >= -1
    /\
    a0 + a2 + a3 - b0 - b3 >= -1;

predicate sbox_lin_star(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, 
                        var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, 
                        var 0..1: d) = 
    - a0 + b0 + b1 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a2 + a3 - b0 + b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 - b0 - b1 - b2 >= -6
    /\
    - a0 - a1 + a2 + a3 - b0 + b1 - b2 >= -3
    /\
    - a0 + a1 - a2 + a3 - b0 + b1 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 + b0 + b1 - b3 >= -4
    /\
    - a0 - a1 - a3 - b0 - b1 - b2 + b3 >= -5
    /\
    - a0 - a1 - a2 + a3 + b1 - b2 + b3 >= -3
    /\
    a0 + a1 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 + a1 - a2 + a3 + b0 + b2 + b3 >= 0
    /\
    - a0 - a1 - a2 - b0 - b1 + b2 + b3 >= -4
    /\
    a1 + a2 + a3 + b0 - b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a2 - a3 - b0 + b1 + b2 - b3 >= -4
    /\
    - a0 + a1 + a2 - a3 + b0 - b1 + b2 - b3 >= -3
    /\
    a0 - a1 - a2 - a3 + b0 - b1 - b2 + b3 >= -4
    /\
    - a0 - a1 + a2 - a3 + b1 + b2 - b3 >= -3
    /\
    a1 + a2 - a3 - b0 - b1 + b2 + b3 >= -2
    /\
    - a0 + a1 - a2 - b0 - b1 + b2 - b3 >= -4
    /\
    a0 - a1 - a3 + b0 - b1 + b2 - b3 >= -3
    /\
    a0 + a1 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a1 + a2 - a3 + b0 - b1 - b2 >= -4
    /\
    - a0 + a1 + a2 + b0 + b1 - b2 - b3 >= -2
    /\
    - a0 - a1 - a2 + a3 - b1 - b3 >= -4
    /\
    - a0 + a2 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    a0 - a1 - a2 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b0 - b2 + b3 >= -2
    /\
    a1 - a2 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    a0 + a2 + a3 + b0 - b1 - b2 + b3 >= -1
    /\
    a1 - a2 - a3 + b1 + b2 + b3 >= -1
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a1 - a2 - a3 + b1 + b3 >= -2
    /\
    a1 + a2 - b0 - b1 - b2 - b3 >= -3
    /\
    a0 - b0 - b1 - b2 - b3 >= -3
    /\
    a0 + a1 + a3 - b0 - b2 >= -1
    /\
    a0 - a1 + a2 - a3 + b1 + b3 >= -1
    /\
    a0 + a1 - b0 - b1 - b2 >= -2
    /\
    - a1 - a2 + a3 + b2 + b3 >= -1
    /\
    - a0 - a1 + a3 + b0 + b2 >= -1
    /\
    a0 + a2 + a3 - b0 - b3 >= -1
    /\
    a0 + a1 + a2 - b0 - b1 >= -1
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a1 + a2 + a3 - b2 - b3 >= -1
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1]) then d = 1
    else d = 0 endif;

predicate sbox_diff_m(var -1..1: a0, var -1..1: a1, var -1..1: a2, var -1..1: a3, 
                      var -1..1: b0, var -1..1: b1, var -1..1: b2, var -1..1: b3,
                      var 0..1: d) = 
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0)
    else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1)
    endif
    /\
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then d = 0 else d = 1 endif;

predicate sbox_lin_m(var -1..1: a0, var -1..1: a1, var -1..1: a2, var -1..1: a3, 
                     var -1..1: b0, var -1..1: b1, var -1..1: b2, var -1..1: b3,
                     var 0..1: d) = 
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0)
    elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1)
    else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1)
    endif
    /\
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then d = 0 else d = 1 endif;
