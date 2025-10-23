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

int: RB;
int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: WU;
int: WM;
int: WL;
int: offset;

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

array[0..1, 0..127] of int: inv_bit_permutation = array2d(0..1, 0..127,
    [37, 70, 67, 77, 14, 96, 0, 126, 40, 59, 18, 26, 91, 48, 7, 63, 127, 116, 60, 101, 64, 25, 108, 122, 34, 16, 74, 69, 6, 98, 92, 82, 95, 103, 112, 15, 8, 39, 29, 123, 47, 17, 23, 99, 105, 65, 1, 43, 107, 38, 31, 79, 5, 85, 114, 56, 46, 117, 53, 62, 71, 109, 2, 88, 28, 35, 41, 27, 110, 75, 4, 113, 10, 104, 21, 66, 19, 68, 50, 102, 118, 86, 83, 11, 78, 32, 51, 20, 24, 57, 125, 100, 52, 90, 93, 13, 120, 58, 106, 36, 45, 55, 87, 119, 115, 111, 12, 73, 94, 33, 81, 97, 30, 84, 42, 72, 61, 49, 89, 76, 22, 121, 124, 80, 54, 9, 3, 44,
    126, 33, 24, 88, 70, 110, 54, 116, 104, 8, 120, 16, 45, 97, 49, 6, 31, 51, 65, 62, 0, 12, 79, 92, 74, 101, 23, 40, 107, 56, 122, 98, 10, 35, 119, 5, 39, 42, 108, 87, 99, 89, 93, 115, 25, 75, 125, 81, 90, 73, 64, 100, 37, 61, 83, 43, 22, 76, 118, 19, 84, 71, 3, 121, 47, 102, 7, 68, 30, 114, 85, 20, 60, 67, 2, 124, 38, 117, 77, 59, 36, 66, 44, 13, 32, 9, 57, 95, 69, 109, 21, 50, 127, 55, 41, 123, 91, 86, 17, 63, 78, 105, 34, 26, 111, 48, 15, 29, 52, 112, 80, 72, 46, 94, 28, 18, 82, 11, 106, 4, 113, 27, 1, 58, 53, 103, 96, 14]);

array[0..1, 0..31] of int: nibble_permutation = array2d(0..1, 0..31, 
    [10, 27, 5, 1, 30, 23, 16, 13, 21, 31, 6, 14, 0, 25, 11, 18, 15, 28, 19, 24, 7, 8, 22, 3, 4, 29, 9, 2, 26, 20, 12, 17,
     26, 13, 7, 11, 29, 0, 17, 21, 23, 5, 18, 25, 12, 10, 28, 2, 14, 19, 24, 22, 1, 8, 4, 31, 15, 6, 27, 9, 16, 30, 20, 3]);


array[0..1, 0..127] of int: ksch_bit_permutation = array2d(0..1, 0..127, 
    [0,53,87,73,22,95,99,48,61,36,108,1,24,67,119,93,54,103,69,112,16,111,94,122,31,66,33,83,47,3,65,62,123,9,101,19,5,58,89,37,38,51,28,106,82,76,121,4,70,7,42,92,104,80,45,75,114,17,2,97,46,107,63,18,109,15,127,43,13,59,29,125,77,11,50,30,12,90,118,64,20,35,57,10,124,56,68,91,116,21,84,98,52,81,126,34,105,27,120,74,6,85,40,72,113,41,23,49,79,55,102,8,117,39,88,26,25,110,14,32,115,100,86,71,78,44,96,60,
     76,30,53,35,31,46,2,79,11,125,110,87,39,91,14,101,97,118,36,48,29,80,57,115,49,18,74,85,61,82,105,126,70,12,47,111,51,17,66,1,60,96,116,71,81,114,104,15,42,124,100,4,113,44,75,89,23,0,84,107,32,26,88,8,69,121,38,94,37,86,54,21,62,123,41,10,16,95,117,65,45,50,72,20,109,58,7,67,108,28,3,55,92,103,24,5,77,9,27,102,122,6,106,22,99,34,90,56,43,83,120,64,78,59,119,93,40,98,52,68,112,33,63,25,19,73,127,13]);

array[0..1, 0..127] of int: inv_ksch_bit_permutation = array2d(0..1, 0..127,
    [0, 11, 58, 29, 47, 36, 100, 49, 111, 33, 83, 73, 76, 68, 118, 65, 20, 57, 63, 35, 80, 89, 4, 106, 12, 116, 115, 97, 42, 70, 75, 24, 119, 26, 95, 81, 9, 39, 40, 113, 102, 105, 50, 67, 125, 54, 60, 28, 7, 107, 74, 41, 92, 1, 16, 109, 85, 82, 37, 69, 127, 8, 31, 62, 79, 30, 25, 13, 86, 18, 48, 123, 103, 3, 99, 55, 45, 72, 124, 108, 53, 93, 44, 27, 90, 101, 122, 2, 114, 38, 77, 87, 51, 15, 22, 5, 126, 59, 91, 6, 121, 34, 110, 17, 52, 96, 43, 61, 10, 64, 117, 21, 19, 104, 56, 120, 88, 112, 78, 14, 98, 46, 23, 32, 84, 71, 94, 66,
     57, 39, 6, 90, 51, 95, 101, 86, 63, 97, 75, 8, 33, 127, 14, 47, 76, 37, 25, 124, 83, 71, 103, 56, 94, 123, 61, 98, 89, 20, 1, 4, 60, 121, 105, 3, 18, 68, 66, 12, 116, 74, 48, 108, 53, 80, 5, 34, 19, 24, 81, 36, 118, 2, 70, 91, 107, 22, 85, 113, 40, 28, 72, 122, 111, 79, 38, 87, 119, 64, 32, 43, 82, 125, 26, 54, 0, 96, 112, 7, 21, 44, 29, 109, 58, 27, 69, 11, 62, 55, 106, 13, 92, 115, 67, 77, 41, 16, 117, 104, 50, 15, 99, 93, 46, 30, 102, 59, 88, 84, 10, 35, 120, 52, 45, 23, 42, 78, 17, 114, 110, 65, 100, 73, 49, 9, 31, 126]);

array[0..offset,0..1, 0..127] of var int: inv_kperm_at_round;
constraint forall (b in 0..1, i in 0..127) (inv_kperm_at_round[0, b, i] = inv_ksch_bit_permutation[b, i]);
constraint forall (n in 1..offset, b in 0..1, i in 0..127) (inv_kperm_at_round[n, b, i] = inv_ksch_bit_permutation[b, inv_kperm_at_round[n - 1, b, i]]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____                                 _      ____                                                _          _    _               __  __              _                _  __            
% |  _ \  _ __  ___ __   __ ___  _ __  | |_   / ___| ___   _ __ ___   _ __ ___    ___   _ __      / \    ___ | |_ (_)__   __ ___  |  \/  |  __ _  ___ | |_  ___  _ __  | |/ / ___  _   _ 
% | |_) || '__|/ _ \\ \ / // _ \| '_ \ | __| | |    / _ \ | '_ ` _ \ | '_ ` _ \  / _ \ | '_ \    / _ \  / __|| __|| |\ \ / // _ \ | |\/| | / _` |/ __|| __|/ _ \| '__| | ' / / _ \| | | |
% |  __/ | |  |  __/ \ V /|  __/| | | || |_  | |___| (_) || | | | | || | | | | || (_) || | | |  / ___ \| (__ | |_ | | \ V /|  __/ | |  | || (_| |\__ \| |_|  __/| |    | . \|  __/| |_| |
% |_|    |_|   \___|  \_/  \___||_| |_| \__|  \____|\___/ |_| |_| |_||_| |_| |_| \___/ |_| |_| /_/   \_\\___| \__||_|  \_/  \___| |_|  |_| \__,_||___/ \__|\___||_|    |_|\_\\___| \__, |
%                                                                                                                                                                                  |___/ 
%  ____   _  _        
% | __ ) (_)| |_  ___ 
% |  _ \ | || __|/ __|
% | |_) || || |_ \__ \
% |____/ |_| \__||___/
% Prevent common active master key bits
array[0..127] of var 0..1: k0;
array[0..127] of var 0..1: k1;

constraint forall(i in 0..127)(k0[inv_kperm_at_round[offset, 0, i]] = bool2int(xb0[0, i] != 0));
constraint forall(i in 0..127)(k1[inv_kperm_at_round[offset, 1, i]] = bool2int(xb1[0, i] != 0));
% constraint sum(i in 0..127)(bool2int(k0[i] + k1[i] == 2)) >= 3;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  ____  
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| __ ) 
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  |  _ \ 
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_) |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||____/ 
% Extend EU backward by one linear layer

array[0..RB, 0..127] of var -1..1: xb0;
array[0..(RB - 1), 0..127] of var -1..1: yb0;
array[0..(RB - 1), 0..127] of var -1..1: zb0;
array[0..(RB - 1), 0..31] of var 0..1: ab0;

% Link the end of EB to the output 
constraint forall(i in 0..127)(xb0[RB, i] = xu0[0, i]);

% S-box layer
constraint forall(round in 0..(RB - 1), i in 0..31)
(
    sbox_diff_m(yb0[round, 4*i], yb0[round, 4*i + 1], yb0[round, 4*i + 2], yb0[round, 4*i + 3], 
                xb0[round, 4*i], xb0[round, 4*i + 1], xb0[round, 4*i + 2], xb0[round, 4*i + 3],
                ab0[round, i])
);

% Permutation layer
constraint forall(round in 0..(RB - 1), i in 0..31, j in 0..3)
(
    if (offset + round) < 4 then 
    (
        zb0[round, bit_permutation[0, 4*i + j]] = yb0[round, 4*i + j]
    ) else 
    (
        zb0[round, 4*nibble_permutation[0, i] + j] = yb0[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RB - 1), i in 0..7, j in 0..3)
(
    xor3_pr1(xb0[round + 1, 16*i + 4*1 + j], xb0[round + 1, 16*i + 4*2 + j], xb0[round + 1, 16*i + 4*3 + j], zb0[round, 16*i + 4*0 + j])
    /\
    xor3_pr1(xb0[round + 1, 16*i + 4*0 + j], xb0[round + 1, 16*i + 4*2 + j], xb0[round + 1, 16*i + 4*3 + j], zb0[round, 16*i + 4*1 + j])
    /\
    xor3_pr1(xb0[round + 1, 16*i + 4*0 + j], xb0[round + 1, 16*i + 4*1 + j], xb0[round + 1, 16*i + 4*3 + j], zb0[round, 16*i + 4*2 + j])
    /\
    xor3_pr1(xb0[round + 1, 16*i + 4*0 + j], xb0[round + 1, 16*i + 4*1 + j], xb0[round + 1, 16*i + 4*2 + j], zb0[round, 16*i + 4*3 + j])
);

% #############################################################################################################################################

array[0..RB, 0..127] of var -1..1: xb1;
array[0..(RB - 1), 0..127] of var -1..1: yb1;
array[0..(RB - 1), 0..127] of var -1..1: zb1;
array[0..(RB - 1), 0..31] of var 0..1: ab1;

% Link the end of EB to the output
constraint forall(i in 0..127)(xb1[RB, i] = xu1[0, i]);

% S-box layer
constraint forall(round in 0..(RB - 1), i in 0..31)
(
    sbox_diff_m(yb1[round, 4*i], yb1[round, 4*i + 1], yb1[round, 4*i + 2], yb1[round, 4*i + 3], 
                xb1[round, 4*i], xb1[round, 4*i + 1], xb1[round, 4*i + 2], xb1[round, 4*i + 3],
                ab1[round, i])
);

% Permutation layer
constraint forall(round in 0..(RB - 1), i in 0..31, j in 0..3)
(
    if (offset + round) < 4 then 
    (
        zb1[round, bit_permutation[1, 4*i + j]] = yb1[round, 4*i + j]
    ) else 
    (
        zb1[round, 4*nibble_permutation[1, i] + j] = yb1[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RB - 1), i in 0..7, j in 0..3)
(
    xor3_pr1(xb1[round + 1, 16*i + 4*1 + j], xb1[round + 1, 16*i + 4*2 + j], xb1[round + 1, 16*i + 4*3 + j], zb1[round, 16*i + 4*0 + j])
    /\
    xor3_pr1(xb1[round + 1, 16*i + 4*0 + j], xb1[round + 1, 16*i + 4*2 + j], xb1[round + 1, 16*i + 4*3 + j], zb1[round, 16*i + 4*1 + j])
    /\
    xor3_pr1(xb1[round + 1, 16*i + 4*0 + j], xb1[round + 1, 16*i + 4*1 + j], xb1[round + 1, 16*i + 4*3 + j], zb1[round, 16*i + 4*2 + j])
    /\
    xor3_pr1(xb1[round + 1, 16*i + 4*0 + j], xb1[round + 1, 16*i + 4*1 + j], xb1[round + 1, 16*i + 4*2 + j], zb1[round, 16*i + 4*3 + j])
);

% the activeness pattern of ab0 and ab1 must be the same
constraint forall(i in 0..31)(ab0[0, i] = ab1[0, i]);

%     _           _    _   __  _        _         _    ____                    _                _         _        
%    / \    _ __ | |_ (_) / _|(_)  ___ (_)  __ _ | |  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___ 
%   / _ \  | '__|| __|| || |_ | | / __|| | / _` || | | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __|
%  / ___ \ | |   | |_ | ||  _|| || (__ | || (_| || | | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \
% /_/   \_\|_|    \__||_||_|  |_| \___||_| \__,_||_|  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/
                                                                                                                  
% some artificial constraints to make the two paths compatible

constraint sum(i in 0..31)(ab0[0, i]) == 2;
% constraint sum(i in 0..31)(dml0[RM - 1, i]) = 1;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..127] of var 0..1: xu0;
array[0..(RU - 1), 0..127] of var 0..1: yu0;
array[0..(RU - 1), 0..127] of var 0..1: zu0;
array[0..(RU - 1), 0..31, 0..1] of var 0..1: pu0;

% S-box layer
constraint forall(round in 0..(RU - 1), i in 0..31)
(
    sbox_diff(xu0[round, 4*i], xu0[round, 4*i + 1], xu0[round, 4*i + 2], xu0[round, 4*i + 3], 
              yu0[round, 4*i], yu0[round, 4*i + 1], yu0[round, 4*i + 2], yu0[round, 4*i + 3],
              pu0[round, i, 0], pu0[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RU - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + round) < 4 then
    (
        zu0[round, bit_permutation[0, 4*i + j]] = yu0[round, 4*i + j]
    ) else
    (
        zu0[round, 4*nibble_permutation[0, i] + j] = yu0[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RU - 1), i in 0..7, j in 0..3)
(
    xor3(zu0[round, 16*i + 4*1 + j], zu0[round, 16*i + 4*2 + j], zu0[round, 16*i + 4*3 + j], xu0[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zu0[round, 16*i + 4*0 + j], zu0[round, 16*i + 4*2 + j], zu0[round, 16*i + 4*3 + j], xu0[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zu0[round, 16*i + 4*0 + j], zu0[round, 16*i + 4*1 + j], zu0[round, 16*i + 4*3 + j], xu0[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zu0[round, 16*i + 4*0 + j], zu0[round, 16*i + 4*1 + j], zu0[round, 16*i + 4*2 + j], xu0[round  +1, 16*i + 4*3 + j])
);

% Exclude the all-zero input/output difference 
constraint sum(i in 0..127) (xu0[0, i]) != 0;

% Compute the probability of the differential transition
var 0..128: PU0;
constraint PU0 = sum(r in 0..(RU - 1), i in 0..31)(3*pu0[r, i, 0] + 2*pu0[r, i, 1]);

% #############################################################################################################################################

array[0..RU, 0..127] of var 0..1: xu1;
array[0..(RU - 1), 0..127] of var 0..1: yu1;
array[0..(RU - 1), 0..127] of var 0..1: zu1;
array[0..(RU - 1), 0..31, 0..1] of var 0..1: pu1;

% S-box layer
constraint forall(round in 0..(RU - 1), i in 0..31)
(
    sbox_diff(xu1[round, 4*i], xu1[round, 4*i + 1], xu1[round, 4*i + 2], xu1[round, 4*i + 3], 
              yu1[round, 4*i], yu1[round, 4*i + 1], yu1[round, 4*i + 2], yu1[round, 4*i + 3],
              pu1[round, i, 0], pu1[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RU - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + round) < 4 then 
    (
        zu1[round, bit_permutation[1, 4*i + j]] = yu1[round, 4*i + j]
    ) else 
    (
        zu1[round, 4*nibble_permutation[1, i] + j] = yu1[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RU - 1), i in 0..7, j in 0..3)
(
    xor3(zu1[round, 16*i + 4*1 + j], zu1[round, 16*i + 4*2 + j], zu1[round, 16*i + 4*3 + j], xu1[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zu1[round, 16*i + 4*0 + j], zu1[round, 16*i + 4*2 + j], zu1[round, 16*i + 4*3 + j], xu1[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zu1[round, 16*i + 4*0 + j], zu1[round, 16*i + 4*1 + j], zu1[round, 16*i + 4*3 + j], xu1[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zu1[round, 16*i + 4*0 + j], zu1[round, 16*i + 4*1 + j], zu1[round, 16*i + 4*2 + j], xu1[round  +1, 16*i + 4*3 + j])
);

% Exclude the all-zero input/output difference 
constraint sum(i in 0..127) (xu1[0, i]) != 0;

% Compute the probability of the differential transition
var 0..128: PU1;
constraint PU1 = sum(r in 0..(RU - 1), i in 0..31)(3*pu1[r, i, 0] + 2*pu1[r, i, 1]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU                                                                                                                                                             

array[0..RM, 0..127] of var -1..1: xmu0;
array[0..(RM - 1), 0..127] of var -1..1: ymu0;
array[0..(RM - 1), 0..127] of var -1..1: zmu0;
array[0..(RM - 1), 0..31] of var 0..1: dmu0;

% link the input of EM to the output of EU
constraint forall(i in 0..127)(xmu0[0, i] = xu0[RU, i]);

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu0[round, 4*i], xmu0[round, 4*i + 1], xmu0[round, 4*i + 2], xmu0[round, 4*i + 3], 
                       ymu0[round, 4*i], ymu0[round, 4*i + 1], ymu0[round, 4*i + 2], ymu0[round, 4*i + 3],
                       dmu0[round, i])
    ) else 
    (
        sbox_diff_m(xmu0[round, 4*i], xmu0[round, 4*i + 1], xmu0[round, 4*i + 2], xmu0[round, 4*i + 3], 
                    ymu0[round, 4*i], ymu0[round, 4*i + 1], ymu0[round, 4*i + 2], ymu0[round, 4*i + 3],
                    dmu0[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + round) < 4 then 
    (
        zmu0[round, bit_permutation[0, 4*i + j]] = ymu0[round, 4*i + j]
    ) else 
    (
        zmu0[round, 4*nibble_permutation[0, i] + j] = ymu0[round, 4*i + j]
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < RMU then 
    (
        xor3(zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*2 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*0 + j])
        /\
        xor3(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*2 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*1 + j])
        /\
        xor3(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*2 + j])
        /\
        xor3(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*2 + j], xmu0[round + 1, 16*i + 4*3 + j])
    ) else 
    (
        xor3_pr1(zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*2 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*0 + j])
        /\
        xor3_pr1(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*2 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*1 + j])
        /\
        xor3_pr1(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*3 + j], xmu0[round + 1, 16*i + 4*2 + j])
        /\
        xor3_pr1(zmu0[round, 16*i + 4*0 + j], zmu0[round, 16*i + 4*1 + j], zmu0[round, 16*i + 4*2 + j], xmu0[round + 1, 16*i + 4*3 + j])
    ) endif
    
);

% #############################################################################################################################################

array[0..RM, 0..127] of var -1..1: xmu1;
array[0..(RM - 1), 0..127] of var -1..1: ymu1;
array[0..(RM - 1), 0..127] of var -1..1: zmu1;
array[0..(RM - 1), 0..31] of var 0..1: dmu1;

% link the input of EM to the output of EU
constraint forall(i in 0..127)(xmu1[0, i] = xu1[RU, i]);

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu1[round, 4*i], xmu1[round, 4*i + 1], xmu1[round, 4*i + 2], xmu1[round, 4*i + 3], 
                       ymu1[round, 4*i], ymu1[round, 4*i + 1], ymu1[round, 4*i + 2], ymu1[round, 4*i + 3],
                       dmu1[round, i])
    ) else 
    (
        sbox_diff_m(xmu1[round, 4*i], xmu1[round, 4*i + 1], xmu1[round, 4*i + 2], xmu1[round, 4*i + 3], 
                    ymu1[round, 4*i], ymu1[round, 4*i + 1], ymu1[round, 4*i + 2], ymu1[round, 4*i + 3],
                    dmu1[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + round) < 4 then 
    (
        zmu1[round, bit_permutation[1, 4*i + j]] = ymu1[round, 4*i + j]
    ) else 
    (
        zmu1[round, 4*nibble_permutation[1, i] + j] = ymu1[round, 4*i + j]
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < RMU then 
    (
        xor3(zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*2 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*0 + j])
        /\
        xor3(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*2 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*1 + j])
        /\
        xor3(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*2 + j])
        /\
        xor3(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*2 + j], xmu1[round + 1, 16*i + 4*3 + j])
    ) else 
    (
        xor3_pr1(zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*2 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*0 + j])
        /\
        xor3_pr1(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*2 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*1 + j])
        /\
        xor3_pr1(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*3 + j], xmu1[round + 1, 16*i + 4*2 + j])
        /\
        xor3_pr1(zmu1[round, 16*i + 4*0 + j], zmu1[round, 16*i + 4*1 + j], zmu1[round, 16*i + 4*2 + j], xmu1[round + 1, 16*i + 4*3 + j])
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

array[0..RM, 0..127] of var -1..1: xml0;
array[0..(RM - 1), 0..127] of var -1..1: yml0;
array[0..(RM - 1), 0..127] of var -1..1: zml0;
array[0..(RM - 1), 0..31] of var 0..1: dml0;

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml0[round, 4*i], yml0[round, 4*i + 1], yml0[round, 4*i + 2], yml0[round, 4*i + 3], 
                   xml0[round, 4*i], xml0[round, 4*i + 1], xml0[round, 4*i + 2], xml0[round, 4*i + 3],
                   dml0[round, i])
    ) else
    (
        sbox_lin_star(yml0[round, 4*i], yml0[round, 4*i + 1], yml0[round, 4*i + 2], yml0[round, 4*i + 3], 
                      xml0[round, 4*i], xml0[round, 4*i + 1], xml0[round, 4*i + 2], xml0[round, 4*i + 3],
                      dml0[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + round) < 4 then 
    (
        zml0[round, bit_permutation[0, 4*i + j]] = yml0[round, 4*i + j]
    ) else 
    (
        zml0[round, 4*nibble_permutation[0, i] + j] = yml0[round, 4*i + j]
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < (RM - RML) then 
    (
        xor3_pr1(xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*2 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*0 + j])
        /\
        xor3_pr1(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*2 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*1 + j])
        /\
        xor3_pr1(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*2 + j])
        /\
        xor3_pr1(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*2 + j], zml0[round, 16*i + 4*3 + j])
    ) else
    (
        xor3(xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*2 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*0 + j])
        /\
        xor3(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*2 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*1 + j])
        /\
        xor3(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*3 + j], zml0[round, 16*i + 4*2 + j])
        /\
        xor3(xml0[round + 1, 16*i + 4*0 + j], xml0[round + 1, 16*i + 4*1 + j], xml0[round + 1, 16*i + 4*2 + j], zml0[round, 16*i + 4*3 + j])
    ) endif
);

% #############################################################################################################################################

array[0..RM, 0..127] of var -1..1: xml1;
array[0..(RM - 1), 0..127] of var -1..1: yml1;
array[0..(RM - 1), 0..127] of var -1..1: zml1;
array[0..(RM - 1), 0..31] of var 0..1: dml1;

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..31)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml1[round, 4*i], yml1[round, 4*i + 1], yml1[round, 4*i + 2], yml1[round, 4*i + 3], 
                   xml1[round, 4*i], xml1[round, 4*i + 1], xml1[round, 4*i + 2], xml1[round, 4*i + 3],
                   dml1[round, i])
    ) else
    (
        sbox_lin_star(yml1[round, 4*i], yml1[round, 4*i + 1], yml1[round, 4*i + 2], yml1[round, 4*i + 3], 
                      xml1[round, 4*i], xml1[round, 4*i + 1], xml1[round, 4*i + 2], xml1[round, 4*i + 3],
                      dml1[round, i])
    ) endif
);

% Permutation layer
constraint forall(round in 0..(RM - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + round) < 4 then 
    (
        zml1[round, bit_permutation[1, 4*i + j]] = yml1[round, 4*i + j]
    ) else 
    (
        zml1[round, 4*nibble_permutation[1, i] + j] = yml1[round, 4*i + j]
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), i in 0..7, j in 0..3)
(
    if round < (RM - RML) then 
    (
        xor3_pr1(xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*2 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*0 + j])
        /\
        xor3_pr1(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*2 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*1 + j])
        /\
        xor3_pr1(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*2 + j])
        /\
        xor3_pr1(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*2 + j], zml1[round, 16*i + 4*3 + j])
    ) else
    (
        xor3(xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*2 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*0 + j])
        /\
        xor3(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*2 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*1 + j])
        /\
        xor3(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*3 + j], zml1[round, 16*i + 4*2 + j])
        /\
        xor3(xml1[round + 1, 16*i + 4*0 + j], xml1[round + 1, 16*i + 4*1 + j], xml1[round + 1, 16*i + 4*2 + j], zml1[round, 16*i + 4*3 + j])
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

array[0..RL, 0..127] of var 0..1: xl0;
array[0..(RL - 1), 0..127] of var 0..1: yl0;
array[0..(RL - 1), 0..127] of var 0..1: zl0;
array[0..(RL - 1), 0..31, 0..1] of var 0..1: cl0;

% link the input of EL to the output of EML
constraint forall(i in 0..127)(xl0[0, i] = xml0[RM, i]);

% S-box layer
constraint forall(round in 0..(RL - 1), i in 0..31)
(
    sbox_lin(yl0[round, 4*i], yl0[round, 4*i + 1], yl0[round, 4*i + 2], yl0[round, 4*i + 3],
             xl0[round, 4*i], xl0[round, 4*i + 1], xl0[round, 4*i + 2], xl0[round, 4*i + 3],
             cl0[round, i, 0], cl0[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RL - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + RM + round) < 4 then 
    (
        zl0[round, bit_permutation[0, 4*i + j]] = yl0[round, 4*i + j]
    ) else 
    (
        zl0[round, 4*nibble_permutation[0, i] + j] = yl0[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RL - 1), i in 0..7, j in 0..3)
(
    xor3(zl0[round, 16*i + 4*1 + j], zl0[round, 16*i + 4*2 + j], zl0[round, 16*i + 4*3 + j], xl0[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zl0[round, 16*i + 4*0 + j], zl0[round, 16*i + 4*2 + j], zl0[round, 16*i + 4*3 + j], xl0[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zl0[round, 16*i + 4*0 + j], zl0[round, 16*i + 4*1 + j], zl0[round, 16*i + 4*3 + j], xl0[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zl0[round, 16*i + 4*0 + j], zl0[round, 16*i + 4*1 + j], zl0[round, 16*i + 4*2 + j], xl0[round + 1, 16*i + 4*3 + j])
);

% Exclude the all-zero output mask
constraint sum(i in 0..127) (xl0[RL, i]) != 0;

% Compute the squared correlation of the linear trail
var 0..128: CL0;
constraint CL0 = sum(r in 0..(RL - 1), i in 0..31)(4*cl0[r, i, 0] + 2*cl0[r, i, 1]);

% #############################################################################################################################################
array[0..RL, 0..127] of var 0..1: xl1;
array[0..(RL - 1), 0..127] of var 0..1: yl1;
array[0..(RL - 1), 0..127] of var 0..1: zl1;
array[0..(RL - 1), 0..31, 0..1] of var 0..1: cl1;

% link the input of EL to the output of EML
constraint forall(i in 0..127)(xl1[0, i] = xml1[RM, i]);

% S-box layer
constraint forall(round in 0..(RL - 1), i in 0..31)
(
    sbox_lin(yl1[round, 4*i], yl1[round, 4*i + 1], yl1[round, 4*i + 2], yl1[round, 4*i + 3],
             xl1[round, 4*i], xl1[round, 4*i + 1], xl1[round, 4*i + 2], xl1[round, 4*i + 3], 
             cl1[round, i, 0], cl1[round, i, 1])
);

% Permutation layer
constraint forall(round in 0..(RL - 1), i in 0..31, j in 0..3)
(
    if (offset + RB + RU + RM + round) < 4 then 
    (
        zl1[round, bit_permutation[1, 4*i + j]] = yl1[round, 4*i + j]
    ) else 
    (
        zl1[round, 4*nibble_permutation[1, i] + j] = yl1[round, 4*i + j]
    ) endif
);

% MixColumns layer
constraint forall(round in 0..(RL - 1), i in 0..7, j in 0..3)
(
    xor3(zl1[round, 16*i + 4*1 + j], zl1[round, 16*i + 4*2 + j], zl1[round, 16*i + 4*3 + j], xl1[round + 1, 16*i + 4*0 + j])
    /\
    xor3(zl1[round, 16*i + 4*0 + j], zl1[round, 16*i + 4*2 + j], zl1[round, 16*i + 4*3 + j], xl1[round + 1, 16*i + 4*1 + j])
    /\
    xor3(zl1[round, 16*i + 4*0 + j], zl1[round, 16*i + 4*1 + j], zl1[round, 16*i + 4*3 + j], xl1[round + 1, 16*i + 4*2 + j])
    /\
    xor3(zl1[round, 16*i + 4*0 + j], zl1[round, 16*i + 4*1 + j], zl1[round, 16*i + 4*2 + j], xl1[round + 1, 16*i + 4*3 + j])
);

% Exclude the all-zero output mask
constraint sum(i in 0..127) (xl1[RL, i]) != 0;

% Compute the squared correlation of the linear trail
var 0..128: CL1;
constraint CL1 = sum(r in 0..(RL - 1), i in 0..31)(4*cl1[r, i, 0] + 2*cl1[r, i, 1]);

% output masks of branch 0 and branch 1 must be equal
constraint forall(i in 0..127)(xl0[RL, i] = xl1[RL, i]);

% assign a one-dimensional array to xl0[RL]
% constraint yml0[RM - 1, 0..127] = array1d(0..127, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);


% #############################################################################################################################################
% #############################################################################################################################################

var 0..RM*(128): NASM0;
var 0..RM*(128): CM0;
constraint CM0 = 
                1*sum(r in 0..(RM - 1), i in 0..127)(
                bool2int(xmu0[r, i] = -1 /\ xml0[r, i] = -1)
                +
                bool2int(xmu0[r, i] = -1 /\ xml0[r, i] = 1)
                +                 
                bool2int(xmu0[r, i] = 1 /\ xml0[r, i] = -1));

constraint NASM0 = sum(r in 0..(RM - 1), i in 0..31)(bool2int(dmu0[r, i] + dml0[r, i] == 2));

% #############################################################################################################################################

var 0..RM*(128): NASM1;
var 0..RM*(128): CM1;
constraint CM1 = 
                1*sum(r in 0..(RM - 1), i in 0..127)(
                bool2int(xmu1[r, i] = -1 /\ xml1[r, i] = -1)
                +
                bool2int(xmu1[r, i] = -1 /\ xml1[r, i] = 1)
                +                 
                bool2int(xmu1[r, i] = 1 /\ xml1[r, i] = -1));

constraint NASM1 = sum(r in 0..(RM - 1), i in 0..31)(bool2int(dmu1[r, i] + dml1[r, i] == 2));

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
%     int_search([yl0[RL - 1, 0, i] | i in 0..63], indomain_min, complete),
%     int_search([yl0[RL - 1, 1, i] | i in 0..63], indomain_min, complete),
%     int_search([yl0[RL - 1, 2, i] | i in 0..63], indomain_min, complete),
%     int_search([yl0[RL - 1, 3, i] | i in 0..63], indomain_min, complete),
%     int_search([yl0[RL - 1, 3, i] | i in 0..63], indomain_min, complete)]);

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

solve minimize (WU*PU0 + WM*CM0 + WL*CL0 + WU*PU1 + WM*CM1 + WL*CL1);
% solve minimize (WU*PU0 + WM*NASM0 + WL*CL0 + 
%                 WU*PU1 + WM*NASM1 + WL*CL1);

%artificial constraints for key recovery:
% constraint (WU*PU0 + WM*CM0 + WL*CL0) <= 30; 
% constraint (WU*PU1 + WM*CM1 + WL*CL1) <= 30;
% solve minimize sum(r in 0..(RB - 1), i in 0..31)(ab0[r, i] + ab1[r, i]);

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
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0 /\ d = 0)
    else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ d = 1) endif;

predicate sbox_lin_m(var -1..1: a0, var -1..1: a1, var -1..1: a2, var -1..1: a3, 
                     var -1..1: b0, var -1..1: b1, var -1..1: b2, var -1..1: b3,
                     var 0..1: d) = 
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0)
    elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1)
    else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1)
    endif
    /\
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then d = 0 else d = 1 endif;