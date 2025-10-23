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


import sys

def trim(docstring):
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)

class DrawDL():
    """
    Draw the shape of a given differential-linear distinguisher
    """

    def __init__(self, dlobject, output_file_name="output.tex", bn=0):
        self.model = dlobject.model
        self.RB = dlobject.RB
        self.RU = dlobject.RU
        self.RM = dlobject.RM
        self.RL = dlobject.RL
        self.RD = dlobject.RU + dlobject.RM + dlobject.RL
        self.offset = dlobject.offset
        self.branchtype = bn
        self.output_file_name = output_file_name
        if bn == 0:
            self.upper_trail = dlobject.upper_trail_0
            self.lower_trail = dlobject.lower_trail_0
        else:
            self.upper_trail = dlobject.upper_trail_1
            self.lower_trail = dlobject.lower_trail_1
        self.attack_summary = dlobject.attack_summary
        self.bit_permutation = [[6,46,62,126,70,52,28,14,36,125,72,83,106,95,4,35,
                                25,41,10,76,87,74,120,42,88,21,11,67,64,38,112,50,
                                85,109,24,65,99,0,49,37,8,66,114,47,127,100,56,40,
                                13,117,78,86,92,58,124,101,55,89,97,9,18,116,59,15,
                                20,45,75,2,77,27,1,60,115,107,26,69,119,3,84,51,
                                123,110,31,82,113,53,81,102,63,118,93,12,30,94,108,32,
                                5,111,29,43,91,19,79,33,73,44,98,48,22,61,68,105,
                                34,71,54,104,17,57,80,103,96,121,23,39,122,90,7,16],
                                [20,122,74,62,119,35,15,66,9,85,32,117,21,83,127,106,
                                11,98,115,59,71,90,56,26,2,44,103,121,114,107,68,16,
                                84,1,102,33,80,52,76,36,27,94,37,55,82,12,112,64,
                                105,14,91,17,108,124,6,93,29,86,123,79,72,53,19,99,
                                50,18,81,73,67,88,4,61,111,49,24,45,57,78,100,22,
                                110,47,116,54,60,70,97,39,3,41,48,96,23,42,113,87,
                                126,13,31,40,51,25,65,125,8,101,118,28,38,89,5,104,
                                109,120,69,43,7,77,58,34,10,63,30,95,75,46,0,92]]
        
        self.nibble_permutation = [[10, 27, 5, 1, 30, 23, 16, 13, 21, 31, 6, 14, 0, 25, 11, 18, 15, 28, 19, 24, 7, 8, 22, 3, 4, 29, 9, 2, 26, 20, 12, 17],
                                   [26, 13, 7, 11, 29, 0, 17, 21, 23, 5, 18, 25, 12, 10, 28, 2, 14, 19, 24, 22, 1, 8, 4, 31, 15, 6, 27, 9, 16, 30, 20, 3]]

    def generate_distinguisher_shape(self):
        """
        Draw the figure of the Rectangle distinguisher
        """

        contents = ""
        # head lines
        contents += trim(r"""
                        \documentclass{standalone}
                        \usepackage{orthros}
                        \usepackage{tugcolors}
                        \usepackage{comment}
                        \begin{document}
                        \begin{tikzpicture}[spn]
                        \orthrosinit""") + "\n\n"
        if self.branchtype == 0:
            contents += trim(r"""\orthrosbranchonetrue % define permutation variant (branch 1=left or 2=right?)
                                 \spnlinktrue""") + "\n"
        else:
            contents += trim(r"""\orthrosbranchonefalse % define permutation variant (branch 1=left or 2=right?)
                                 \spnlinktrue""") + "\n"
        # draw EB
        for r in range(self.RB):
            if self.offset + r < 4:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundbits""") + "\n"
            else:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundnibble""") + "\n"
            x_one = [str(i) for i in range(128) if int(self.model.getVarByName(f"x{self.branchtype}_{r}_{i}_{1}").Xn) == 1]
            x_unknown = [str(i) for i in range(128) if int(self.model.getVarByName(f"x{self.branchtype}_{r}_{i}_0").Xn) == 1]
            if self.offset + r < 4:
                yz_one = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if int(self.model.getVarByName(f"y{self.branchtype}_{r}_{i}_{1}").Xn) == 1]
                yz_unknown = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if int(self.model.getVarByName(f"y{self.branchtype}_{r}_{i}_{0}").Xn) == 1]
            else:
                yz_one = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if int(self.model.getVarByName(f"y{self.branchtype}_{r}_{i}_{1}").Xn) == 1]
                yz_unknown = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if int(self.model.getVarByName(f"y{self.branchtype}_{r}_{i}_{0}").Xn) == 1]
            diff_active_sboxes = [str(i) for i in range(32) if int(self.model.getVarByName(f"active{self.branchtype}_{r}_{i}").Xn) >= 1]
            contents += trim(r"""\orthrosmarkbits[keyone]""") + "\n"
            contents += "{{{}}}".format(",".join(x_one)) + "\n"
            contents += "{{{}}}".format(",".join(diff_active_sboxes)) + "\n"
            contents += "{{{}}}".format(",".join(yz_one)) + "\n"
            contents += trim(r"""\orthrosmarkbits[keyunknown]""") + "\n"
            contents += "{{{}}}".format(",".join(x_unknown)) + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(yz_unknown)) + "\n"

        # draw EU
        for r in range(self.RU):
            if self.offset + self.RB + r < 4:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                     \orthrosroundbits""") + "\n"
            else: 
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                     \orthrosroundnibble""") + "\n"
            x_one = [str(i) for i in range(128) if self.upper_trail["x"][r][i] == 1]
            x_unknown = [str(i) for i in range(128) if self.upper_trail["x"][r][i] == -1]
            if self.offset + self.RB + r < 4:
                yz_one = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.upper_trail["y"][r][i] == 1]
                yz_unknown = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.upper_trail["y"][r][i] == -1]   
            else:
                yz_one = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.upper_trail["y"][r][i] == 1]
                yz_unknown = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.upper_trail["y"][r][i] == -1]
            diff_active_sboxes = [str(i) for i in range(32) if sum([int(self.model.getVarByName(f"pu{self.branchtype}_{r}_{i}_{j}").Xn) for j in range(2)]) >= 1]
            contents += trim(r"""\orthrosmarkbits[diffone]""")
            contents += "{{{}}}".format(",".join(x_one)) + "\n"
            contents += "{{{}}}".format(",".join(diff_active_sboxes)) + "\n"
            contents += "{{{}}}".format(",".join(yz_one)) + "\n"
            contents += trim(r"""\orthrosmarkbits[diffunknown]""")
            contents += "{{{}}}".format(",".join(x_unknown)) + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(yz_unknown)) + "\n"

        # draw EM
        for r in range(self.RM):
            if self.offset + self.RB + self.RU + r < 4:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundbits""") + "\n"
            else: 
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundnibble""") + "\n"
            x_diffone = [str(i) for i in range(128) if self.upper_trail["x"][self.RU + r][i] == 1]
            x_diffunknown = [str(i) for i in range(128) if self.upper_trail["x"][self.RU + r][i] == -1]
            x_linone = [str(i) for i in range(128) if self.lower_trail["x"][r][i] == 1]
            x_linunknown = [str(i) for i in range(128) if self.lower_trail["x"][r][i] == -1]
            if self.offset + self.RB + self.RU + r < 4:
                yz_diffone = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.upper_trail["y"][self.RU + r][i] == 1]
                yz_diffunknown = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.upper_trail["y"][self.RU + r][i] == -1]
                yz_linone = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.lower_trail["y"][r][i] == 1]
                yz_linunknown = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.lower_trail["y"][r][i] == -1]
            else:
                yz_diffone = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.upper_trail["y"][self.RU + r][i] == 1]
                yz_diffunknown = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.upper_trail["y"][self.RU + r][i] == -1]
                yz_linone = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.lower_trail["y"][r][i] == 1]
                yz_linunknown = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.lower_trail["y"][r][i] == -1]
            diff_active_sboxes = [str(i) for i in range(32) if int(self.model.getVarByName(f"active_mu{self.branchtype}_{r}_{i}").Xn) >= 1]
            lin_active_sboxes = [str(i) for i in range(32) if int(self.model.getVarByName(f"active_ml{self.branchtype}_{r}_{i}").Xn) >= 1]
            common_active_sboxes = [str(i) for i in range(32) if int(self.model.getVarByName(f"active_mu{self.branchtype}_{r}_{i}").Xn) + int(self.model.getVarByName(f"active_ml{self.branchtype}_{r}_{i}").Xn) == 2]
            contents += trim(r"""\orthrosmarkbits[diffone]""") + "\n"
            contents += "{{{}}}".format(",".join(x_diffone)) + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(yz_diffone)) + "\n"
            contents += trim(r"""\orthrosmarkbits[diffunknown]""") + "\n"
            contents += "{{{}}}".format(",".join(x_diffunknown)) + "\n"
            contents += "{{{}}}".format(",".join(diff_active_sboxes)) + "\n"
            contents += "{{{}}}".format(",".join(yz_diffunknown)) + "\n"
            contents += trim(r"""\orthrosmarkbits[linone]""") + "\n"
            contents += "{{{}}}".format(",".join(x_linone)) + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(yz_linone)) + "\n"
            contents += trim(r"""\orthrosmarkbits[linunknown]""") + "\n"
            contents += "{{{}}}".format(",".join(x_linunknown)) + "\n"
            contents += "{{{}}}".format(",".join(lin_active_sboxes)) + "\n"
            contents += "{{{}}}".format(",".join(yz_linunknown)) + "\n"
            contents += trim(r"""\orthrosmarkbits[common]""") + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(common_active_sboxes)) + "\n"
            contents += r"{}" + "\n"
            
        # draw EL
        for r in range(self.RL):
            if self.offset + self.RB + self.RU + self.RM + r < 4:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundbits""") + "\n"
            else:
                contents += trim(r"""\spnwordtobitswitch % switch from nibble view to bit view
                                    \orthrosroundnibble""") + "\n"
                
            x_one = [str(i) for i in range(128) if self.lower_trail["x"][self.RM + r][i] == 1]
            x_unknown = [str(i) for i in range(128) if self.lower_trail["x"][self.RM + r][i] == -1]
            if self.offset + self.RB + self.RU + self.RM + r < 4:
                yz_one = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.lower_trail["y"][self.RM + r][i] == 1]
                yz_unknown = [f"{i}/{self.bit_permutation[self.branchtype][i]}" for i in range(128) if self.lower_trail["y"][self.RM + r][i] == -1]
            else:
                yz_one = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.lower_trail["y"][self.RM + r][i] == 1]
                yz_unknown = [f"{i}/{4*self.nibble_permutation[self.branchtype][i//4] + i%4}" for i in range(128) if self.lower_trail["y"][self.RM + r][i] == -1]
            lin_active_sboxes = [str(i) for i in range(32) if sum([int(self.model.getVarByName(f"pl{self.branchtype}_{r}_{i}_{j}").Xn) for j in range(2)]) >= 1]
            contents += trim(r"""\orthrosmarkbits[linone]""") + "\n"
            contents += "{{{}}}".format(",".join(x_one)) + "\n"
            contents += "{{{}}}".format(",".join(lin_active_sboxes)) + "\n"
            contents += "{{{}}}".format(",".join(yz_one)) + "\n"
            contents += trim(r"""\orthrosmarkbits[linunknown]""") + "\n"
            contents += "{{{}}}".format(",".join(x_unknown)) + "\n"
            contents += r"{}" + "\n"
            contents += "{{{}}}".format(",".join(yz_unknown)) + "\n"
        contents += "\n\n" + r"""\begin{comment}""" + "\n"
        contents += self.attack_summary
        contents += r"""\end{comment}""" + "\n"
        contents += r"""\end{tikzpicture}""" + "\n"
        contents += trim(r"""\end{document}""")
        with open(self.output_file_name, "w") as output_file:
            output_file.write(contents)
