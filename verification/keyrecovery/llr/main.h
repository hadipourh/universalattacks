/*
Copyright (C) 2025 Yosuke Todo and Hosein Hadipour
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

#pragma once
#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

const int modeLeftBr = 0;
const int modeRightBr = 1;
const int modePRF = 2;

class ORTHROS {
public:
   vector<vector<unsigned char>> rkL;
   vector<vector<unsigned char>> rkR;

   void keySchedule(vector<unsigned char> key);
   vector<unsigned char> KEYSCHEDULING(vector<unsigned char> rk, int BRANCH);
   vector<unsigned char> BITPERM(vector<unsigned char> state, int BRANCH);
   vector<unsigned char> NIBBLEPERM(vector<unsigned char> state, int BRANCH);
   vector<unsigned char> Process_Left_Branch(vector<unsigned char> state, int offset, int roundNum);
   vector<unsigned char> Process_Right_Branch(vector<unsigned char> state, int offset, int roundNum);
   vector<unsigned char> encryption(vector<unsigned char> pt, int offset, int roundNum, int mode);
};


void testVector(void);
vector<int> keyIndexUpdateLeft(vector<int> in);
vector<int> keyIndexUpdateRight(vector<int> in);