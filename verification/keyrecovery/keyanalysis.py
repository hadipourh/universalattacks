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

This program alaizes the valid key space for differential-linear attacks on Orthros-PRF. 
Note that, this script is applicable to key recovery attacks with only 2 active S-boxes in the key recovery part. 
"""

from itertools import product
import time
from math import log2

inv_kp_left = [0, 11, 58, 29, 47, 36, 100, 49, 111, 33, 83, 73, 76, 68, 118, 65, 20, 57, 63, 35, 80, 89, 4, 106, 12, 116, 115, 97, 42, 70, 75, 24, 119, 26, 95, 81, 9, 39, 40, 113, 102, 105, 50, 67, 125, 54, 60, 28, 7, 107, 74, 41, 92, 1, 16, 109, 85, 82, 37, 69, 127, 8, 31, 62, 79, 30, 25, 13, 86, 18, 48, 123, 103, 3, 99, 55, 45, 72, 124, 108, 53, 93, 44, 27, 90, 101, 122, 2, 114, 38, 77, 87, 51, 15, 22, 5, 126, 59, 91, 6, 121, 34, 110, 17, 52, 96, 43, 61, 10, 64, 117, 21, 19, 104, 56, 120, 88, 112, 78, 14, 98, 46, 23, 32, 84, 71, 94, 66]
inv_kp_right = [57, 39, 6, 90, 51, 95, 101, 86, 63, 97, 75, 8, 33, 127, 14, 47, 76, 37, 25, 124, 83, 71, 103, 56, 94, 123, 61, 98, 89, 20, 1, 4, 60, 121, 105, 3, 18, 68, 66, 12, 116, 74, 48, 108, 53, 80, 5, 34, 19, 24, 81, 36, 118, 2, 70, 91, 107, 22, 85, 113, 40, 28, 72, 122, 111, 79, 38, 87, 119, 64, 32, 43, 82, 125, 26, 54, 0, 96, 112, 7, 21, 44, 29, 109, 58, 27, 69, 11, 62, 55, 106, 13, 92, 115, 67, 77, 41, 16, 117, 104, 50, 15, 99, 93, 46, 30, 102, 59, 88, 84, 10, 35, 120, 52, 45, 23, 42, 78, 17, 114, 110, 65, 100, 73, 49, 9, 31, 126]

sb4 = [0x1,0x0,0x2,0x4,0x3,0x8,0x6,0xd,0x9,0xa,0xb,0xe,0xf,0xc,0x7,0x5]
sb8 = [(sb4[a] << 4)&0xf0|sb4[b] for (a, b) in product(range(16), repeat=2)]

n = 8
m = 8

def derive_involved_key_positions(offset, active_indices):
    """
    Given the offset and activeness pattern of S-boxes in 1-round key reccovery, 
    this function return the involved key bits in two branches as well as the intersection of them.
    """

    active_indices = [4*i + j for i in active_indices for j in range(4)]
    inv_kp_left_at_roundr = [[0 for i in range(128)] for _ in range(offset + 1)]
    inv_kp_right_at_roundr = [[0 for i in range(128)] for _ in range(offset + 1)]
    inv_kp_left_at_roundr[0] = [inv_kp_left[i] for i in range(128)]
    inv_kp_right_at_roundr[0] = [inv_kp_right[i] for i in range(128)]

    for r in range(1, offset + 1):
        for i in range(128):
            inv_kp_left_at_roundr[r][i] = inv_kp_left[inv_kp_left_at_roundr[r - 1][i]]
            inv_kp_right_at_roundr[r][i] = inv_kp_right[inv_kp_right_at_roundr[r - 1][i]]

    involved_key_bits_left = []
    involved_key_bits_right = []

    for i in range(128):
        if i in active_indices:
            involved_key_bits_left += [inv_kp_left_at_roundr[offset][i]]
            involved_key_bits_right += [inv_kp_right_at_roundr[offset][i]]

    temp_left = set(involved_key_bits_left)
    temp_right = set(involved_key_bits_right)
    intersection = temp_left.intersection(temp_right)
    common_indices = dict()
    for element in intersection:
        left_index = involved_key_bits_left.index(element)
        right_index = involved_key_bits_right.index(element)
        common_indices[left_index] = right_index
    return involved_key_bits_left, involved_key_bits_right, intersection, common_indices


def check_candidate(k1, k2, dx, dy1, dy2):
    """
    Derive all possible good pairs for a input difference dx and output differences dy1 and dy2.
    """
    # n = 8
    good_pairs = []
    for x in range(2**n):
        if (sb8[x ^ k1] ^ sb8[x ^ k1 ^ dx] == dy1) and (sb8[x ^ k2] ^ sb8[x ^ k2 ^ dx] == dy2):
            x1 = hex(x)[2:].zfill(2)
            x2 = hex(x ^ dx)[2:].zfill(2)
            if not (x1[0], x2[0], x1[1], x2[1]) in good_pairs:
                good_pairs.append((x1[0], x2[0], x1[1], x2[1]))
    return good_pairs

if __name__ == "__main__":
    # Specify the offset and active indices of the S-boxes in the key recovery part.
    # dy1: output difference in the left branch
    # dy2: output difference in the right branch
    #------------------------------------------------------------------------------------
    # offset = 0
    # active_indices = [2, 27]
    # dy1 = 0x44
    # dy2 = 0x81
    #------------------------------------------------------------------------------------
    # 8-round attack version 1:
    # offset = 0
    # active_indices = [27, 29]
    # dy1 = 0x44
    # dy2 = 0x81
    #------------------------------------------------------------------------------------
    # 8-round attack version 5:
    # offset = 0
    # active_indices = [20, 21]
    # dy1 = 0x12
    # dy2 = 0x41
    #------------------------------------------------------------------------------------
    # 8-round attack version 6:
    # offset = 0
    # active_indices = [6, 24]
    # dy1 = 0x42
    # dy2 = 0x41
    #------------------------------------------------------------------------------------
    # 8-round attack version 7:
    # offset = 0
    # active_indices = [20, 21]
    # dy1 = 0x12
    # dy2 = 0x41 # This one is essentially the same as version 5
    #------------------------------------------------------------------------------------
    # 8-round attack version 8:
    # offset = 0
    # active_indices = [15, 18]
    # dy1 = 0x82
    # dy2 = 0x44
    #------------------------------------------------------------------------------------
    # 6-round attack version 2:
    offset = 2
    active_indices = [15, 18]
    dy1 = 0x82
    dy2 = 0x44
    #------------------------------------------------------------------------------------
    # 5-round attack version 1:
    # offset = 2
    # active_indices = [13, 28]
    # dy1 = 0x42
    # dy2 = 0x44
    #------------------------------------------------------------------------------------
    # 5-round attack version 2:
    # offset = 2
    # active_indices = [6, 13]
    # dy1 = 0x88
    # dy2 = 0x82
    #####################################################################################
    left_ik, right_ik, common_ik, common_indices = derive_involved_key_positions(offset=offset, active_indices=active_indices)    
    total_number_of_keys = 2**(len(left_ik) + len(right_ik) - len(common_ik))
    
    # Print configuration before starting
    print("="*85)
    print("6-Round Orthros Key Analysis (Python version)")
    print("="*85)
    print(f"Configuration:")
    print(f"  offset = {offset}")
    print(f"  active_indices = {active_indices}")
    print(f"  dy1 = 0x{dy1:02x}")
    print(f"  dy2 = 0x{dy2:02x}")
    print(f"\nInvolved key bits:")
    print(f"  left involved key bits:  {left_ik}")
    print(f"  right involved key bits: {right_ik}")
    print(f"  Common key bits: {common_ik}")
    print(f"  Total number of keys: {total_number_of_keys}")
    print(f"\nStarting key iteration (this may take a while)...")
    print("="*85)
    
    start_time = time.time()
    cnt = 0
    key_dictionary = {}
    
    for k1, k2 in product(range(2**n), repeat=2):
        good_pairs = []
        # Check the key conditions (dependency between the key bits)
        flag = True
        for item in common_indices.items():
            ls = 8 - item[0] - 1
            rs = 8 - item[1] - 1
            if (k1 >> ls)&0x1 != (k2 >> rs)&0x1:                
                flag = False
                break
        if flag:
            for dx in range(1, 2**n):
                good_pairs.extend(check_candidate(k1, k2, dx, dy1, dy2))
            k1, k2 = (((k1>>4)&0xf)<<0x4|(k2>>4)&0xf), ((k1&0xf)<<0x4|k2&0xf)
            if good_pairs != []:
                good_pairs.sort()
                key_dictionary[tuple(good_pairs)] = key_dictionary.get(tuple(good_pairs), []) + [(k1, k2)]
                cnt += 1
    
    # Print weak-key classes
    print("\nWeak-key classes:")
    print("=" * 85)
    sum_Xi = 0
    for i, (key, value) in enumerate(key_dictionary.items()):
        value = set(value)
        num_pairs = len(key)  # |X_i| is the number of good pairs in the class
        sum_Xi += num_pairs
        print(f"Class {i}: Good pairs: {key}, |X_i| = {num_pairs}, num_keys = {len(value)}")

    elapsed_time = time.time() - start_time
    print("#"*85)
    print("Attack summary:")
    print("Elapsed time: {:.2f} seconds".format(elapsed_time))
    print(f"left involved key bits: {left_ik}")
    print(f"right involved key bits: {right_ik}")
    print(f"Common key bits: {common_ik}")
    print(f"Total number of keys: {total_number_of_keys}")
    print(f"Number of weak keys: {cnt}")
    print(f"Number of unique indices for the hash table (indices are sorted good pairs): {len(key_dictionary)}")
    print(f"Sum of |X_i| across all classes: {sum_Xi} (log2 ≈ {log2(sum_Xi):.2f})")
    print(f"The conditions holds for {cnt} keys out of {total_number_of_keys} keys.")
    log2_weak_key_space_size = log2(cnt/total_number_of_keys)
    learned_bits_in_the_worst_case = abs(log2((total_number_of_keys - cnt)/total_number_of_keys))
    print("The attack works for 2^({:.2f}) portion of the key space.".format(log2_weak_key_space_size))
    print("Number of learned bits in the worst case (when the key is not weak): {:0.2f}".format(learned_bits_in_the_worst_case))