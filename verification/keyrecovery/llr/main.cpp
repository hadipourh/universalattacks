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

#include "main.h"

#include <string.h>

#include <string>

#include <algorithm>
#include <bitset>
#include <boost/math/distributions/normal.hpp>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <set>
#include <unordered_map>

// This function may have a logical flaw
using namespace std;

namespace {
void logConflictMessage() {
  static int conflictLogCount = 0;
  constexpr int conflictLogLimit = 8;
  if (conflictLogCount < conflictLogLimit) {
    cerr << "Skipping inconsistent weak-key templates: conflicting bit assignments detected." << endl;
    conflictLogCount++;
    if (conflictLogCount == conflictLogLimit) {
      cerr << "Further inconsistent template messages suppressed." << endl;
    }
  }
}
}  // namespace

// the length is always equal
bool mergeString(string str1, string str2, string& out) {
  out = str1;
  for (size_t i = 0; i < str1.size(); ++i) {
    char c1 = str1[i];
    char c2 = str2[i];

    //
    if ((c1 == '0' || c1 == '1') && (c2 == '0' || c2 == '1')) {
      if (c1 != c2) {
        logConflictMessage();
        return false;
      }
    }

    //
    else if ((c1 == '0' || c1 == '1') && c2 == '-') {
      out[i] = c1;
    } else if ((c2 == '0' || c2 == '1') && c1 == '-') {
      out[i] = c2;
    }
  }
  return true;
}

// reverse dictionary
map<vector<int>, vector<int>> reverseDictionary(const map<vector<int>, vector<int>>& input) {
  map<vector<int>, vector<int>> reversed;
  for (const auto& [key, value] : input) {
    if (reversed.find(value) != reversed.end()) {
      throw logic_error("Duplicate values found in the input map, cannot reverse.");
    }
    reversed[value] = key;
  }
  return reversed;
}

//
string convert(int key, vector<int> index) {
  string key128(128, '-');
  for (int i = 0; i < index.size(); i++) {
    if ((key >> i) & 1)
      key128[index[i]] = '1';
    else
      key128[index[i]] = '0';
  }
  return key128;
}
bool checkConvert(int key, vector<int> index, string& out) {
  for (int i = 0; i < index.size(); i++) {
    if ((out[index[i]] == '1') && (((key >> i) & 1) == 0)) {
      return false;
    }
    if ((out[index[i]] == '0') && (((key >> i) & 1) == 1)) {
      return false;
    }

    if ((key >> i) & 1) {
      out[index[i]] = '1';
    } else {
      out[index[i]] = '0';
    }
  }
  return true;
}

// Function to collect good pairs for a given SBOX and delta values
map<vector<string>, vector<int>> collectGoodPairs(const unsigned char SBOX[], int deltaLTarget, int deltaRTarget, vector<int> indexKey, vector<int>& isWeak, int unordered) {
  map<vector<int>, vector<string>> dict;
  for (int kL = 0; kL < 16; kL++) {
    for (int kR = 0; kR < 16; kR++) {
      int valid = true;

      string key(128, '-');
      if (checkConvert((kR << 4) ^ kL, indexKey, key) == false) {
        // invalid key
        continue;
      }

      vector<int> goodPairs;
      for (int x = 0; x < 16; x++) {
        if (unordered == true) {
          for (int y = x; y < 16; y++) {
            int deltaL = SBOX[x ^ kL] ^ SBOX[y ^ kL];
            int deltaR = SBOX[x ^ kR] ^ SBOX[y ^ kR];

            if ((deltaL == deltaLTarget) && (deltaR == deltaRTarget)) {
              int data = (x << 4) ^ y;
              goodPairs.push_back(data);
            }
          }
        } else {
          for (int y = 0; y < 16; y++) {
            int deltaL = SBOX[x ^ kL] ^ SBOX[y ^ kL];
            int deltaR = SBOX[x ^ kR] ^ SBOX[y ^ kR];

            if ((deltaL == deltaLTarget) && (deltaR == deltaRTarget)) {
              int data = (x << 4) ^ y;
              goodPairs.push_back(data);
            }
          }
        }
      }

      dict[goodPairs].push_back(key);

      if (goodPairs.size() > 0)
        isWeak[(kR << 4) ^ kL] = 1;
    }
  }
  map<vector<string>, vector<int>> reverseDict;
  for (const auto& [key, value] : dict) {
    reverseDict[value] = key;
  }

  return reverseDict;
}

bool is_included(vector<string> guessedKey, string key) {
  for (int i = 0; i < guessedKey.size(); i++) {
    if (key == guessedKey[i])
      return true;
  }
  return false;
}
std::vector<int> findDuplicates(std::vector<int> vec) {
  std::sort(vec.begin(), vec.end());

  std::vector<int> result;
  for (size_t i = 1; i < vec.size(); ++i) {
    if (vec[i] == vec[i - 1] && (result.empty() || result.back() != vec[i])) {
      result.push_back(vec[i]);
    }
  }
  return result;
}

// class to store and compute the information related to each distinguisher
class DISTINGUISHER {
 public:
  int offset;
  vector<double> cor;
  vector<int> index;
  vector<int> deltaL;
  vector<int> deltaR;
  vector<int> indexKey0;
  vector<int> indexKey1;
  map<vector<string>, vector<int>> dict0;
  map<vector<string>, vector<int>> dict1;
  vector<int> isWeak0;
  vector<int> isWeak1;

  // to spped up for experiments
  vector<int> numStrongKeysCondition;
  vector<int> numTotalKeysCondition;

  void setUp(int offset, vector<double> cor, vector<int> index, vector<int> deltaL, vector<int> deltaR) {
    // parameter
    this->offset = offset;
    this->cor = cor;
    this->index = index;
    this->deltaL = deltaL;
    this->deltaR = deltaR;

    // key schedule
    vector<int> rkL(128);
    vector<int> rkR(128);
    for (int i = 0; i < 128; i++) {
      rkL[i] = i;
      rkR[i] = i;
    }
    for (int i = 0; i <= offset; i++) {
      rkL = keyIndexUpdateLeft(rkL);
      rkR = keyIndexUpdateRight(rkR);
    }
    for (int j = 0; j < 4; j++) {
      indexKey0.push_back(rkL[4 * index[0] + 3 - j]);
      indexKey1.push_back(rkL[4 * index[1] + 3 - j]);
    }
    for (int j = 0; j < 4; j++) {
      indexKey0.push_back(rkR[4 * index[0] + 3 - j]);
      indexKey1.push_back(rkR[4 * index[1] + 3 - j]);
    }

    for (int i = 0; i < indexKey0.size(); i++) {
      cout << indexKey0[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < indexKey1.size(); i++) {
      cout << indexKey1[i] << " ";
    }
    cout << endl;

    //
    unsigned char SBOX[16] = {0x1, 0x0, 0x2, 0x4, 0x3, 0x8, 0x6, 0xd, 0x9, 0xa, 0xb, 0xe, 0xf, 0xc, 0x7, 0x5};

    // the first active nibble
    isWeak0.resize(256, 0);
    dict0 = collectGoodPairs(SBOX, deltaL[0], deltaR[0], indexKey0, isWeak0, true);

    // the second active nibble
    isWeak1.resize(256, 0);
    dict1 = collectGoodPairs(SBOX, deltaL[1], deltaR[1], indexKey1, isWeak1, false);
  }
  double capacity(void) {
    double C = 0;
    for (int i = 1; i < cor.size(); i++) {
      C += pow(cor[i], 2);
    }
    return C;
  }

  // this function assumes everything are independent
  double getEntropy(vector<int> restricIndex, int restricValue) {
    // parameter
    int numStrongKeys = 0;
    int numWeakKeys = 0;

    // first, we identify the total number of keys satisfying this restriction
    auto it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        //
        int numKeys = 0;
        for (int i = 0; i < key0.size(); i++) {
          for (int j = 0; j < key1.size(); j++) {
            string mergedKey;
            if (mergeString(key0[i], key1[j], mergedKey) == false) {
              continue;
            }

            bool valid = true;
            for (int k = 0; k < restricIndex.size(); k++) {
              char c1 = mergedKey[restricIndex[k]];
              if ((c1 == '0') && (((restricValue >> k) & 1) == 1)) {
                valid = false;
              }
              if ((c1 == '1') && (((restricValue >> k) & 1) == 0)) {
                valid = false;
              }
            }
            if (valid == false)
              continue;

            numKeys++;
          }
        }

        // if there is no consistent weak key in the class, continue
        if (numKeys == 0) {
          it1++;
          continue;
        }

        // if there is no good pairs, it's strong key and continue
        if (((*it0).second.size() == 0) || ((*it1).second.size() == 0)) {
          numStrongKeys += numKeys;
          it1++;
          continue;
        }

        // number of weak keys
        numWeakKeys += numKeys;

        it1++;
      }
      it0++;
    }
    int numTotalKeys = numStrongKeys + numWeakKeys;

    cout << log2(numTotalKeys) << endl;

    // compute the conditional entropy
    double entropy = 0;
    it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        //
        int numKeys = 0;
        for (int i = 0; i < key0.size(); i++) {
          for (int j = 0; j < key1.size(); j++) {
            string mergedKey;
            if (mergeString(key0[i], key1[j], mergedKey) == false) {
              continue;
            }

            bool valid = true;
            for (int k = 0; k < restricIndex.size(); k++) {
              char c1 = mergedKey[restricIndex[k]];
              if ((c1 == '0') && (((restricValue >> k) & 1) == 1)) {
                valid = false;
              }
              if ((c1 == '1') && (((restricValue >> k) & 1) == 0)) {
                valid = false;
              }
            }
            if (valid == false)
              continue;

            numKeys++;
          }
        }

        // if there is no consistent weak key in the class, continue
        if (numKeys == 0) {
          it1++;
          continue;
        }

        // if there is no good pairs, it's strong key and continue
        if (((*it0).second.size() == 0) || ((*it1).second.size() == 0)) {
          it1++;
          continue;
        }

        //
        double prob = (double)(numKeys) / (double)(numTotalKeys);
        entropy -= prob * log2(prob);

        it1++;
      }
      it0++;
    }
    double prob = (double)(numStrongKeys) / (double)(numTotalKeys);
    entropy -= prob * log2(prob);

    // cout << entropy << endl;

    return entropy;
  }
  double getEntropy(void) {
    // pick overlapped index
    vector<int> mergeIndex;
    mergeIndex.insert(mergeIndex.end(), indexKey0.begin(), indexKey0.end());
    mergeIndex.insert(mergeIndex.end(), indexKey1.begin(), indexKey1.end());

    mergeIndex = findDuplicates(mergeIndex);
    printf("overlapped index : ");
    for (int i = 0; i < mergeIndex.size(); i++) {
      cout << mergeIndex[i] << " ";
    }
    cout << endl;

    int numTotalKeys = 1 << (16 - mergeIndex.size());

    // We identify the total number of keys satisfying this restriction
    double entropy = 0;
    int numStrongKeys = 0;
    auto it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        // if the key is not valid, continue
        int numKeys = 0;
        for (int i = 0; i < key0.size(); i++) {
          for (int j = 0; j < key1.size(); j++) {
            string mergedKey;
            if (mergeString(key0[i], key1[j], mergedKey) == false)
              continue;

            numKeys++;
          }
        }

        // if there is no consistent weak key in the class, continue
        if (numKeys == 0) {
          it1++;
          continue;
        }

        // if there is no good pairs, it's strong key and continue
        if (((*it0).second.size() == 0) || ((*it1).second.size() == 0)) {
          numStrongKeys += numKeys;
          it1++;
          continue;
        }

        // average entropy
        double prob = (double)(numKeys) / (double)(numTotalKeys);
        entropy += prob * log2(1.0 / prob);

        it1++;
      }
      it0++;
    }

    // entropy for the strong key
    double prob = (double)(numStrongKeys) / (double)(numTotalKeys);
    entropy += prob * log2(1.0 / prob);

    return entropy;
  }
  void evaluateStatistics(double N, double success_probability) {
    // parameter
    int numStrongKeys = 0;
    double expected_false_positive_number = 0;
    int numWeakKeys = 0;
    double time = 0;
    int weakKeyClassCount = 0;

    //
    auto it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        //
        int numKeys = 0;
        for (int i = 0; i < key0.size(); i++) {
          for (int j = 0; j < key1.size(); j++) {
            string mergedKey;
            if (mergeString(key0[i], key1[j], mergedKey) == true) {
              numKeys++;
            }
          }
        }

        // if there is no consistent weak key in the class, continue
        if (numKeys == 0) {
          it1++;
          continue;
        }

        // if there is no good pairs, it's strong key and continue
        if (((*it0).second.size() == 0) || ((*it1).second.size() == 0)) {
          numStrongKeys += numKeys;
          it1++;
          continue;
        }

        // number of weak keys
        numWeakKeys += numKeys;
        weakKeyClassCount++;

        //
        vector<int> data0 = (*it0).second;
        vector<int> data1 = (*it1).second;
        double sample = N * data0.size() * data1.size();
        time += data0.size() * data1.size();

        // right
        double rightAvg = 0.5 * sample * capacity();
        double rightVar = sample * capacity();
        boost::math::normal distRight(rightAvg, sqrt(rightVar));
        double th = boost::math::quantile(complement(distRight, success_probability));

        // wrong
        double wrongAvg = -0.5 * sample * capacity();
        double wrongVar = sample * capacity();
        boost::math::normal distWrong(wrongAvg, sqrt(wrongVar));
        double false_positive = boost::math::cdf(complement(distWrong, th));
        expected_false_positive_number += false_positive;

        //
        it1++;
      }
      it0++;
    }
    double entropy = getEntropy();

    cout << "# of weak key classes         : " << weakKeyClassCount << endl;
    cout << "# of weak keys (with ks)      : " << numWeakKeys << endl;
    cout << "# of strong keys (with ks)    : " << numStrongKeys << endl;
    cout << "capacity (log 2)              : " << log2(capacity()) << endl;
    cout << "data complexity (log 2)       : " << log2(256 * N) << endl;
    cout << "entropy                       : " << entropy << endl;
    cout << "sum_i |X_i| (log 2)           : " << log2(time) << endl;
    cout << "success probability           : " << success_probability << endl;
    cout << "expected # of false positives : " << expected_false_positive_number << endl;
    cout << endl;
  }

  void precomputeConditionalNumers(vector<int> restricIndex, int restricValue) {
    // parameter
    int numStrongKeys = 0;
    int numWeakKeys = 0;

    // first, we identify the total number of keys satisfying this restriction
    auto it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        //
        int numKeys = 0;
        for (int i = 0; i < key0.size(); i++) {
          for (int j = 0; j < key1.size(); j++) {
            string mergedKey;
            mergeString(key0[i], key1[j], mergedKey);

            bool valid = true;
            for (int k = 0; k < restricIndex.size(); k++) {
              char c1 = mergedKey[restricIndex[k]];
              if ((c1 == '0') && (((restricValue >> k) & 1) == 1)) {
                valid = false;
              }
              if ((c1 == '1') && (((restricValue >> k) & 1) == 0)) {
                valid = false;
              }
            }
            if (valid == false)
              continue;

            numKeys++;
          }
        }

        // if there is no good pairs, it's strong key and continue
        if (((*it0).second.size() == 0) || ((*it1).second.size() == 0)) {
          numStrongKeys += numKeys;
        } else {
          // number of weak keys
          numWeakKeys += numKeys;
        }
        it1++;
      }
      it0++;
    }
    int numTotalKeys = numStrongKeys + numWeakKeys;

    numStrongKeysCondition.push_back(numStrongKeys);
    numTotalKeysCondition.push_back(numTotalKeys);

    cout << restricValue << "\t" << numStrongKeys << "\t" << numTotalKeys << endl;
  }
  long double getConditionalSelfInformationExperimentaly(vector<int> restricIndex, int restricValue, vector<unsigned char> key, bool debug = false) {
    // key schedule
    ORTHROS cipher;
    cipher.keySchedule(key);

    // key
    int rk0 = (cipher.rkR[offset][index[0]] << 4) ^ cipher.rkL[offset][index[0]];
    int rk1 = (cipher.rkR[offset][index[1]] << 4) ^ cipher.rkL[offset][index[1]];

    // string mergedKey;
    // mergeString(convert(rk0, indexKey0), convert(rk1, indexKey1), mergedKey);
    // cout << mergedKey << endl;

    // compute the conditional entropy with given key
    auto it0 = dict0.begin();
    while (it0 != dict0.end()) {
      auto it1 = dict1.begin();
      while (it1 != dict1.end()) {
        // check one weak-key class
        vector<string> key0 = (*it0).first;
        vector<string> key1 = (*it1).first;

        // this key
        if ((is_included(key0, convert(rk0, indexKey0)) && is_included(key1, convert(rk1, indexKey1)))) {
          // if strong key, return the probability
          vector<int> data0 = (*it0).second;
          vector<int> data1 = (*it1).second;
          if (data0.size() * data1.size() == 0) {
            long double p = (long double)(numStrongKeysCondition[restricValue]) / (long double)(numTotalKeysCondition[restricValue]);
            return p;
          }

          // count the number of real keys in the same weak-key class
          int numKeys = 0;
          for (int i = 0; i < key0.size(); i++) {
            for (int j = 0; j < key1.size(); j++) {
              string mergedKey;
              if (mergeString(key0[i], key1[j], mergedKey) == false)
                continue;

              bool valid = true;
              for (int k = 0; k < restricIndex.size(); k++) {
                char c1 = mergedKey[restricIndex[k]];
                if ((c1 == '0') && (((restricValue >> k) & 1) == 1)) {
                  valid = false;
                }
                if ((c1 == '1') && (((restricValue >> k) & 1) == 0)) {
                  valid = false;
                }
              }
              if (valid == false)
                continue;

              numKeys++;
            }
          }

          // no valid key, it's probabiliy 0
          if (numKeys == 0) {
            return 0;
          }

          // weak key
          long double p = (long double)(numKeys) / (long double)(numTotalKeysCondition[restricValue]);
          return p;
        }

        it1++;
      }
      it0++;
    }

    cerr << "error" << endl;
    return 0;
  }

  // check distinguisher experimentally
  void experiments(int N, double success_probability) {
    // parameter
    int round = 6;
    int m = 4;
    int M = 1 << m;

    // correlation
    vector<double> bias(M, 0);
    for (int eta = 0; eta < M; eta++) {
      for (int a = 1; a < M; a++) {
        bias[eta] += pow(2, -m) * pow(-1, __builtin_popcountl(a & eta) % 2) * cor[a];
      }
    }

    //
    int num_find_solution = 0;
    int num_weak_key_trial = 0;
    int num_strong_key_trial = 0;
    vector<double> sum_correlation(16, 0);
    for (int loop = 0; loop < 1000; loop++) {
      std::random_device seed_gen;
      std::mt19937 mtk(seed_gen());

      // set random key
      vector<unsigned char> key(32, 0);
      for (int x = 0; x < 32; x++) {
        key[x] = mtk() & 0xF;
      }

      // set cipher
      ORTHROS cipher;
      cipher.keySchedule(key);

      // key
      int rk0 = (cipher.rkR[offset][index[0]] << 4) ^ cipher.rkL[offset][index[0]];
      int rk1 = (cipher.rkR[offset][index[1]] << 4) ^ cipher.rkL[offset][index[1]];

      // key
      int is_weak = isWeak0[rk0] * isWeak1[rk1];
      if (is_weak == 0) {
        loop--;
        continue;

        printf("%4d ", loop + 1);
        printf("key0 = %02X ", rk0);
        printf("key1 = %02X ", rk1);
        cout << "\t";
        cout << "strong  \t";
        num_strong_key_trial++;
        cout << endl;
      } else {
        // loop--;
        // continue;

        printf("%4d ", loop + 1);
        printf("key0 = %02X ", rk0);
        printf("key1 = %02X ", rk1);
        cout << "\t";
        cout << "weak  \t";
        num_weak_key_trial++;
        cout << endl;
      }
      cout << endl;

      // query and store
      auto start = std::chrono::high_resolution_clock::now();

      vector<vector<unsigned char>> table(256, vector<unsigned char>(N));
#pragma omp parallel
      {
        std::random_device seed_gen;
        std::mt19937 mt(seed_gen());

#pragma omp for
        for (int i = 0; i < N; i++) {
          vector<unsigned char> plaintext(32, 0);
          for (int x = 0; x < 32; x++)
            plaintext[x] = mt() & 0xF;

          for (int x15 = 0; x15 < 16; x15++) {
            for (int x18 = 0; x18 < 16; x18++) {
              plaintext[15] = x15;
              plaintext[18] = x18;

              vector<unsigned char> ct = cipher.encryption(plaintext, offset, round, modePRF);
              table[(x15 << 4) ^ x18][i] = ct[29];
            }
          }
        }
      }

      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      fprintf(stdout, "data collection  %lf[s]\t", duration / 1e9);
      cout << endl;

      // check statistics
      int num_detection = 0;
      auto it0 = dict0.begin();
      while (it0 != dict0.end()) {
        auto it1 = dict1.begin();
        while (it1 != dict1.end()) {
          // candidate of good pairs
          vector<int> data0 = (*it0).second;
          vector<int> data1 = (*it1).second;

          if (data0.size() * data1.size() == 0) {
            it1++;
            continue;
          }

          // guessed keys
          vector<string> guessedKey0 = (*it0).first;
          vector<string> guessedKey1 = (*it1).first;

          /*
          // for debug : only check the correct key
          if ((is_included(guessedKey0, convert(rk0, indexKey0)) & is_included(guessedKey1, convert(rk1, indexKey1)))) {
            for(int i=0;i<guessedKey0.size();i++){
              for(int j=0;j<guessedKey1.size();j++){
                cout << guessedKey0[i] << " " << guessedKey1[j] << endl;
               }
            }
          } else {
            it1++;
            continue;
          }
          */

          //
          double mu0 = 0;
          for (int i = 0; i < M; i++) {
            mu0 += (pow(2, -m) + bias[i]) * bias[i] / pow(2, -m);
          }
          int sample = N * data0.size() * data1.size();
          double rightAvg = sample * mu0;
          double rightVar = sample * capacity();

          boost::math::normal distRight(rightAvg, sqrt(rightVar));
          double th = boost::math::quantile(complement(distRight, success_probability));

          //
          vector<long long int> array(256, 0);
          for (int x0 = 0; x0 < data0.size(); x0++) {
            for (int x1 = 0; x1 < data1.size(); x1++) {
              int a = (data0[x0] >> 4) & 0xF;
              int b = (data0[x0] >> 0) & 0xF;

              int c = (data1[x1] >> 4) & 0xF;
              int d = (data1[x1] >> 0) & 0xF;

#pragma omp parallel
              {
                vector<long long int> arrayTmp(256, 0);

#pragma omp for
                for (int i = 0; i < N; i++) {
                  unsigned char x = table[(a << 4) ^ c][i] ^ table[(b << 4) ^ d][i];
                  arrayTmp[x]++;
                }

#pragma omp critical
                {
                  for (int i = 0; i < 256; i++) {
                    array[i] += arrayTmp[i];
                  }
                }
              }
            }
          }

          // compute LLR statistics
          double LLR = 0;
          for (int i = 0; i < 16; i++) {
            LLR += array[i] * bias[i] * M;
          }

          if (LLR > th) {
            num_detection++;

            if ((is_included(guessedKey0, convert(rk0, indexKey0)) & is_included(guessedKey1, convert(rk1, indexKey1)))) {
              num_find_solution++;
              cout << "\t\tfind" << "\t";
            } else {
              cout << "\t\t    " << "\t";
            }

            printf("{15 : ");
            for (int i = 0; i < data0.size(); i++) {
              printf("%02X ", data0[i]);
            }
            printf("}");

            printf("{18 : ");
            for (int i = 0; i < data1.size(); i++) {
              printf("%02X ", data1[i]);
            }
            printf("}");
            cout << endl;
          }

          it1++;
        }
        it0++;
      }
      double ex_success_prob = (double)num_find_solution / (double)num_weak_key_trial;

      cout << "\t\t # detection            : " << num_detection << endl;
      cout << "\t\t observed success prob. : " << ex_success_prob << endl;
      cout << endl;
    }
  }
};

void getSummary(vector<long double> listSi){
  int size = listSi.size();
  cout << size << "\t";

  sort(listSi.begin(), listSi.end());
  cout << "worst : " << listSi[0] << "\t";
  cout << "1/4 : " << listSi[(size - 1) / 4] << "\t";
  cout << "1/2 : " << listSi[(size - 1) / 2] << "\t";
  cout << "3/4 : " << listSi[3 * (size - 1) / 4] << "\t";
  cout << "best : " << listSi[size - 1] << "\t";

  long double sum = 0;
  for(int i=0;i<size;i++){
    sum += listSi[i];
  }
  sum /= size;
  cout << "average : " << sum << endl;

}
double getMergeEntropy(vector<DISTINGUISHER> dises) {
  // pick overlapped index
  vector<int> mergeIndex;
  for (int i = 0; i < dises.size(); i++) {
    vector<int> tmp = dises[i].indexKey0;
    tmp.insert(tmp.end(), dises[i].indexKey1.begin(), dises[i].indexKey1.end());
    sort(tmp.begin(), tmp.end());
    tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());

    mergeIndex.insert(mergeIndex.end(), tmp.begin(), tmp.end());
  }

  mergeIndex = findDuplicates(mergeIndex);
  printf("overlapped index : ");
  for (int i = 0; i < mergeIndex.size(); i++) {
    cout << mergeIndex[i] << " ";
  }
  cout << endl;

  int outerLoopSize = 1 << mergeIndex.size();
  long double mergeEntropy = 0;
  for (int i = 0; i < outerLoopSize; i++) {
    for (int j = 0; j < dises.size(); j++) {
      mergeEntropy += dises[j].getEntropy(mergeIndex, i);
    }
  }
  mergeEntropy /= outerLoopSize;

  cout << "merge entropy : " << mergeEntropy << endl;

  return mergeEntropy;
}
void checkMergeEntropyExperimentally(vector<DISTINGUISHER> dises, const string& exportPath = "") {
  // random
  std::random_device seed_gen;
  std::mt19937_64 mt(seed_gen());

  // pick overlapped index
  vector<int> mergeIndex;
  for (int i = 0; i < dises.size(); i++) {
    vector<int> tmp = dises[i].indexKey0;
    tmp.insert(tmp.end(), dises[i].indexKey1.begin(), dises[i].indexKey1.end());
    sort(tmp.begin(), tmp.end());
    tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());

    mergeIndex.insert(mergeIndex.end(), tmp.begin(), tmp.end());
  }
  mergeIndex = findDuplicates(mergeIndex);
  printf("overlapped index : ");
  for (int i = 0; i < mergeIndex.size(); i++) {
    cout << mergeIndex[i] << " ";
  }
  cout << endl;
  // mergeIndex.clear();

  // precompute each size
  int outerLoopSize = 1 << mergeIndex.size();
  for (int i = 0; i < outerLoopSize; i++) {
    for (int j = 0; j < dises.size(); j++) {
      dises[j].precomputeConditionalNumers(mergeIndex, i);
    }
  }

  // get self information
  vector<long double> listSi;
  ofstream histOut;
  if (!exportPath.empty()) {
    histOut.open(exportPath, ios::out | ios::trunc);
  }
  for (int loop = 0; loop < 10000; loop++) {
    // randomly choose the key
    vector<unsigned char> key(32);
    for (int i = 0; i < 32; i++) {
      key[i] = mt() & 0xF;
    }

    ORTHROS cipher;
    cipher.keySchedule(key);


    int outerLoopSize = 1 << mergeIndex.size();
    long double mergeProbability = 0;
    for (int i = 0; i < outerLoopSize; i++) {
      long double conditionalProbability = 1;
      for (int j = 0; j < dises.size(); j++) {
        conditionalProbability *= dises[j].getConditionalSelfInformationExperimentaly(mergeIndex, i, key);
      }
      mergeProbability += conditionalProbability;
    }
    mergeProbability /= outerLoopSize;

    // cout << mergeProbability << endl;
    double mergeSi = log2(1 / mergeProbability);
    listSi.push_back(mergeSi);

    if (histOut.is_open()) {
      histOut << fixed << setprecision(8) << mergeSi << '\n';
    }

    if ((loop % 100) == 99){
      getSummary(listSi);
    }
  }

  if (histOut.is_open()) {
    histOut.close();
  }

}

// evaluate the attack
void test1(const string& histogramPath) {
  // parameter
  int m = 4;
  double N = pow(2.0, 32.0);
  double success_probability = 0.99;

  // distinguisher 0
  cout << "distinguisher 0" << endl;
  DISTINGUISHER dis0;
  vector<double> cor0 = {1, pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11)};
  dis0.setUp(2, cor0, {15, 18}, {8, 2}, {4, 4});
  dis0.evaluateStatistics(N, success_probability);

  // distinguisher 1
  cout << "distinguisher 1" << endl;
  DISTINGUISHER dis1;
  vector<double> cor1 = {1, pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15), pow(2, -15)};
  dis1.setUp(2, cor1, {20, 21}, {1, 2}, {4, 1});
  dis1.evaluateStatistics(N, success_probability);

  // distinguisher 2
  cout << "distinguisher 2" << endl;
  DISTINGUISHER dis2;
  vector<double> cor2 = {
      1, pow(2, -9.04 - 5), pow(2, -8.86 - 5), pow(2, -9.16 - 5), pow(2, -8.84 - 5), pow(2, -8.84 - 5), pow(2, -8.78 - 5), pow(2, -8.93 - 5), pow(2, -8.61 - 5), pow(2, -9.13 - 5), pow(2, -9.13 - 5), pow(2, -9.33 - 5), pow(2, -8.63 - 5), pow(2, -9.34 - 5), pow(2, -8.46 - 5), pow(2, -9.20 - 5)};
  dis2.setUp(2, cor2, {23, 31}, {4, 4}, {4, 4});
  dis2.evaluateStatistics(N, success_probability);

  // distinguisher 3
  cout << "distinguisher 3" << endl;
  DISTINGUISHER dis3;
  vector<double> cor3 = {1, pow(2, -15.74), pow(2, -16.09), pow(2, -15.95), pow(2, -15.88), pow(2, -15.92), pow(2, -15.91), pow(2, -15.79), pow(2, -16.08), pow(2, -15.84), pow(2, -15.84), pow(2, -15.74), pow(2, -15.48), pow(2, -16.22), pow(2, -15.07), pow(2, -16.33)};
  dis3.setUp(2, cor3, {5, 10}, {4, 4}, {4, 4});
  dis3.evaluateStatistics(N, success_probability);

  //cout << "dis0, dis2" << endl;
  //getMergeEntropy({dis0, dis2});
  //checkMergeEntropyExperimentally({dis0, dis2});

  cout << "dis0, dis1, dis2, dis3" << endl;
  getMergeEntropy({dis0, dis1, dis2, dis3});
  checkMergeEntropyExperimentally({dis0, dis1, dis2, dis3}, histogramPath);
  return;
}

// attack experiment
void test0(void) {
  // parameter
  int m = 4;
  int N = 1 << 20;
  double success_probability = 0.7;

  // distinguisher 0
  cout << "distinguisher 0" << endl;
  DISTINGUISHER dis0;
  vector<double> cor0 = {1, pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11), pow(2, -11)};
  dis0.setUp(2, cor0, {15, 18}, {8, 2}, {4, 4});
  dis0.evaluateStatistics((double)N, success_probability);
  dis0.experiments(N, success_probability);
}

void test7(const string& histogramPath) {
  double N = pow(2.0, 104);
  double success_probability = 0.99;

  cout << "distinguisher 0" << endl;
  DISTINGUISHER dis0;
  vector<double> cor0 = {
      1.0,
      pow(2.0, -44.02),
      pow(2.0, -43.30),
      pow(2.0, -43.20),
      pow(2.0, -44.19),
      pow(2.0, -43.68),
      pow(2.0, -43.81),
      pow(2.0, -43.80),
      pow(2.0, -43.42),
      pow(2.0, -44.15),
      pow(2.0, -43.70),
      pow(2.0, -43.05),
      pow(2.0, -43.28),
      pow(2.0, -43.77),
      pow(2.0, -43.22),
      pow(2.0, -43.89)};
  dis0.setUp(0, cor0, {27, 29}, {0x4, 0x4}, {0x8, 0x1});
  dis0.evaluateStatistics(N, success_probability);

  cout << "distinguisher 1" << endl;
  DISTINGUISHER dis1;
  vector<double> cor1 = {
      1.0,
      -pow(2.0, -48.02),
      -pow(2.0, -47.93),
      -pow(2.0, -48.05),
      -pow(2.0, -47.81),
      -pow(2.0, -47.83),
      -pow(2.0, -47.75),
      -pow(2.0, -47.85),
      -pow(2.0, -47.57),
      -pow(2.0, -48.12),
      -pow(2.0, -48.17),
      -pow(2.0, -48.25),
      -pow(2.0, -47.66),
      -pow(2.0, -48.25),
      -pow(2.0, -47.45),
      -pow(2.0, -48.22)};
  dis1.setUp(0, cor1, {20, 21}, {0x1, 0x2}, {0x4, 0x1});
  dis1.evaluateStatistics(N, success_probability);

  cout << "distinguisher 2" << endl;
  DISTINGUISHER dis2;
  vector<double> cor2 = {
      1.0,
      0.0,
      pow(2.0, -49.58),
      pow(2.0, -50.13),
      pow(2.0, -50.35),
      pow(2.0, -49.46),
      pow(2.0, -50.25),
      pow(2.0, -50.42),
      pow(2.0, -48.75),
      pow(2.0, -51.73),
      pow(2.0, -50.67),
      pow(2.0, -49.93),
      pow(2.0, -50.04),
      pow(2.0, -48.83),
      0.0,
      pow(2.0, -48.71)};
  dis2.setUp(0, cor2, {6, 24}, {0x4, 0x2}, {0x4, 0x1});
  dis2.evaluateStatistics(N, success_probability);

  cout << "distinguisher 3" << endl;
  DISTINGUISHER dis3;
  vector<double> cor3 = {
      1.0,
      -pow(2.0, -52.14),
      -pow(2.0, -51.62),
      -pow(2.0, -51.81),
      -pow(2.0, -52.29),
      -pow(2.0, -51.92),
      -pow(2.0, -52.00),
      -pow(2.0, -52.05),
      -pow(2.0, -51.40),
      -pow(2.0, -52.51),
      -pow(2.0, -52.21),
      -pow(2.0, -51.52),
      -pow(2.0, -51.50),
      -pow(2.0, -51.53),
      -pow(2.0, -51.43),
      -pow(2.0, -51.73)};
  dis3.setUp(0, cor3, {15, 18}, {0x8, 0x2}, {0x4, 0x4});
  dis3.evaluateStatistics(N, success_probability);

  cout << "dis0, dis1, dis2, dis3" << endl;
  getMergeEntropy({dis0, dis1, dis2, dis3});
  checkMergeEntropyExperimentally({dis0, dis1, dis2, dis3}, histogramPath);
}

namespace {
void printUsage(const char* programName) {
  cout << "Usage: " << programName << " <command> [options]\n\n";
  cout << "Commands:\n";
  cout << "  verify-6r          Reproduce the 6-round verification experiment (test0).\n";
  cout << "  histogram-5r       Sample bit information for the four 5-round distinguishers (test1).\n";
  cout << "  histogram-7r       Sample bit information for the four 7-round distinguishers (test7).\n";
  cout << "  help               Show this message.\n\n";
  cout << "Histogram options (for histogram-5r and histogram-7r):\n";
  cout << "  --histogram <file>  Write histogram samples to <file>.\n";
  cout << "  --histogram=<file>  Same as above.\n";
  cout << "  --no-histogram      Skip writing histogram samples to disk.\n\n";
}

bool parseHistogramOptions(int argc, char** argv, int startIndex, const string& defaultPath, string& histogramPath) {
  histogramPath = defaultPath;
  const string histogramPrefix = "--histogram=";
  for (int i = startIndex; i < argc; ++i) {
    string arg = argv[i];
    if ((arg == "--histogram" || arg == "-o") && (i + 1 < argc)) {
      histogramPath = argv[++i];
    } else if (arg.rfind(histogramPrefix, 0) == 0) {
      histogramPath = arg.substr(histogramPrefix.size());
    } else if (arg == "--no-histogram") {
      histogramPath.clear();
    } else {
      cerr << "Unknown option: " << arg << "\n";
      return false;
    }
  }
  return true;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    printUsage(argv[0]);
    return 1;
  }

  string command = argv[1];
  if (command == "help" || command == "--help" || command == "-h") {
    printUsage(argv[0]);
    return 0;
  }

  if (command == "verify-6r") {
    cout << "Running 6-round verification experiment (test0)" << endl;
    test0();
    return 0;
  }

  if (command == "histogram-5r") {
    string histogramPath;
    if (!parseHistogramOptions(argc, argv, 2, "histogram-5.csv", histogramPath)) {
      printUsage(argv[0]);
      return 1;
    }

    cout << "Running 5-round distinguisher histogram experiment (test1)" << endl;
    test1(histogramPath);
    return 0;
  }

  if (command == "histogram-7r") {
    string histogramPath;
    if (!parseHistogramOptions(argc, argv, 2, "histogram-7.csv", histogramPath)) {
      printUsage(argv[0]);
      return 1;
    }

    cout << "Running 7-round distinguisher histogram experiment (test7)" << endl;
    test7(histogramPath);
    return 0;
  }

  if (command == "attack") {
    cerr << "Command 'attack' has been renamed to 'verify-6r'." << endl;
  } else if (command == "merge" || command == "merge6") {
    cerr << "Command '" << command << "' has been renamed to 'histogram-5r'." << endl;
  } else if (command == "merge8") {
    cerr << "Command 'merge8' has been renamed to 'histogram-7r'." << endl;
  }

  cerr << "Unknown command: " << command << "\n";
  printUsage(argv[0]);
  return 1;
}