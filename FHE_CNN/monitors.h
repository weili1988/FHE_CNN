#pragma once

#include <fstream>
#include <string>
#include "seal.h"

using namespace std;
using namespace seal;

// return plaintxt as string or sring vectors
// overload
string returnPlain(Plaintext pt);
vector<string> returnPlain(vector<Plaintext> pt);
vector<vector<string>> returnPlain(vector<vector<Plaintext>> ct);

// save plaintxt to txt files
// overload
void savePlain(Plaintext ct, string file);
void savePlain(vector<Plaintext> ct, string file);
void savePlain(vector<vector<Plaintext>> ct, string file);

// return ciphertxt as 1D, 2D and 3D string vectors
// overload
vector<string> returnCipher(Ciphertext ct);
vector<vector<string>> returnCipher(vector<Ciphertext> ct);
vector<vector<vector<string>>> returnCipher(vector<vector<Ciphertext>> ct);

// save ciphertxt to txt files
// overload
void saveCipher(Ciphertext ct, string file);
void saveCipher(vector<Ciphertext> ct, string file);
void saveCipher(vector<vector<Ciphertext>> ct, string file);
