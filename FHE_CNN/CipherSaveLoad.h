#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "seal.h"
#include "bigpoly.h"

using namespace seal;
using namespace std;

// save the 2D ciphertest array into binary files. Each file corresponds to one point.
void ciphertextSave2D(vector<vector<Ciphertext>> &x, string path);
//path: path+ "prename". The file names will be prename00, prename01...
void ciphertextLoad2D(vector<vector<Ciphertext>> &x, string path);
// save keys
void saveKeys(Ciphertext publicKey, Plaintext secretKey, EvaluationKeys evaluationKey, string path);
// load keys
void loadKeys(Ciphertext &publicKey, Plaintext &secretKey, EvaluationKeys &evaluationKey, string path);
// save parms
void saveParms(EncryptionParameters parms, string path);
// load parms
void loadParms(EncryptionParameters &parms, string path);
