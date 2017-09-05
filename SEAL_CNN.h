#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include "seal.h"
#include "bigpoly.h"
#include "AveragePooling.h"
#include "CipherSaveLoad.h"
#include "utility.h"

using namespace seal;
using namespace std;

vector<Ciphertext> cnnCalMnist(vector<vector<Ciphertext>> x, FractionalEncoder encoder, Evaluator evaluator, Decryptor decrptor, string path);
vector<Ciphertext> cnnCalMnist_slim(vector<vector<Ciphertext>> x, FractionalEncoder encoder, Evaluator evaluator, Decryptor decrptor, string path);

vector<vector<vector<vector<Plaintext> > > > convertTensor(vector<vector<double> > &mat, FractionalEncoder encoder, int sizeX, int sizeY, int sizeZ, int nb_nodes);
Ciphertext weightedSum2D(double row, double col, vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator);
vector<vector<Ciphertext> > cnn2D(vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator);
vector<vector<Ciphertext> > compressCNN(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<Plaintext> > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, double nodeBias, bool relu, bool relu_relin);
Ciphertext relu_appr(Ciphertext x, Evaluator evaluator, FractionalEncoder encoder, bool relin);
vector<vector<vector<Ciphertext> > > cnn3D(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<vector<Plaintext> > > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, int nb_nodes, vector<vector<double> > b, bool relu, bool relu_relin);
vector<Ciphertext> flatten3D_cipher(vector<vector<vector<Ciphertext> > > &matrix);
Ciphertext weightedSum1D_ciphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, Plaintext bias, Evaluator evaluator, int col);
vector<Ciphertext> fullyConnectCiphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, vector<vector<Plaintext>> &b, Evaluator evaluator, FractionalEncoder encodor, bool relu);
void printFlat1D(vector<Ciphertext> flat, FractionalEncoder encoder, Decryptor decryptor);
void printCipher2D(vector<vector<Ciphertext> > mat, FractionalEncoder encoder, Decryptor decryptor);