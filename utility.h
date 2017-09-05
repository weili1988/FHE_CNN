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

using namespace seal;
using namespace std;

vector<vector<double> > readCSV(string pathFile);
vector<vector<Ciphertext>> encrypt(vector<vector<double> > loadValues, FractionalEncoder encoder, Encryptor encryptor);
void print_example_banner(string title);
vector<vector<Plaintext>> encodeWB(vector<vector<double>> weights, FractionalEncoder encoder);
vector<vector<Ciphertext>> fullyConnect(vector<vector<Ciphertext>> X, vector<vector<Plaintext>> w, vector<vector<Plaintext>> b, FractionalEncoder encoder, Evaluator evaluator);
vector<vector<Ciphertext>> fc_BN_Relu(vector<vector<Ciphertext>> X, vector<vector<Plaintext>> w, vector<vector<Plaintext>> b, FractionalEncoder encoder, Evaluator evaluator);
void testLoading(vector<vector<double> > dataX);
vector<vector<Plaintext>> transpose(vector<vector<Plaintext>> w);
Ciphertext weightedSum(vector<Ciphertext> x, vector<Plaintext> w, vector<Plaintext> b, Evaluator evaluator);
vector<vector<double>> decypt2D(vector<vector<Ciphertext>> encrypted, Decryptor decryptor, FractionalEncoder encoder);
Ciphertext relu(Ciphertext x, FractionalEncoder encoder, Evaluator evaluator, int deg);