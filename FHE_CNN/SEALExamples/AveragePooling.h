#pragma once
#include <string>
#include "seal.h"

using namespace std;
using namespace seal;

Ciphertext avgSum(vector<vector<Ciphertext> > &x, int poolSize, int row, int col, FractionalEncoder encoder, Evaluator evaluator);
vector<vector<Ciphertext> > averagePooling2D(vector<vector<Ciphertext> > &x, int poolSize, int strideSize, FractionalEncoder encoder, Evaluator evaluator);
vector<vector<vector<Ciphertext> > > averagePooling3D(vector<vector<vector<Ciphertext> > > &x, int poolSize, int strideSize, FractionalEncoder encoder, Evaluator evaluator);
