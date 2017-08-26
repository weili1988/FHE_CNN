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

using namespace seal;
using namespace std;

vector<vector<double> > readCSV(string pathFile);
vector<vector<Ciphertext>> encrypt(vector<vector<double> > loadValues, FractionalEncoder encoder, Encryptor encryptor);
vector<vector<Plaintext>> encodeWB(vector<vector<double>> weights, FractionalEncoder encoder);
vector<Ciphertext> cnnCalMnist(vector<vector<Ciphertext>> x, FractionalEncoder encoder, Evaluator evaluator, Decryptor decrptor);
void CNN_test();
vector<vector<vector<vector<Plaintext> > > > convertTensor(vector<vector<double> > &mat, FractionalEncoder encoder, int sizeX, int sizeY, int sizeZ, int nb_nodes);
Ciphertext weightedSum2D(double row, double col, vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator);
vector<vector<Ciphertext> > cnn2D(vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator);
vector<vector<Ciphertext> > compressCNN(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<Plaintext> > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, double nodeBias, bool relu);
Ciphertext relu_appr(Ciphertext x, Evaluator evaluator, FractionalEncoder encoder);
vector<vector<vector<Ciphertext> > > cnn3D(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<vector<Plaintext> > > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, int nb_nodes, vector<vector<double> > b, bool relu);
vector<Ciphertext> flatten3D_cipher(vector<vector<vector<Ciphertext> > > &matrix);
Ciphertext weightedSum1D_ciphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, Plaintext bias, Evaluator evaluator, int col);
vector<Ciphertext> fullyConnectCiphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, vector<vector<Plaintext>> &b, Evaluator evaluator, FractionalEncoder encodor, bool relu);

void CNN_test(){
	string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights";
	string x0, w1, b1, w2, b2, fcW, fc_b;

	x0.append(path);
	x0.append("/7.csv");

	vector<vector<double> > input = readCSV(x0); // read MNIST number
												 // Create encryption parameters
	EncryptionParameters parms;
	parms.set_poly_modulus("1x^4096 + 1");
	parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
	parms.set_plain_modulus(1 << 6);
	//parms.set_decomposition_bit_count(16); // this number needs to be small to decrease noise;
	parms.set_decomposition_bit_count(32); // this number needs to be small to decrease noise;
	parms.validate(); // parms.validate() is the end of parms creation

	// Generate keys.
	cout << "Generating keys ..." << endl;
	KeyGenerator generator(parms);
	generator.generate(4); // int parameter is the evaluation key count 
	cout << "... key generation complete" << endl;
	Ciphertext public_key = generator.public_key();
	Plaintext secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();
	/*
	We will need a fractional encoder for dealing with the rational numbers. Here we reserve
	64 coefficients of the polynomial for the integral part (low-degree terms) and expand the
	fractional part to 32 terms of precision (base 3) (high-degree terms).
	*/
	FractionalEncoder encoder(parms.plain_modulus(), parms.poly_modulus(), 64, 32, 3);
	// Create encryptor, evaluator, decryptor
	Encryptor encryptor(parms, public_key);
	Evaluator evaluator(parms, evaluation_keys);
	Decryptor decryptor(parms, secret_key);

	cout << "encrypting x: " << endl;
	vector<vector<Ciphertext>> x = encrypt(input, encoder, encryptor); // encryption

	cout << "CNN calculation: " << endl;
	vector<Ciphertext> cnnEncrptRes = cnnCalMnist(x, encoder, evaluator, decryptor); // calculation // !!!!!! remember to take out the decryptor after development is done!
	
	// Print size and noise budget of result. 
	cout << "Size of weightedSum without relinearization: " << cnnEncrptRes[0].size() << endl;
	cout << "Noise budget in weightedSum without relinearization: " << decryptor.invariant_noise_budget(cnnEncrptRes[0]) << " bits" << endl;
	// decrypt
	vector<Plaintext> decryptedCnnRes;
	for (int i = 0; i < cnnEncrptRes.size(); i++) {
		decryptedCnnRes.emplace_back(decryptor.decrypt(cnnEncrptRes[i]));
	}
	// decode
	cout << "---------------- Final results: -------------------------" << endl;
	for (int i = 0; i < decryptedCnnRes.size(); i++) {
		cout<< encoder.decode(decryptedCnnRes[i]) << endl;
	}

}
vector<Ciphertext> cnnCalMnist(vector<vector<Ciphertext>> x, FractionalEncoder encoder, Evaluator evaluator, Decryptor decryptor) { // remember to delete this decrypter after development done
	// This is server
	// has encrypted x, encoder, evaluator, (and weights in server); output the vecotor (10 ciphter text) back to user who can decrypt the ciphertext back to plaintext
	string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights";
	string w1, b1, w2, b2, fcW, fc_b;
	w1.append(path);
	w1.append("/CNN1_W.txt");
	b1.append(path);
	b1.append("/CNN1_b.txt");
	w2.append(path);
	w2.append("/CNN2_W.txt");
	b2.append(path);
	b2.append("/CNN2_b.txt");
	fcW.append(path);
	fcW.append("/FC_W.txt");
	fc_b.append(path);
	fc_b.append("/FC_b.txt");

	vector<vector<double> > cnn1W = readCSV(w1);
	vector<vector<double> > cnn1b = readCSV(b1);
	vector<vector<double> > cnn2W = readCSV(w2);
	vector<vector<double> > cnn2b = readCSV(b2);
	vector<vector<double> > fcWeight = readCSV(fcW);
	vector<vector<double> > fc_bias = readCSV(fc_b);
	
	vector<vector<vector<vector<Plaintext> > > > converted_cnn1W = convertTensor(cnn1W, encoder, 3, 3, 1, 8); // convert to tensor with shape (8, 1, 3, 3) and finished encoding
	
	int nb_nodesCNN1W = 8;
	int depth_cnn1W = 1;
	vector<vector<vector<Ciphertext> > > input;
	input.emplace_back(x); // need to make x into tensor with x, y, z dimension
	vector<vector<vector<Ciphertext> > > cnn1_product = cnn3D(input, converted_cnn1W, evaluator, encoder, depth_cnn1W, nb_nodesCNN1W, cnn1b, true);

	int poolsize1 = 2;
	int stridesize1 = 2;
	vector<vector<vector<Ciphertext> > > cnn1_pooled = averagePooling3D(cnn1_product, poolsize1, stridesize1, encoder, evaluator);

	int nb_nodesCNN2W = 16;
	int depth_cnn2W = 8;
	vector<vector<vector<vector<Plaintext> > > > converted_cnn2W = convertTensor(cnn2W, encoder, 3, 3, 8, 16);
	vector<vector<vector<Ciphertext> > > cnn2_product = cnn3D(cnn1_pooled, converted_cnn2W, evaluator, encoder, depth_cnn2W, nb_nodesCNN2W, cnn2b, true);

	int poolsize2 = 4;
	int stridesize2 = 4;
	vector<vector<vector<Ciphertext> > > cnn2_pooled = averagePooling3D(cnn2_product, poolsize2, stridesize2, encoder, evaluator);

	//This is to flat the cnn result
	vector<Ciphertext> flat = flatten3D_cipher(cnn2_pooled);
	//fully connect layer
	vector<vector<Plaintext>> encodedFcW = encodeWB(fcWeight, encoder);
	vector<vector<Plaintext>> encodedFcBias = encodeWB(fc_bias, encoder);
	vector<Ciphertext> fcResult = fullyConnectCiphertext(flat, encodedFcW, encodedFcBias, evaluator, encoder, false); // no relu at last layer
	return fcResult;
	/*------------------------------ below is for testing only-----------------------------------------------------------*/
	/*
	//simple test for data loading and weightedSum2D //
	vector<vector<Ciphertext> > x0;
	for (int i = 0; i < 3; i++) {
		vector<Ciphertext> x0Row;
		for (int j = 0; j < 3; j++) {
			x0Row.emplace_back(x[i][j]);
		}
		x0.emplace_back(x0Row);
	}
	vector<vector<Plaintext> > cnnW = converted_cnn1W[0][0];
	Ciphertext sum = weightedSum2D(0, 0, x0, cnnW, evaluator);

	Plaintext sum_plain = decryptor.decrypt(sum);
	double sum_double =	encoder.decode(sum_plain);

	cout << "weightedSum: " << sum_double << endl;
	*/
	/*
	// test cnn2D
	vector<vector<Plaintext> > cnnW = converted_cnn1W[0][0];
	vector<vector<Ciphertext> > cnn2D_w100 = cnn2D(x, cnnW, evaluator);
	for (int i = 0; i < cnn2D_w100.size(); i++) {
	for (int j = 0; j < cnn2D_w100[0].size(); j++) {
	Plaintext tmpElement = decryptor.decrypt(cnn2D_w100[i][j]);
	cout << encoder.decode(tmpElement) << ((j < cnn2D_w100.at(i).size()) ? ", " : ".\n");
	}
	}
	*/
	/*
	// test cnn3D for 1 node
	vector<vector<vector<Plaintext> > > cnnW = converted_cnn1W[0]; // use the first CNN weigth 1X3X3 (depth = 1)
	int depth = 1;
	double nodeBias = cnn1b[0][0];

	vector<vector<vector<Ciphertext> > > input;
	input.emplace_back(x); // need to make x into tensor with x, y, z dimension

	cout << input.size() << endl;
	cout << input[0].size() << endl;
	cout << input[0][0].size() << endl;

	vector<vector<Ciphertext> > cnn2D_w10 = compressCNN(input, cnnW, evaluator, encoder, depth, nodeBias, false); // no relu
	for (int i = 0; i < cnn2D_w10.size(); i++) {
		for (int j = 0; j < cnn2D_w10[0].size(); j++) {
			Plaintext tmpElement = decryptor.decrypt(cnn2D_w10[i][j]);
			cout << encoder.decode(tmpElement) << ((j < cnn2D_w10.at(i).size()) ? ", " : ".\n");
		}
		cout<< " row " << i << " finished" << endl;
	}
	for (int dp = 0; dp < depth; dp++) {
		for (int i = 0; i < cnnW[dp].size(); i++) {
			for (int j = 0; j < cnnW[dp][0].size(); j++) {
				//cout << encoder.decode(cnnW[i][j]) << " " << endl;
				cout << encoder.decode(cnnW[dp][i][j]) << ((j < cnnW.at(dp).at(i).size()) ? ", " : ".\n");
			}
			cout << " row " << i << " finished" << endl;
		}
		cout << " depth " << dp << " finished" << endl;
	}
	*/

	/*
	// testing CNN3D with 8 nodes plus Fully connect layer
	int depth = 1;
	int nb_nodes = 8;
	vector<vector<vector<Ciphertext> > > input;
	input.emplace_back(x); // need to make x into tensor with x, y, z dimension
	// CNN3D -> output 2,2,8
	vector<vector<vector<Ciphertext> > > cnn3D_w1 = cnn3D(input, converted_cnn1W, evaluator, encoder, depth, nb_nodes, cnn1b, false); //no relu
	//Flatten _> output 32
	vector<Ciphertext> flat = flatten3D_cipher(cnn3D_w1);
	//prepare weights for 32rows and 10 cols
	int FcW_nb_rows = 32;
	vector<vector<double>> testWeight; // nb_rows = FcW_nb_rows; nb_cols == 10 
	for (int i = 0; i < FcW_nb_rows; i++) {
		testWeight.emplace_back(fcWeight[i]); // all cols needs to be copied;
	}
	vector<vector<Plaintext>> encodedTestFcW = encodeWB(testWeight, encoder);
	vector<vector<Plaintext>> encodedTestFcBias = encodeWB(fc_bias, encoder);
	// Fully connect
	vector<Ciphertext> fcResult = fullyConnectCiphertext(flat, encodedTestFcW, encodedTestFcBias, evaluator, false); // no relu at last layer

	cout << "printing CNN3D result:" << endl;
	for (int dp = 0; dp < cnn3D_w1.size(); dp++) {
		for (int i = 0; i < cnn3D_w1[dp].size(); i++) {
			for (int j = 0; j < cnn3D_w1[dp][0].size(); j++) {
				Plaintext cnnEleDecrypt = decryptor.decrypt(cnn3D_w1[dp][i][j]);
				cout << encoder.decode(cnnEleDecrypt) << ((j < cnn3D_w1.at(dp).at(i).size()) ? ", " : ".\n");
			}
			cout << " row " << i << " finished" << endl;
		}
		cout << " depth " << dp << " finished" << endl;
	}

	cout << "printing weights" << endl;
	for (int node = 0; node < nb_nodes; node++) {
		for (int dp = 0; dp < depth; dp++) {
			for (int i = 0; i < converted_cnn1W[node][dp].size(); i++) {
				for (int j = 0; j < converted_cnn1W[node][dp][0].size(); j++) {
					//cout << encoder.decode(cnnW[i][j]) << " " << endl;
					cout << encoder.decode(converted_cnn1W[node][dp][i][j]) << ((j < converted_cnn1W.at(node).at(dp).at(i).size()) ? ", " : ".\n");
				}
				cout << " row " << i << " finished" << endl;
			}
			cout << " depth " << dp << " finished" << endl;
		}
		cout << "node:" << node << "finished" << endl;
	}
	return fcResult;
	*/
}
vector<vector<vector<vector<Plaintext> > > > convertTensor(vector<vector<double> > &mat, FractionalEncoder encoder, int sizeX, int sizeY, int sizeZ, int nb_nodes) {
	// this function convert 2D tensor with (sizeX*sizeY * sizeZ , nb_nodes) into 4D tensor (nb_nodes, sizeZ, sizeX, sizeY)
	vector<vector<vector<vector<Plaintext> > > > output4D = {};// node, z, x, y
	cout << "has to be zeros, otherwise means matrix dimension doesn't match: " << mat.size() * mat[0].size() - sizeX * sizeY * sizeZ * nb_nodes << endl;
	for (int node = 0; node < nb_nodes; node++) {
		vector<vector<vector<Plaintext> > > tensor3D = {};
		for (int dp = 0; dp < sizeZ; dp++) {
			vector<vector<Plaintext> > tensor2D = {};
			for (int row = 0; row < sizeY; row++) {
				vector<Plaintext> tensor1D = {};
				for (int i = 0; i < sizeX; i++) {
					//int offset = dp * sizeX * sizeY + row * sizeX + i; // find the right idx from 2D to 4D;
					int offset = (row * sizeY + i) * sizeZ + dp;
					Plaintext encoded_number = encoder.encode(mat[offset][node]); // encode the element A[i][j]
					tensor1D.emplace_back(encoded_number);
				}
				tensor2D.emplace_back(tensor1D);
			}
			tensor3D.emplace_back(tensor2D);
		}
		output4D.emplace_back(tensor3D);
	}
	cout << output4D.size() << endl;
	cout << output4D[0].size() << endl;
	cout << output4D[0][0].size() << endl;
	cout << output4D[0][0][0].size() << endl;
	return output4D;
}
Ciphertext weightedSum2D(double row, double col, vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator) {
	// input x is a bigger matrix than w; both 2D, we want (row, col) result from these two
	vector<Ciphertext> encrypted_products;
	for (int i = row; i < row + w.size(); i++) {
		for (int j = col; j < col + w[0].size(); j++) {
			//sum = sum + x[i][j] * w[i - row][j - col];
			Ciphertext enc_plain_product = evaluator.multiply_plain(x[i][j], w[i - row][j - col]);
			encrypted_products.emplace_back(enc_plain_product);
		}
	}
	return evaluator.add_many(encrypted_products);
}
vector<vector<Ciphertext> > cnn2D(vector<vector<Ciphertext> > &x, vector<vector<Plaintext> > &w, Evaluator evaluator) {
	// this function only calculate 2D with depth = 1
	// initialize output cnn product
	int pSizeX = x.size() - w.size() + 1; // stride == 1
	int pSizeY = x[0].size() - w[0].size() + 1;
	vector<vector<Ciphertext> > product = {};

	for (int i = 0; i < pSizeX; i++) {
		vector<Ciphertext> pRow = {};
		for (int j = 0; j < pSizeY; j++) {
			pRow.emplace_back(weightedSum2D(i, j, x, w, evaluator));
		}
		product.emplace_back(pRow);
	}
	/*
	for (int row = 0; row < product.size(); row++) {
	//vector<double> productRow;
	for (int col = 0; col < product[0].size(); col++) {
	// Multiply the row of A by the column of B to get the row, column of product.
	product[row][col] = weightedSum2D(row, col, x, w, evaluator);
	}
	}
	*/
	return product;
}
vector<vector<Ciphertext> > compressCNN(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<Plaintext> > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, double nodeBias, bool relu) {
	vector<vector<vector<Ciphertext> > > tmp = {}; // calculate 2D CNN for every depth, push into tmp
	for (int i = 0; i < depth; i++) {
		tmp.emplace_back(cnn2D(x[i], w[i], evaluator));
	}
	vector<vector<Ciphertext> > cmp = tmp[0];
	
	Plaintext encoded_bias = encoder.encode(nodeBias);
	for (int row = 0; row < cmp.size(); row++) {
		for (int col = 0; col < cmp[0].size(); col++) {
			for (int i = 1; i < depth; i++) {
				//cmp[row][col] = cmp[row][col] + tmp[i][row][col]; // sum along depth;
				cmp[row][col] = evaluator.add(cmp[row][col], tmp[i][row][col]);
			}
			//cmp[row][col] = cmp[row][col] + nodeBias; // adding bias after all depth finished .* 
			cmp[row][col] = evaluator.add_plain(cmp[row][col], encoded_bias);
			if (relu) {
				//cmp[row][col] = max(0.0, cmp[row][col]);
				// relu function x^2
				cmp[row][col] = relu_appr(cmp[row][col], evaluator, encoder); // relu function, input ciphertext, output ciphertext
			}
		}
	}

	return cmp; // return compressed along z - axis
}

Ciphertext relu_appr(Ciphertext x, Evaluator evaluator, FractionalEncoder encoder) {
	// 0.1992 + 0.5002x + 0.1997x^2

	Plaintext degree0 = encoder.encode(0.1992);
	Ciphertext degree1 = evaluator.multiply_plain(x, encoder.encode(0.5002));
	Ciphertext degree01 = evaluator.add_plain(degree1, degree0);

	Ciphertext degree2 = evaluator.exponentiate(x, 2);
	degree2 = evaluator.multiply_plain(degree2, encoder.encode(0.1997));
	degree2 = evaluator.relinearize(degree2);

	vector<Ciphertext> encrypted_degrees;
	encrypted_degrees.emplace_back(degree01);
	encrypted_degrees.emplace_back(degree2);

	Ciphertext result = evaluator.add_many(encrypted_degrees);
	result = evaluator.relinearize(result);
	//cout << "ciphertext size" << to_string(result.size()) << endl;
	return result;

	//return evaluator.square(x); // befor just use square
}

vector<vector<vector<Ciphertext> > > cnn3D(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<vector<Plaintext> > > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, int nb_nodes, vector<vector<double> > b, bool relu) {
	// input x: 3D tensor (depth, x, y) ; w: 4D tensor (nb_nodes, depth, x, y) 
	// output is (nb_nodes, x, y)
	// for one CNN layer, for loop every node in weights matrix;
	vector<vector<vector<Ciphertext> > > product = {}; // NOTE: the depth of cnn output determined ONLY by number of CNN node; NOT x or w 's depth in other words, product.size() == number of nodes
	if (b[0].size() != 1) {
		cout << "nb_colomns for bias has to be 1, but it is = " << b[0].size() << endl;
	}
	for (int i = 0; i < nb_nodes; i++) {
		double nodeBias = b[i][0];
		vector<vector<Ciphertext> > cmp = compressCNN(x, w[i], evaluator, encoder, depth, nodeBias, relu); // compressCNN calculate ith node output;
		product.emplace_back(cmp);
	}
	cout << "One cnn3D is done" << endl;
	return product;
}
vector<Ciphertext> flatten3D_cipher(vector<vector<vector<Ciphertext> > > &matrix) {
	// match with keras flatten layer;
	int depth = matrix.size();
	int rowSize = matrix[0].size();
	int colSize = matrix[0][0].size();
	vector<Ciphertext> flat;
	for (int i = 0; i < rowSize; i++) {
		for (int j = 0; j < colSize; j++)
		{
			for (int k = 0; k < depth; k++)
			{
				flat.emplace_back(matrix[k][i][j]);
			}
		}
	}
	return flat;
}

vector<Ciphertext> fullyConnectCiphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, vector<vector<Plaintext>> &b, Evaluator evaluator, FractionalEncoder encodor, bool relu) {
	// x.shape 1 X m ; w.shape m X n ; b.shape n X 1 ; use b[i][0]
	if (x.size() != w.size()) {
		cout << "x size has to be equal to w.size; and w[0].size is the nb_neurons in this FC layer. This X.size = " << x.size() << "; w.size = " << w.size() << endl;
	}
	cout << "x.size = " << x.size() << endl;
	cout << "w.size" << w.size() << endl;
	cout << "w[0].size " << w[0].size() << endl;
	cout << "b.size: " << b.size() << endl;
	cout << "b[0].size " << b[0].size() << endl;
	vector<Ciphertext> dotProduct; // size will be nb_neurons
	int nb_neurons = w[0].size();
	for (int i = 0; i < nb_neurons; i++) {
		Plaintext bias_neuron = b[i][0];
		Ciphertext element = weightedSum1D_ciphertext(x, w, bias_neuron, evaluator, i); // ith col of w , .* x , then + b;
		if (relu) {
			element = relu_appr(element, evaluator, encodor);
		}
		dotProduct.emplace_back(element);
	}
	return dotProduct;
}
Ciphertext weightedSum1D_ciphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, Plaintext bias, Evaluator evaluator, int col) {
	vector<Ciphertext> products;
	for (int i = 0; i < x.size(); i++) {
		//sum = sum + x[i] * w[i][col];
		products.emplace_back(evaluator.add_plain(x[i], w[i][col]));
	}	
	Ciphertext sum = evaluator.add_many(products);
	//sum = sum + b[col][0];
	sum = evaluator.add_plain(sum, bias);
	return sum;
}