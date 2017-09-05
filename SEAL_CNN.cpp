#include "SEAL_CNN.h"

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
	
	//real test of 2 layers CNN3D + 1 FC
	int nb_nodesCNN1W = 8;
	int depth_cnn1W = 1;
	vector<vector<vector<Ciphertext> > > input;
	input.emplace_back(x); // need to make x into tensor with x, y, z dimension
	vector<vector<vector<Ciphertext> > > cnn1_product = cnn3D(input, converted_cnn1W, evaluator, encoder, depth_cnn1W, nb_nodesCNN1W, cnn1b, true, true); // relu, with relin

	int poolsize1 = 2;
	int stridesize1 = 2;
	vector<vector<vector<Ciphertext> > > cnn1_pooled = averagePooling3D(cnn1_product, poolsize1, stridesize1, encoder, evaluator);
	printCipher2D(cnn1_pooled.at(0), encoder, decryptor);
	int nb_nodesCNN2W = 16;
	int depth_cnn2W = 8;
	vector<vector<vector<vector<Plaintext> > > > converted_cnn2W = convertTensor(cnn2W, encoder, 3, 3, 8, 16);
	vector<vector<vector<Ciphertext> > > cnn2_product = cnn3D(cnn1_pooled, converted_cnn2W, evaluator, encoder, depth_cnn2W, nb_nodesCNN2W, cnn2b, true, false); // relu, withou relin

	int poolsize2 = 4;
	int stridesize2 = 4;
	vector<vector<vector<Ciphertext> > > cnn2_pooled = averagePooling3D(cnn2_product, poolsize2, stridesize2, encoder, evaluator);

	//This is to flat the cnn result
	vector<Ciphertext> flat = flatten3D_cipher(cnn2_pooled);
	printFlat1D(flat, encoder, decryptor);
	//fully connect layer
	vector<vector<Plaintext>> encodedFcW = encodeWB(fcWeight, encoder);
	vector<vector<Plaintext>> encodedFcBias = encodeWB(fc_bias, encoder);
	vector<Ciphertext> fcResult = fullyConnectCiphertext(flat, encodedFcW, encodedFcBias, evaluator, encoder, false); // no relu at last layer, no relin either
	return fcResult;
}
vector<Ciphertext> cnnCalMnist_slim(vector<vector<Ciphertext>> x, FractionalEncoder encoder, Evaluator evaluator, Decryptor decryptor) { // remember to delete this decrypter after development done
																																	// This is server
																																	// has encrypted x, encoder, evaluator, (and weights in server); output the vecotor (10 ciphter text) back to user who can decrypt the ciphertext back to plaintext
	string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights_test";
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

	vector<vector<vector<vector<Plaintext> > > > converted_cnn1W = convertTensor(cnn1W, encoder, 3, 3, 1, 1); // convert to tensor with shape (8, 1, 3, 3) and finished encoding

																											  //real test of 2 layers CNN3D + 1 FC
	int nb_nodesCNN1W = 1;
	int depth_cnn1W = 1;
	vector<vector<vector<Ciphertext> > > input;
	input.emplace_back(x); // need to make x into tensor with x, y, z dimension
	vector<vector<vector<Ciphertext> > > cnn1_product = cnn3D(input, converted_cnn1W, evaluator, encoder, depth_cnn1W, nb_nodesCNN1W, cnn1b, true, true); // relu, with relin

	int poolsize1 = 2;
	int stridesize1 = 2;
	vector<vector<vector<Ciphertext> > > cnn1_pooled = averagePooling3D(cnn1_product, poolsize1, stridesize1, encoder, evaluator);
	printCipher2D(cnn1_pooled.at(0), encoder, decryptor);
	int nb_nodesCNN2W = 1;
	int depth_cnn2W = 1;
	vector<vector<vector<vector<Plaintext> > > > converted_cnn2W = convertTensor(cnn2W, encoder, 3, 3, 1, 1);
	vector<vector<vector<Ciphertext> > > cnn2_product = cnn3D(cnn1_pooled, converted_cnn2W, evaluator, encoder, depth_cnn2W, nb_nodesCNN2W, cnn2b, true, false); // relu, withou relin

	int poolsize2 = 2;
	int stridesize2 = 2;
	vector<vector<vector<Ciphertext> > > cnn2_pooled = averagePooling3D(cnn2_product, poolsize2, stridesize2, encoder, evaluator);

	//This is to flat the cnn result
	vector<Ciphertext> flat = flatten3D_cipher(cnn2_pooled);
	printFlat1D(flat, encoder, decryptor);
	//fully connect layer
	vector<vector<Plaintext>> encodedFcW = encodeWB(fcWeight, encoder);
	vector<vector<Plaintext>> encodedFcBias = encodeWB(fc_bias, encoder);
	vector<Ciphertext> fcResult = fullyConnectCiphertext(flat, encodedFcW, encodedFcBias, evaluator, encoder, false); // no relu at last layer, no relin either
	return fcResult;
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
vector<vector<Ciphertext> > compressCNN(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<Plaintext> > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, double nodeBias, bool relu, bool relu_relin) {
	vector<vector<vector<Ciphertext> > > tmp = {}; // calculate 2D CNN for every depth, push into tmp
	for (int i = 0; i < depth; i++) {
		tmp.emplace_back(cnn2D(x[i], w[i], evaluator));
		cout << "finished 1 CNN2D" << endl;
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
				cmp[row][col] = relu_appr(cmp[row][col], evaluator, encoder, relu_relin); // relu function, input ciphertext, output ciphertext
			}
		}
	}

	return cmp; // return compressed along z - axis
}

Ciphertext relu_appr(Ciphertext x, Evaluator evaluator, FractionalEncoder encoder, bool relin) {
	/*
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
	*/
	Ciphertext result = evaluator.square(x);
	if (relin) {
		result = evaluator.relinearize(result);
	}
	return result; // befor just use square
}

vector<vector<vector<Ciphertext> > > cnn3D(vector<vector<vector<Ciphertext> > > &x, vector<vector<vector<vector<Plaintext> > > > &w, Evaluator evaluator, FractionalEncoder encoder, int depth, int nb_nodes, vector<vector<double> > b, bool relu, bool relu_relin) {
	// input x: 3D tensor (depth, x, y) ; w: 4D tensor (nb_nodes, depth, x, y) 
	// output is (nb_nodes, x, y)
	// for one CNN layer, for loop every node in weights matrix;
	vector<vector<vector<Ciphertext> > > product = {}; // NOTE: the depth of cnn output determined ONLY by number of CNN node; NOT x or w 's depth in other words, product.size() == number of nodes
	if (b[0].size() != 1) {
		cout << "nb_colomns for bias has to be 1, but it is = " << b[0].size() << endl;
	}
	for (int i = 0; i < nb_nodes; i++) {
		double nodeBias = b[i][0];
		vector<vector<Ciphertext> > cmp = compressCNN(x, w[i], evaluator, encoder, depth, nodeBias, relu, relu_relin); // compressCNN calculate ith node output;
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
			element = relu_appr(element, evaluator, encodor, false); // no need for relinearization at last layer
		}
		dotProduct.emplace_back(element);
	}
	return dotProduct;
}
Ciphertext weightedSum1D_ciphertext(vector<Ciphertext> &x, vector<vector<Plaintext> > &w, Plaintext bias, Evaluator evaluator, int col) {
	vector<Ciphertext> products;
	for (int i = 0; i < x.size(); i++) {
		//sum = sum + x[i] * w[i][col];
		products.emplace_back(evaluator.multiply_plain(x[i], w[i][col])); // .*
	}	
	Ciphertext sum = evaluator.add_many(products);
	//sum = sum + b[col][0];
	sum = evaluator.add_plain(sum, bias);
	return sum;
}
void printFlat1D(vector<Ciphertext> flat, FractionalEncoder encoder, Decryptor decryptor) {
	int size = flat.size();
	cout << "printing flat, the input of FC" << endl;
	for (int i = 0; i < size; i++) {
		Plaintext tmp = decryptor.decrypt(flat.at(i));
		cout << encoder.decode(tmp) << " ";
	}
	cout << "finished printing 1D flat" << endl;
}
void printCipher2D(vector<vector<Ciphertext> > mat, FractionalEncoder encoder, Decryptor decryptor) {
	int n_rows = mat.size();
	int n_cols = mat.at(0).size();
	for (int i = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++) {
			Plaintext tmp = decryptor.decrypt(mat.at(i).at(j));
			cout << encoder.decode(tmp) << " ";
		}
		cout << "finished " << i << "row" << endl;
	}
	cout << "finished 2D print" << endl;
}