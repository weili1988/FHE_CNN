#include "AveragePooling.h"

Ciphertext avgSum(vector<vector<Ciphertext> > &x, int poolSize, int row, int col, FractionalEncoder encoder, Evaluator evaluator) {
	Ciphertext sum = Ciphertext::Ciphertext(BigPolyArray());
	sum = x[row][col];
	for (int i = 0; i < poolSize; i++) {
		for (int j = 0; j < poolSize; j++) {
			if ((i == 0) && (j == 0))
				continue;
			sum = evaluator.add(sum, x[row + i][col + j]);
			//cout << "i:" << i << " j:" << j << endl;
		}
	}
	double numerator = 1.0 / (poolSize * poolSize);
	Plaintext numPlain = encoder.encode(numerator);
	sum = evaluator.multiply_plain(sum, numPlain);
	return sum;
}

vector<vector<Ciphertext> > averagePooling2D(vector<vector<Ciphertext> > &x, int poolSize, int strideSize, FractionalEncoder encoder, Evaluator evaluator) {
	int rsize = (x.size() - poolSize) / strideSize + 1;
	int csize = (x[0].size() - poolSize) / strideSize + 1;
	vector<vector<Ciphertext> > pooled;
	for (int row = 0; row < rsize; row++) {
		vector<Ciphertext> pooledRow;
		for (int col = 0; col < csize; col++) {
			Ciphertext avg = avgSum(x, poolSize, row * strideSize, col * strideSize, encoder, evaluator);
			pooledRow.emplace_back(avg);
		}
		pooled.emplace_back(pooledRow);
	}
	return pooled;
}

vector<vector<vector<Ciphertext> > > averagePooling3D(vector<vector<vector<Ciphertext> > > &x, int poolSize, int strideSize, FractionalEncoder encoder, Evaluator evaluator) {
	// input x has (depth, x, y)
	// output has depth untouched, size of x = (original - poolSize) / strideSize + 1; Assuming square pool and stride
	vector<vector<vector<Ciphertext> > > res;
	for (int i = 0; i < x.size(); i++) {
		res.emplace_back(averagePooling2D(x[i], poolSize, strideSize, encoder, evaluator));
	}
	return res;
}