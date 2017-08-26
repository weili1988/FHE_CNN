#include <iostream>
#include <algorithm>    // std::max
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
vector<vector<double> > cnn2D(vector<vector<double> > &x, vector<vector<double> > &w);
void print2D(vector<vector<double> > &nums);
void print1D(vector<double> & nums);
double weightedSum(double row, double col, vector<vector<double> > &x, vector<vector<double> > &w);
//vector<vector<double> > compressCNN(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth);
//vector<vector<double> > compressCNN(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth, double nodeBias);
vector<vector<double> > compressCNN(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth, double nodeBias, bool relu);
//vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth, double nb_nodes);
//vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<vector<double> > > > &w, double depth, double nb_nodes);
//vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<vector<double> > > > &w, double depth, double nb_nodes, vector<vector<double> > b);
vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<vector<double> > > > &w, double depth, double nb_nodes, vector<vector<double> > b, bool relu);
vector<vector<double> > readCSV(string pathFile);
vector<vector<vector<vector<double> > > > convertTensor(vector<vector<double> > &mat, int sizeX, int sizeY, int sizeZ, int nb_nodes);
vector<vector<vector<double> > > averagePooling3D(vector<vector<vector<double> > > &x, int poolSize, int strideSize);
vector<vector<double> > averagePooling2D(vector<vector<double> > &x, int poolSize, int strideSize);
double avgSum(vector<vector<double> > &x, int poolSize, int row, int col);
vector<double> flatten3D(vector<vector<vector<double> > > &matrix);
vector<double> fullyConnect(vector<double> &x, vector<vector<double> > &w, vector<vector<double>> &b, bool relu);
double weightedSum1D(vector<double> &x, vector<vector<double> > &w, vector<vector<double>> &b, int col);
int main()
{

	string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights";
	string x0, w1, b1, w2, b2, fcW, fc_b;
	x0.append(path);
	x0.append("/7.csv");
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
	vector<vector<double> > input = readCSV(x0);
	vector<vector<vector<double> > > x = {};
	x.emplace_back(input);

	vector<vector<double> > cnn1W = readCSV(w1);
	vector<vector<double> > cnn1b = readCSV(b1);
	vector<vector<double> > cnn2W = readCSV(w2);
	vector<vector<double> > cnn2b = readCSV(b2);
	vector<vector<double> > fcWeight = readCSV(fcW);
	vector<vector<double> > fc_bias = readCSV(fc_b);
	
	vector<vector<vector<vector<double> > > > converted_cnn1W = convertTensor(cnn1W, 3, 3, 1, 8);
	int nb_nodesCNN1W = 8;
	int depth_cnn1W = 1;
	vector<vector<vector<double> > > cnn1_product = cnn3D(x, converted_cnn1W, depth_cnn1W, nb_nodesCNN1W, cnn1b, true);
	
	int poolsize1 = 2;
	int stridesize1 = 2;
	vector<vector<vector<double> > > cnn1_pooled = averagePooling3D(cnn1_product, poolsize1, stridesize1);

	int nb_nodesCNN2W = 16;
	int depth_cnn2W = 8;
	vector<vector<vector<vector<double> > > > converted_cnn2W = convertTensor(cnn2W, 3, 3, 8, 16);
	
	vector<vector<vector<double> > > cnn2_product = cnn3D(cnn1_pooled, converted_cnn2W, depth_cnn2W, nb_nodesCNN2W, cnn2b, true);
	
	int poolsize2 = 4;
	int stridesize2 = 4;
	vector<vector<vector<double> > > cnn2_pooled = averagePooling3D(cnn2_product, poolsize2, stridesize2);

	vector<double> flat = flatten3D(cnn2_pooled);

	vector<double> result = fullyConnect(flat, fcWeight, fc_bias, false);

	print2D(cnn1_product[0]);
	print2D(cnn1_pooled[0]);
	print2D(cnn2_product[0]);

	print1D(flat);
	print1D(result);

	return 0;
}

vector<vector<vector<vector<double> > > > convertTensor(vector<vector<double> > &mat, int sizeX, int sizeY, int sizeZ, int nb_nodes) {
	// this function convert 2D tensor with (sizeX*sizeY * sizeZ , nb_nodes) into 4D tensor (nb_nodes, sizeZ, sizeX, sizeY)
	vector<vector<vector<vector<double> > > > output4D = {};// node, z, x, y
	cout << "has to be zeros, otherwise means matrix dimension doesn't match: " << mat.size() * mat[0].size() - sizeX * sizeY * sizeZ * nb_nodes << endl;
	for (int node = 0; node < nb_nodes; node++) {
		vector<vector<vector<double> > > tensor3D = {};
		for (int dp = 0; dp < sizeZ; dp++) {
			vector<vector<double> > tensor2D = {};
			for (int row = 0; row < sizeY; row++) {
				vector<double> tensor1D = {};
				for (int i = 0; i < sizeX; i++) {
					//int offset = dp * sizeX * sizeY + row * sizeX + i; // find the right idx from 2D to 4D;
					int offset = (row * sizeY + i) * sizeZ + dp;
					tensor1D.emplace_back(mat[offset][node]);
					//output4D[node][dp][row][i] = dp * sizeZ + ;
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
/*
vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth, double nb_nodes) {
// input x: 3D tensor (depth, x, y) ; w: 3D tensor (depth, x, y)
// for one CNN layer, for loop every node in weights matrix;
vector<vector<vector<double> > > product = {}; // NOTE: the depth of cnn output determined ONLY by number of CNN node; NOT x or w 's depth
// in other words, product.size() == number of nodes
for (int i = 0; i < nb_nodes; i++) {
//vector<vector<double> > cmp = compriseCNN(vector<vector<vector<double> > > x, vector<vector<vector<double> > > w[i], double depth) // compress
vector<vector<double> > cmp = compressCNN(x, w, depth); // needs to adjust to w[i] i for ith node;
product.emplace_back(cmp);
}

return product;
}
*/

vector<vector<vector<double> > > cnn3D(vector<vector<vector<double> > > &x, vector<vector<vector<vector<double> > > > &w, double depth, double nb_nodes, vector<vector<double> > b, bool relu) {
	// input x: 3D tensor (depth, x, y) ; w: 4D tensor (nb_nodes, depth, x, y) 
	// output is (nb_nodes, x, y)
	// for one CNN layer, for loop every node in weights matrix;
	vector<vector<vector<double> > > product = {}; // NOTE: the depth of cnn output determined ONLY by number of CNN node; NOT x or w 's depth
												   // in other words, product.size() == number of nodes
	if (b[0].size() != 1) {
		cout << "nb_colomns for bias has to be 1, but it is = " << b[0].size() << endl;
	}
	
	for (int i = 0; i < nb_nodes; i++) {
		double nodeBias = b[i][0];
		vector<vector<double> > cmp = compressCNN(x, w[i], depth, nodeBias, relu); // compressCNN calculate ith node output;
		product.emplace_back(cmp);
	}

	return product;
}

vector<vector<double> > compressCNN(vector<vector<vector<double> > > &x, vector<vector<vector<double> > > &w, double depth, double nodeBias, bool relu) {
	vector<vector<vector<double> > > tmp = {}; // calculate 2D CNN for every depth, push into tmp
	for (int i = 0; i < depth; i++) {
		tmp.emplace_back(cnn2D(x[i], w[i]));
	}
	vector<vector<double> > cmp = tmp[0];
	/*
	for (int i = 1; i < depth; i++) {
		for (int row = 0; row < cmp.size(); row++) {
			for (int col = 0; col < cmp[0].size(); col++) {
				cmp[row][col] = cmp[row][col] + tmp[i][row][col];
			}
		}
	}
	*/
	for (int row = 0; row < cmp.size(); row++) {
		for (int col = 0; col < cmp[0].size(); col++) {
			for (int i = 1; i < depth; i++) {
				cmp[row][col] = cmp[row][col] + tmp[i][row][col]; // sum along depth;
			}
			cmp[row][col] = cmp[row][col] + nodeBias; // adding bias after all depth finished .* 
			if (relu) {
				cmp[row][col] = max(0.0, cmp[row][col]);
			}
		}
	}

	return cmp; // return compressed along z - axis
}

vector<vector<double> > cnn2D(vector<vector<double> > &x, vector<vector<double> > &w) {
	// this function only calculate 2D with depth = 1
	
	// initialize output cnn product
	int pSizeX = x.size() - w.size() + 1; // stride == 1
	int pSizeY = x[0].size() - w[0].size() + 1;
	vector<vector<double> > product = {};

	for (int i = 0; i < pSizeX; i++) {
		vector<double> pRow = {};
		for (int j = 0; j < pSizeY; j++) {
			pRow.emplace_back(0);
		}
		product.emplace_back(pRow);
	}

	for (int row = 0; row < product.size(); row++) {
		//vector<double> productRow;
		for (int col = 0; col < product[0].size(); col++) {
			// Multiply the row of A by the column of B to get the row, column of product.
			product[row][col] = weightedSum(row, col, x, w);
		}
	}
	return product;
}
double weightedSum(double row, double col, vector<vector<double> > &x, vector<vector<double> > &w) {
	double res = 0;
	for (int i = row; i < row + w.size(); i++) {
		for (int j = col; j < col + w[0].size(); j++) {
			res += x[i][j] * w[i - row][j - col];
		}
	}		
	return res;
}
void print2D(vector<vector<double> > &nums) {
	double sizeX = nums.size();
	double sizeY = nums[0].size();

	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < sizeY; j++) {
			cout << nums[i][j] << " ";
		}
		cout << "\n";
	}
}
void print1D(vector<double> & nums) {
	for (int i = 0; i < nums.size(); i++) {
		cout << nums[i] << " ";
	}
	cout << "\n";
}
vector<vector<double> > readCSV(string pathFile)
{
	// read csv file doubleo values vectors, note that csv file is separated with ',' and note csv is ?2D?
	vector<vector<double> > values;
	vector<double> valueline;
	ifstream fin(pathFile);
	string item;
	for (string line; getline(fin, line); )
	{
		istringstream in(line);

		while (getline(in, item, ','))
		{
			valueline.push_back(atof(item.c_str()));
		}
		values.push_back(valueline);
		valueline.clear();
	}
	return values;
}
vector<vector<vector<double> > > averagePooling3D(vector<vector<vector<double> > > &x, int poolSize, int strideSize) {
	// input x has (depth, x, y)
	// output has depth untouched, size of x = (original - poolSize) / strideSize + 1; Assuming square pool and stride
	vector<vector<vector<double> > > res;
	for (int i = 0; i < x.size(); i++) {
		res.emplace_back(averagePooling2D(x[i], poolSize, strideSize));
	}
	return res;
}
vector<vector<double> > averagePooling2D(vector<vector<double> > &x, int poolSize, int strideSize) {
	int size = (x.size() - poolSize) / strideSize + 1;
	vector<vector<double> > pooled;
	for (int row = 0; row < size; row++) {
		vector<double> pooledRow;
		for (int col = 0; col < size; col++) {
			double avg = avgSum(x, poolSize, row * strideSize, col * strideSize);
			pooledRow.emplace_back(avg);
		}
		pooled.emplace_back(pooledRow);
	}
	return pooled;
}
double avgSum(vector<vector<double> > &x, int poolSize, int row, int col) {
	double sum = 0;
	for (int i = 0; i < poolSize; i++) {
		for (int j = 0; j < poolSize; j++) {
			sum = sum + x[row + i][col + j];
			//cout << "i:" << i << " j:" << j << endl;
		}
	}
	return sum / (poolSize * poolSize);
}
vector<double> flatten3D(vector<vector<vector<double> > > &matrix) {
	int depth = matrix.size();
	int rowSize = matrix[0].size();
	int colSize = matrix[0][0].size();
	vector<double> flat;
	/*
	for (int i = 0; i < depth; i++) {
		for (int j = 0; j < rowSize; j++)
		{
			for (int k = 0; k < colSize; k++)
			{
				flat.emplace_back(matrix[i][j][k]);
			}
		}
	}
	*/
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
vector<double> fullyConnect(vector<double> &x, vector<vector<double> > &w, vector<vector<double>> &b, bool relu) {
	// x.shape 1 X m ; w.shape m X n ; b.shape n X 1 ; use b[i][0]
	if (x.size() != w.size()) {
		cout << "x size has to be equal to w.size; and w[0].size is the nb_neurons in this FC layer. This X.size = " << x.size() << "; w.size = " << w.size() << endl;
	}
	cout << "x.size = " << x.size() << endl;
	cout << "w.size" << w.size() << endl;
	cout << "w[0].size " << w[0].size() << endl;
	cout << "b.size: " << b.size() << endl;
	cout << "b[0].size " << b[0].size() << endl;
	vector<double> dotProduct;
	int nb_neurons = w[0].size();
	for (int i = 0; i < nb_neurons; i++) {
		double element = weightedSum1D(x, w, b, i); // ith col of w , .* x , then + b;
		if (relu) {
			element = max(0.0, element);
		}
		dotProduct.emplace_back(element);
	}
	return dotProduct;
}
double weightedSum1D(vector<double> &x, vector<vector<double> > &w, vector<vector<double>> &b, int col) {
	double sum = 0.0;
	for (int i = 0; i < x.size(); i++) {
		sum = sum + x[i] * w[i][col];
	}
	sum = sum + b[col][0];
	return sum;
}