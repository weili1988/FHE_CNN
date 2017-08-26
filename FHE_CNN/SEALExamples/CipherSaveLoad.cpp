#include "CipherSaveLoad.h"
void ciphertextSave2D(vector<vector<Ciphertext>> &x, string path) {
	int rowSize = x.size();
	int colSize = x[0].size();

	for (int row = 0; row < x.size(); row++)
		for (int col = 0; col < x[0].size(); col++) {
		    string file = path;
			file.append("r");
			file.append(to_string(row));
			file.append("c");
			file.append(to_string(col));
			ofstream output(file, std::ios::out | std::ofstream::binary);
			x[row][col].save(output);
			printf("Data at row%d, col%d saved\n", row, col);
		}

}

void ciphertextSave3D(vector<vector<vector<Ciphertext>>> &x, string path) {
	//input x[depth,x,y]
	path.append("d");
	for (int depth = 0; depth < x.size(); depth++) {
		string file = path;
		file.append(to_string(depth));
		printf("Saving depth %d", depth);
		ciphertextSave2D(x[depth], file);
	}
}



void ciphertextLoad2D(vector<vector<Ciphertext>> &x, string path, int rowSize, int colSize) {
	for (int row = 0; row < rowSize; row++) {
		vector<Ciphertext> cipher1D;
		for (int col = 0; col < colSize; col++) {
		    string file = path;
			file.append("r");
			file.append(to_string(row));
			file.append("c");
			file.append(to_string(col));
			ifstream input(file, std::ios::in | std::ifstream::binary);
			Ciphertext tmp = Ciphertext::Ciphertext(BigPolyArray());
			tmp.load(input);
			cipher1D.emplace_back(tmp);
			printf("Data at row%d, col%d loaded\n", row, col);
		}
		x.emplace_back(cipher1D);
	}

}

void ciphertextLoad3D(vector<vector<vector<Ciphertext>>> &x, string path, int depthSize, int rowSize, int colSize) {
	path.append("d");
	for (int depth = 0; depth < depthSize; depth++) {
		vector<vector<Ciphertext>> cipher2D;
		string file = path;
		file.append(to_string(depth));
		printf("Loading depth %d", depth);
		ciphertextLoad2D(cipher2D, file, rowSize, colSize);
	}
}

void saveKeys(Ciphertext publicKey, Plaintext secretKey, EvaluationKeys evaluationKey, string path) {
	string file1, file2, file3;
	file1 = path;
	file1.append("publickey");
	ofstream output1(file1, std::ios::out | std::ofstream::binary);
	publicKey.save(output1);
	file2 = path;
	file2.append("secretkey");
	ofstream output2(file2, std::ios::out | std::ofstream::binary);
	secretKey.save(output2);
	file3 = path;
	file3.append("evaluationkeys");
	ofstream output3(file3, std::ios::out | std::ofstream::binary);
	evaluationKey.save(output3);
}

void loadKeys(Ciphertext &publicKey, Plaintext &secretKey, EvaluationKeys &evaluationKey, string path) {
	string file1, file2, file3;
	file1 = path;
	file1.append("publickey");
	ifstream input1(file1, std::ios::in | std::ifstream::binary);
	publicKey.load(input1);
	file2 = path;
	file2.append("secretkey");
	ifstream input2(file2, std::ios::in | std::ifstream::binary);
	secretKey.load(input2);
	file3 = path;
	file3.append("evaluationkeys");
	ifstream input3(file3, std::ios::in | std::ifstream::binary);
	evaluationKey.load(input3);
}

void saveParms(EncryptionParameters parms, string path) {
	path.append("parms");
	ofstream output(path, std::ios::out | std::ofstream::binary);
	parms.save(output);
}

void loadParms(EncryptionParameters &parms, string path) {
	path.append("parms");
	ifstream input(path, std::ios::in | std::ifstream::binary);
	parms.load(input);
}