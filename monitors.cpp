#include "monitors.h"
// Plain text
// input: single plaintxt 
// output: a string
string returnPlain(Plaintext pt) {
	string res;
	res = pt.to_string();
	return res;
}

// input: vector of plaintxt
// output: vector of strings
vector<string> returnPlain(vector<Plaintext> pt) {
	vector<string> res1D;

	for (int i = 0; i < pt.size(); i++) {
		res1D.emplace_back(returnPlain(pt[i]));
	}

	return res1D;
}

// input: 2D vector of plaintxt
// output: 2D vector of strings
vector<vector<string>> returnPlain(vector<vector<Plaintext>> pt) {
	vector<vector<string>> res2D;

	for (int i = 0; i < pt.size(); i++) {
		res2D.emplace_back(returnPlain(pt[i]));
	}

	return res2D;
}

// save single plaintxt into a txt file
void savePlain(Plaintext pt, string file) {
	string res;
	res = returnPlain(pt);

	// save to a text file
	ofstream outputFile(file);
	outputFile << res << endl;
}

// save vector of plaintxt
void savePlain(vector<Plaintext> pt, string file) {
	vector<string> res1D;
	res1D = returnPlain(pt);

	// save to a text file
	ofstream outputFile(file);
	for (int row = 0; row < res1D.size(); row++) {
		outputFile << "r" << row << endl;
		outputFile << res1D[row] << endl;
	}
}

// save 2d vector of plaintxt
void savePlain(vector<vector<Plaintext>> pt, string file) {
	vector<vector<string>> res2D;
	res2D = returnPlain(pt);

	// save to a text file
	ofstream outputFile(file);
	for (int row = 0; row < res2D.size(); row++)
		for (int col = 0; col < res2D[row].size(); col++) {
			outputFile << "r" << row << "c" << col << endl;
			outputFile << res2D[row][col] << endl;
		}
}

// Cipher text

// input: single ciphertxt 
// output: a vector of strings
vector<string> returnCipher(Ciphertext ct) {
	vector<string> res;
	BigPolyArray bpact(ct);
	
	for (int i = 0; i < ct.size(); i++) {
		string tmp = bpact[i].to_string();
		res.emplace_back(tmp);
		// cout << i << ": " << tmp << endl;

	}
	return res;
}

// input: vector of ciphertxt
// output: vector of vector of strings
vector<vector<string>> returnCipher(vector<Ciphertext> ct) {
	vector<vector<string>> res1D;

	for (int i = 0; i < ct.size(); i++) {
		res1D.emplace_back(returnCipher(ct[i]));
	}

	return res1D;
}

// input: vector of ciphertxt
// output: vector of vector of strings
vector<vector<vector<string>>> returnCipher(vector<vector<Ciphertext>> ct) {
	vector<vector<vector<string>>> res2D;

	for (int i = 0; i < ct.size(); i++) {
		res2D.emplace_back(returnCipher(ct[i]));
	}

	return res2D;
}

// save single ciphertxt
void saveCipher(Ciphertext ct, string file) {
	vector<string> res;
	res = returnCipher(ct);

	// save to a text file
	ofstream outputFile(file);
	outputFile << "size: " << res.size() << endl;
	for (const auto &row : res)
		outputFile << row << endl;
}

// save vector of ciphertxt
void saveCipher(vector<Ciphertext> ct, string file) {
	vector<vector<string>> res1D;
	res1D = returnCipher(ct);

	// save to a text file
	ofstream outputFile(file);
	for (int row = 0; row < res1D.size(); row++) {
		outputFile << "r" << row << " Size: " << res1D[row].size() << endl;
		for (const auto &col : res1D[row])
			outputFile << col << endl;
	}
}

// save 2d vector of ciphertxt
void saveCipher(vector<vector<Ciphertext>> ct, string file) {
	vector<vector<vector<string>>> res2D;
	res2D = returnCipher(ct);

	// save to a text file
	ofstream outputFile(file);
	for (int row = 0; row < res2D.size(); row++) 
	for (int col = 0; col < res2D[row].size(); col++) {
		outputFile << "r" << row << "c" << col << " Size: " << res2D[row][col].size() << endl;
		for (const auto &col : res2D[row][col])
			outputFile << col << endl;
	}
}