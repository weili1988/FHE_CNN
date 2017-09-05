#include "utility.h"

int logitstic_regression() {

	print_example_banner("Logistic regression");

	string pathFile = "./data/data.csv";
	vector<vector<double> > loadValues = readCSV(pathFile); //load from csv into 2D vector

	//// encryption necessarry encoder and encryptor 
	// Create encryption parameters
	EncryptionParameters parms;
	parms.set_poly_modulus("1x^2048 + 1");
	parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));
	parms.set_plain_modulus(1 << 8);
	parms.validate();

	// Generate keys.
	cout << "Generating keys ..." << endl;
	KeyGenerator generator(parms);
	generator.generate();
	cout << "... key generation complete" << endl;
	Ciphertext public_key = generator.public_key();
	Plaintext secret_key = generator.secret_key();

	//FractionalEncoder to encode
	FractionalEncoder encoder(parms.plain_modulus(), parms.poly_modulus(), 64, 32, 3);

	// Create the rest of the tools
	Encryptor encryptor(parms, public_key);
	Evaluator evaluator(parms);
	Decryptor decryptor(parms, secret_key);

	//// encrypt every number of a 2D array into cipertext 
	vector<vector<Ciphertext>> encrypted2D = encrypt(loadValues, encoder, encryptor); //
	
	//// decrypt 
	// Decrypt to test the accuracy ; 
	cout << "Decrypting 3rd row 3rd col... " << endl;
	Plaintext plain_result = decryptor.decrypt(encrypted2D.at(1).at(1));
	cout << "finished decrypting 3rd row 3rd col ..." << endl;

	// Print the result of 3rd row 3rd col
	double result = encoder.decode(plain_result);
	cout << "3rd row 3rd col: " << result << endl;

	cout << "loadValues at 3rd row 3rd col:";
	cout << loadValues.at(1).at(1) << endl;

	return 1;
}

vector<vector<double> > readCSV(string pathFile)
{
	// read csv file into values vectors, note that csv file is separated with ',' and note csv is ?2D?
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

vector<vector<Ciphertext>> encrypt(vector<vector<double> > loadValues, FractionalEncoder encoder, Encryptor encryptor) {
	// traverse loadValues which is a 2D vector, encrypt all into Ciphertext;
	vector<vector<Ciphertext> > encrypted; // this is to be returned 

	for (int i = 0; i < loadValues.size(); i++) {
		//loadValues.at(i); // encrypt row i

		///
		cout << "Encrypting ... ";
		vector<Ciphertext> encrypted_rationals;
		for (int j = 0; j < loadValues.at(i).size(); j++)
		{
			Plaintext encoded_number = encoder.encode(loadValues.at(i).at(j)); // encode the element A[i][j]
			encrypted_rationals.emplace_back(encryptor.encrypt(encoded_number));
			//cout << to_string(loadValues.at(i).at(j)) << ((j < loadValues.at(i).size()) ? ", " : ".\n");
		}
		encrypted.emplace_back(encrypted_rationals);
		//cout << "finished encrypting " << to_string(i) << " row" << endl;
	}
	return encrypted;
}

void testLoading(vector<vector<double> > dataX) {
	// test if the loading is correct
	int m = dataX.size();
	int n = dataX.at(0).size();

	for(int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << to_string(dataX.at(i).at(j)) << ((j < dataX.at(i).size()) ? ", " : ".\n");
		}
	}
}

int testfullyConnect() {

	print_example_banner("Fully connected layer");
	
	/// loading data
	cout << "loading data" << endl;
	string pathDataFile = "./data/data.csv";
	string pathTargetFile = "./data/target.csv";

	string pathWeightsFile = "./data/weights.csv";
	string pathBiasFile = "./data/bias.csv";

	string pathWeightsOutputFile = "./data/weights_output.csv";
	string pathBiasOutputFile = "./data/bias_output.csv";

	vector<vector<double> > data = readCSV(pathDataFile); //load from csv into 2D vector
	vector<vector<double> > tgt = readCSV(pathTargetFile); //load from csv into 2D vector

	vector<vector<double> > weights = readCSV(pathWeightsFile); //load from csv into 2D vector
	vector<vector<double> > bias = readCSV(pathBiasFile); //load from csv into 2D vector

	vector<vector<double> > weightsOutput = readCSV(pathWeightsOutputFile); //load from csv into 2D vector
	vector<vector<double> > biasOutput = readCSV(pathBiasOutputFile); //load from csv into 2D vector

	//testLoading(data); this is to print the loaded data

	// Create encryption parameters
	EncryptionParameters parms;
	parms.set_poly_modulus("1x^4096 + 1");
	parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
	parms.set_plain_modulus(1 << 6);
	parms.set_decomposition_bit_count(16);
	parms.validate(); // parms.validate() is the end of parms creation

	// Generate keys.
	cout << "Generating keys ..." << endl;
	KeyGenerator generator(parms);
	generator.generate(4); // int parameter is the evaluation key count 
	cout << "... generator finished" << endl;
	Ciphertext public_key = generator.public_key();
	Plaintext secret_key = generator.secret_key();
	EvaluationKeys evaluation_keys = generator.evaluation_keys();
	cout << "... key generation complete: public key, secret key and evaluation key" << endl;

	//FractionalEncoder to encode
	FractionalEncoder encoder(parms.plain_modulus(), parms.poly_modulus(), 64, 32, 3);
	// Create the rest of the tools
	Encryptor encryptor(parms, public_key);
	Evaluator evaluator(parms, evaluation_keys);
	Decryptor decryptor(parms, secret_key);
	
	//// encrypt every number of a 2D array into cipertext 
	vector<vector<Ciphertext>> dataEncrypted = encrypt(data, encoder, encryptor); 

	vector<vector<Plaintext>> weightsEncoded = encodeWB(weights, encoder); // encode weights to do FHE plain_mul calculation
	vector<vector<Plaintext>> biasEncoded = encodeWB(bias, encoder); // encode bias to do FHE plain_plus calculation

	vector<vector<Plaintext>> weightsOptEncoded = encodeWB(weightsOutput, encoder); // encode weights to do FHE plain_mul calculation
	vector<vector<Plaintext>> biasOptEncoded = encodeWB(biasOutput, encoder); // encode bias to do FHE plain_plus calculation

	//FC no activation
	//vector<vector<Ciphertext>> fc1 = fullyConnect(dataEncrypted, weightsEncoded, biasEncoded, encoder, evaluator);
	//vector<vector<Ciphertext>> fcOutput = fullyConnect(fc1, weightsOptEncoded, biasOptEncoded, encoder, evaluator);
	
	//fc_BN_Relu 3rd poly
	//vector<vector<Ciphertext>> fcOutput = fc_BN_Relu(dataEncrypted, weightsEncoded, biasEncoded, encoder, evaluator);
	vector<vector<Ciphertext>> fc1 = fc_BN_Relu(dataEncrypted, weightsEncoded, biasEncoded, encoder, evaluator); // with relu
	vector<vector<Ciphertext>> fcOutput = fullyConnect(fc1, weightsOptEncoded, biasOptEncoded, encoder, evaluator); // without relu
	cout << "fcOutput rows" << to_string(fcOutput.size()) << endl;
	cout << "fcOutput cols" << to_string(fcOutput.at(0).size()) << endl;

	vector<vector<double>> decryptedResults = decypt2D(fcOutput, decryptor, encoder);

	// printing decryptedResults to check
	for (int i = 0; i < decryptedResults.size(); i++) {
		for (int j = 0; j < decryptedResults.at(0).size(); j++) {
			cout << to_string(decryptedResults.at(i).at(j)) << ((j < decryptedResults.at(i).size()) ? ", " : ".\n");
		}
	}

	return 1;
}

vector<vector<double>> decypt2D(vector<vector<Ciphertext>> encrypted, Decryptor decryptor, FractionalEncoder encoder) {
	int nb_rows = encrypted.size();
	int nb_cols = encrypted.at(0).size();

	vector<vector<double>> results;

	// Decrypt
	cout << "Decrypting ... ";
	for (int i = 0; i < nb_rows; i++) {
		vector<double> decryptedRow;
		for (int j = 0; j < nb_cols; j++) {
			Plaintext plain_result = decryptor.decrypt(encrypted.at(i).at(j)); // first decrypt
			double element_result = encoder.decode(plain_result);                      // then decode to print
			decryptedRow.emplace_back(element_result);
		}
		results.emplace_back(decryptedRow);
	}
	cout << " 2D decrypted done." << endl;
	return results;
}


vector<vector<Plaintext>> transpose(vector<vector<Plaintext>> w) {
	int w_rows = w.size(); 
	int w_cols = w.at(0).size();

	vector<vector<Plaintext>> transposed;

	for (int i = 0; i < w_cols; i++) {
		vector<Plaintext> transposedRow;
		for (int j = 0; j < w_rows; j++) {
			transposedRow.emplace_back(w.at(j).at(i));
		}
		transposed.emplace_back(transposedRow);
	}

	return transposed;
}

vector<vector<Ciphertext>> fullyConnect(vector<vector<Ciphertext>> X, vector<vector<Plaintext>> w, vector<vector<Plaintext>> b, FractionalEncoder encoder, Evaluator evaluator) {
	// fully connected layer
	vector<vector<Ciphertext>> fcOutput; // this is to be returned

	vector<vector<Plaintext>> w_transposed = transpose(w);
	vector<vector<Plaintext>> b_transposed = transpose(b);

	int X_rows = X.size();
	int X_cols = X.at(0).size();

	int w_rows = w.size(); // !!!! NOTE w_rows is number_features from last layer
	int w_cols = w.at(0).size(); // !!!! NOTE w_cols is number of neurons in THIS layer

	int b_rows = b.size(); // b.size() == 1
	int b_cols = b.at(0).size(); // b_cols == w_cols = number of neurons in this layer

	int w_rows_tr = w_transposed.size(); // 
	int w_cols_tr = w_transposed.at(0).size();

	int b_rows_tr = b_transposed.size(); // 
	int b_cols_tr = b_transposed.at(0).size();

	cout << "xrows: " << to_string(X_rows) << endl;
	cout << "xcols: " << to_string(X_cols) << endl;

	cout << "wrows: " << to_string(w_rows) << endl;
	cout << "wcols: " << to_string(w_cols) << endl;
	cout << "brows: " << to_string(b_rows) << endl;
	cout << "bcols: " << to_string(b_cols) << endl;

	cout << "wrows^: " << to_string(w_rows_tr) << endl;
	cout << "wcols^: " << to_string(w_cols_tr) << endl;
	cout << "brows^: " << to_string(b_rows_tr) << endl;
	cout << "bcols^: " << to_string(b_cols_tr) << endl;

	for (int i = 0; i < X_rows; i++) {
		cout << "Encrypting ... ";
		vector<Ciphertext> fc_row; // fc_row has length = w_cols = w_row_tr which is nb_neurons
		for (int j = 0; j < w_rows_tr; j++) {
			Ciphertext element = weightedSum(X.at(i), w_transposed.at(j), b_transposed.at(j), evaluator); // this is weighted sum of X row i and w col j
			fc_row.emplace_back(element);
		}
		fcOutput.emplace_back(fc_row); // fcOutpu has length = X_rows which is nb_samples
		cout << "finished calculating " << to_string(i) << " row" << endl;
	}
	return fcOutput;
}

Ciphertext weightedSum(vector<Ciphertext> x, vector<Plaintext> w, vector<Plaintext> b, Evaluator evaluator) {
	
	//Ciphertext enc_plain_product = evaluator.multiply_plain(encrypted_rationals[i], encoded_coefficients[i]);
	vector<Ciphertext> encrypted_products; // x.*w = x(1)*w(1) + x(2)*w(2) + ...
	for (int i = 0; i < w.size(); i++) {
		Ciphertext enc_plain_product = evaluator.multiply_plain(x.at(i), w.at(i)); // this is to cal x(i)*w(i), put them in a vector then add all elements in that vector to get the dot product;
		encrypted_products.emplace_back(enc_plain_product);
	}
	Ciphertext encrypted_dot_product = evaluator.add_many(encrypted_products);
	encrypted_dot_product = evaluator.add_plain(encrypted_dot_product, b.at(0)); // adding bias
	return encrypted_dot_product;
}

vector<vector<Plaintext>> encodeWB(vector<vector<double>> weights, FractionalEncoder encoder) {
	// encode 2D doubles into 2D Plaintext to do the FHE calculation
	vector<vector<Plaintext> > weightsEncoded; // this is to be returned 
	for (int i = 0; i < weights.size(); i++) {
		//loadValues.at(i); // encrypt row i
		//cout << "Encoding ... ";
		vector<Plaintext> encoded_rationals; // every element is encoded then emplace_back to this palintext vector which is then emplaced back to returned vet<vet>
		for (int j = 0; j < weights.at(i).size(); j++)
		{
			Plaintext encoded_number = encoder.encode(weights.at(i).at(j)); // encode the element A[i][j]
			encoded_rationals.emplace_back(encoded_number);
			//cout << to_string(weights.at(i).at(j)) << ((j < weights.at(i).size()) ? ", " : ".\n");
		}
		weightsEncoded.emplace_back(encoded_rationals);
		//cout << "finished encoding " << to_string(i) << " row" << endl;
	}
	
	return weightsEncoded;
}

vector<vector<Ciphertext>> fc_BN_Relu(vector<vector<Ciphertext>> X, vector<vector<Plaintext>> w, vector<vector<Plaintext>> b, FractionalEncoder encoder, Evaluator evaluator) {
	// fully connected layer + BN + Relu
	// input 2D vector has (nb_sample, nb_features)
	// input weights has (nb_features, nb_neurons)
	// input bias has (nb_neurons, 1)
	cout << "----------this layer is fc + bn + relu-------------------" << endl;

	double sigma_ = 0.1; // this is actually 1 divide sigma , just to make it easy, note that type has to be double, cause we used franctional encoder
	int reluDeg = 3;
	vector<vector<Ciphertext>> fcOutput; // this is to be returned

	vector<vector<Plaintext>> w_transposed = transpose(w);
	vector<vector<Plaintext>> b_transposed = transpose(b);

	int X_rows = X.size();
	int X_cols = X.at(0).size();

	int w_rows = w.size(); // !!!! NOTE w_rows is number_features from last layer
	int w_cols = w.at(0).size(); // !!!! NOTE w_cols is number of neurons in THIS layer

	int b_rows = b.size();
	int b_cols = b.at(0).size();

	int w_rows_tr = w_transposed.size(); // 
	int w_cols_tr = w_transposed.at(0).size();

	int b_rows_tr = b_transposed.size(); // 
	int b_cols_tr = b_transposed.at(0).size();

	cout << "xrows: " << to_string(X_rows) << endl;
	cout << "xcols: " << to_string(X_cols) << endl;

	cout << "wrows: " << to_string(w_rows) << endl;
	cout << "wcols: " << to_string(w_cols) << endl;
	cout << "brows: " << to_string(b_rows) << endl;
	cout << "bcols: " << to_string(b_cols) << endl;

	cout << "wrows^: " << to_string(w_rows_tr) << endl;
	cout << "wcols^: " << to_string(w_cols_tr) << endl;
	cout << "brows^: " << to_string(b_rows_tr) << endl;
	cout << "bcols^: " << to_string(b_cols_tr) << endl;

	for (int i = 0; i < X_rows; i++) {
		cout << "Encrypting ... ";
		vector<Ciphertext> fc_row; // fc_row has length = w_cols = w_row_tr which is nb_neurons
		for (int j = 0; j < w_rows_tr; j++) {
			Ciphertext element = weightedSum(X.at(i), w_transposed.at(j), b_transposed.at(j), evaluator); // this is weighted sum of X row i and w col j
			//batch normalization
			//element = evaluator.multiply_plain(element, encoder.encode(sigma_)); // BN
			//relu : 0.15+0.5012*x+0.2981*x*x-0.0004*x*x*x -(0.0388)*x*x*x*x
			element = relu(element, encoder, evaluator, reluDeg);
			fc_row.emplace_back(element);
		}
		fcOutput.emplace_back(fc_row); // fcOutpu has length = X_rows which is nb_samples
		cout << "finished calculating " << to_string(i) << " row" << endl;
	}
	return fcOutput;
}
Ciphertext relu(Ciphertext x, FractionalEncoder encoder, Evaluator evaluator, int deg) {
	if (deg == 4) {
		// relu : 0.15+0.5012*x+0.2981*x*x-0.0004*x*x*x -(0.0388)*x*x*x*x
		int destiSize = 2;

		Plaintext degree0 = encoder.encode(0.15);
		Ciphertext degree1 = evaluator.multiply_plain(x, encoder.encode(0.5002));
		Ciphertext degree01 = evaluator.add_plain(degree1, degree0);

		Ciphertext degree2 = evaluator.exponentiate(x, 2);
		degree2 = evaluator.multiply_plain(degree2, encoder.encode(0.2981));
		degree2 = evaluator.relinearize(degree2);

		Ciphertext degree3 = evaluator.exponentiate(x, 3);
		degree3 = evaluator.multiply_plain(degree3, encoder.encode(-0.0004));
		degree3 = evaluator.relinearize(degree3);

		Ciphertext degree4 = evaluator.exponentiate(x, 4);
		degree4 = evaluator.multiply_plain(degree4, encoder.encode(-0.0388));
		degree4 = evaluator.relinearize(degree4);

		vector<Ciphertext> encrypted_degrees;
		encrypted_degrees.emplace_back(degree01);
		encrypted_degrees.emplace_back(degree2);
		encrypted_degrees.emplace_back(degree3);
		encrypted_degrees.emplace_back(degree4);

		Ciphertext result = evaluator.add_many(encrypted_degrees);
		cout << "ciphertext size" << to_string(result.size()) << endl;
		return result;
	}else if (deg == 3) {
		//// degree3 relu : 0.1995 + 0.5002x + 0.1994x^2 -0.0164X^3
		int destiSize = 2;

		Plaintext degree0 = encoder.encode(0.1995);
		Ciphertext degree1 = evaluator.multiply_plain(x, encoder.encode(0.5002));
		Ciphertext degree01 = evaluator.add_plain(degree1, degree0);

		Ciphertext degree2 = evaluator.exponentiate(x, 2);
		degree2 = evaluator.multiply_plain(degree2, encoder.encode(0.1994));
		degree2 = evaluator.relinearize(degree2);

		Ciphertext degree3 = evaluator.exponentiate(x, 3);
		degree3 = evaluator.multiply_plain(degree3, encoder.encode(-0.0164));
		degree3 = evaluator.relinearize(degree3);
		cout << "degree 3 finished -----------------" << endl;

		vector<Ciphertext> encrypted_degrees;
		encrypted_degrees.emplace_back(degree01);
		encrypted_degrees.emplace_back(degree2);
		encrypted_degrees.emplace_back(degree3);

		Ciphertext result = evaluator.add_many(encrypted_degrees);
		cout << "ciphertext size" << to_string(result.size()) << endl;
		return result;
	}else if (deg == 2) {
		//// degree2 relu : 0.1992 + 0.5002x + 0.1997x^2;
		int destiSize = 2;

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
		cout << "ciphertext size" << to_string(result.size()) << endl;
		return result;
	}

}