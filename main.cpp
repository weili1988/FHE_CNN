/*
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "seal.h"
#include "bigpoly.h"
*/

#include "tests.h"

int main()
{
	// cnn test 8 nodes in 1st layer and 16 nodes in 2nd layer plus 1 FC
	//string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights_test";
	//CNN_test(path);

	// cnn test slim version 1 node on each CNN layer plus 1 FC
	// this is path where weights and test data are stored, this path shall have cnn1W cnn1b; cnn2W, cnn2b, fcW and fcB together with. digit.csv
	string path = "C:/Users/liw5/Downloads/self_learn/infosec/encryption_DL/FHE_MachineLearning/data/MNIST_CNN_weights_test"; 
	CNN_test_slim(path);

    // Wait for ENTER before closing screen.
    cout << "Press ENTER to exit" << endl;
    char ignore;
    cin.get(ignore);

    //return result;
	return -1;
}
