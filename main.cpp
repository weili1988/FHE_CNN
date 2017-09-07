/*
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "seal.h"
#include "bigpoly.h"
*/
#include <iostream>
#include <string>
#include <sstream>
#include "tests.h"

int main()
{
	// cnn test 8 nodes in 1st layer and 16 nodes in 2nd layer plus 1 FC
	//string path = "./CNN_cpp/MNIST_CNN_weights";
	//CNN_test(path);

	// cnn test slim version 1 node on each CNN layer plus 1 FC
	// this is path where weights and test data are stored, this path shall have cnn1W cnn1b; cnn2W, cnn2b, fcW and fcB together with. digit.csv
	//string path = "./CNN_cpp/MNIST_CNN_weights"; 
	//cout << "assuming you have been using default folder, please enter the digit you want to predict, hit enter when done" << endl;
	//int digit; // choose a digit from 0 to 9
	//cin >> digit;
	//cout << "predicting digit:" << to_string(digit) << endl;

	//CNN_test_slim(path, digit);
	// path to store the txt files
	string path = "./CNN_cpp/test";
	encode_test(path);

    //// Wait for ENTER before closing screen.
    //cout << "Press ENTER to exit" << endl;
    //char ignore;
    //cin.get(ignore);

    //return result;
	return -1;
}
