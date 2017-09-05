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
	//CNN_test();

	// cnn test slim version 1 node on each CNN layer plus 1 FC
	CNN_test_slim();

    // Wait for ENTER before closing screen.
    cout << "Press ENTER to exit" << endl;
    char ignore;
    cin.get(ignore);

    //return result;
	return -1;
}
