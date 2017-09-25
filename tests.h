#pragma once

#include <fstream>
#include <string>
#include "seal.h"
#include <cmath>
#include <chrono>
#include "CipherSaveLoad.h"
#include "AveragePooling.h"
#include "monitors.h"
#include "utility.h"
#include "SEAL_CNN.h"

using namespace std;
using namespace seal;

void CNN_test(string path);
void CNN_test_slim(string path, int digit);
int logitstic_regression();
int testfullyConnect();

void save_load_test(string path, string source);
void average_pooling_test(string path);
void monitor_test(string path, string source);
// test the working range for fractional encoder/decoder
bool ifDecodeCorrect(double number, FractionalEncoder encoder, double residue);
void encode_test(string path);
void timing_test(string path);

void print_example_banner(string title);
void example_basics();
void example_weighted_average();
void example_parameter_selection();
void example_batching();
void example_relinearization();
void example_relinearization_part1();
void example_relinearization_part2();
void example_relinearization_part3();
void example_timing();
void example_linear_regression();
