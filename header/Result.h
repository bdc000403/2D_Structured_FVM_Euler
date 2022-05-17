#pragma once

#include "Config.h"
#include "Grid.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

class Result {

public:
	Result() {};
	~Result () {};

	int num_center_x;

	int num_center_y;

	int num_var;

	double Error;

	double Cl;			// Cal_cl계산부터

	double AOA;

	double Free_mach;

	double Cd;

	double Fx;

	double Fy;

	double Lift;

	double Drag;

	std::vector<double> Prev_solution;

	std::vector<double> Surface_length;

	std::vector<double> nx;

	std::vector<double> ny;

	std::string file_name;

	std::vector<double> x;

	std::vector<double> y;

	std::string conv_file_name;

	std::ofstream conv_file;

	std::string cp_file_name;

	void Save(std::vector<double>& Solution, double time, int file_index);

	void Save(std::vector<double>& Solution, double time, std::string file_name);

	void Post(std::vector<double>& Solution, double time, std::string file_name);

	void Save_conv_hist(std::vector<double>& Solution, int step, double time, double Error);

	void Save_cp(std::vector<double>& Solution);

	double Cal_Error(std::vector<double>& Solution);

	void Cal_coefficient(std::vector<double>& Soution);

	void Close_file();

	void Test_save(std::vector<double>& Solution, double time, int file_index);

	void Initialize();

};

