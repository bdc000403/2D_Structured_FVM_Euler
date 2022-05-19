#pragma once

#include "Config.h"
#include "Grid.h"
#include "Result.h"
#include "Initial_condition.h"
#include "Time_integration.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

class Solver {

public:

	Solver() {};
	~Solver() {};

	std::shared_ptr<Timeintegration> Timeintegration_;

	void initialize();

	void Solve_Unsteady();

	void Solve_Steady();

	double End_time;

	double CFL;

	double Tolerance;

	double Time_interval;

	int num_center_x;

	int num_center_y;

	int num_result_file;

	int Max_iter;

	int Save_interval;

	std::vector<double> Solution;

	std::vector<double> Residual;

	int dt;

	std::shared_ptr<Initial> Initial_;

	std::shared_ptr<Result> Result_;
	
	void (Timeintegration::* Update) (int dt); // Spatialdisrete에서 직접 residual 받아서 업데이트

private:



};