#pragma once

#include "Config.h"
#include "Grid.h"

class Initial {

public:

	Initial() {};
	~Initial() {};

private:

	double AOA;

	int num_grid_point_x;

	int num_grid_point_y;

	int num_total_grid_point;

	int num_center_x;

	int num_center_y;

	int num_total_center_point;

	double Mach_number;

public:

	void Free_stream(std::vector<double>& Solution);

	void Test(std::vector<double>& Solution);

};