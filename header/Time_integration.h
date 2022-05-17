#pragma once

#include <vector>
#include "Spatial_discrete.h"

class Timeintegration {
public:

	Timeintegration() {};
	~Timeintegration() {};

	std::shared_ptr<Spatialdiscrete> Spatialdiscrete_;

	void initialize();

	std::vector<double> Residual;

	int test;
	
	int num_center_x;

	int num_center_y;

	int num_total_center;

	int num_total_var;

	virtual void Update(std::vector<double>& Solution, double dt) =0;

	double Cal_dt(std::vector<double>& Solution);

	std::vector<double> intermediate_1;

	std::vector<double> intermediate_2;

	std::vector<double> intermediate_3;

	std::vector<double> intermediate_4;

	std::vector<double> Cell_area;

};

class SSPRK : public Timeintegration {

	void Update(std::vector<double>& Solution, double dt); // &로 받아와서 업데이트된 값 저장해주기

};

class RK3 : public Timeintegration {

	void Update(std::vector<double>& Solution, double dt); // &로 받아와서 업데이트된 값 저장해주기

};