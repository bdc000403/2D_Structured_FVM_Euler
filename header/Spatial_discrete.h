#pragma once

#include <vector>
#include <string>
#include "Grid.h"
#include "Config.h"
#include "Flux.h"
#include "Boundary_condition.h"

class Spatialdiscrete {

public:
	Spatialdiscrete(){}
	~Spatialdiscrete(){}

	std::shared_ptr<FVS> FVS_;
	std::shared_ptr<Boundarycondition> BC_;

	void initialize();			// 함수 종류 저장

	void Cal_Residual(std::vector<double>& Solution, std::vector<double>& Residual);

	void (Spatialdiscrete::*Cal_interpolation_) (std::vector<double>& Solution);

	void Cal_interpolation(std::vector<double>& Solution);

	void Cal_coordinate_transform();

	void Cal_inverse_coordinate_transform();

	void Cal_smoothness();

	double Cal_dt(std::vector<double>& Solution);

	virtual void Cal_phi() = 0;

	std::vector<int> BC;

	std::vector<double> Flux;

private:

	double Mach_number;

	double AOA;


protected:

	int num_grid_point_x;

	int num_grid_point_y;

	int num_total_grid_point;

	int num_center_x;

	int num_center_y;

	int num_total_center_point;

	int num_ghost_x;

	int num_ghost_y;

	double CFL;

	std::vector<double> x_phi;

	std::vector<double> y_phi;

	std::vector<double> x_theta;

	std::vector<double> y_theta;

	std::vector<double> nx;

	std::vector<double> Cell_area;

	std::vector<double> ny;

	std::vector<double> North_solution;		//jmax boundary

	std::vector<double> South_solution;

	std::vector<double> West_solution;

	std::vector<double> East_solution;

	std::vector<double> Ghost_solution;

	std::vector<double> Surface_length;

	std::vector<double> dx;

	std::vector<double> dy;

};


class Minmod_limiter : public Spatialdiscrete {
	


public:

	Minmod_limiter() {};
	~Minmod_limiter() {};

	void Cal_phi();

};

class Off : public Spatialdiscrete {



public:

	Off() {};
	~Off() {};

	void Cal_phi();

};