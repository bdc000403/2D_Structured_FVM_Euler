#pragma once

#include "Grid.h"
#include "Config.h"
#include "Boundary_condition.h"
#include <vector>
#include <string>



class FVS {

public:
	std::vector<double> uL_;

	std::vector<double> uR_;

	std::vector<double> vL_;

	std::vector<double> vR_;

	std::vector<double> pL_;

	std::vector<double> pR_;

	std::vector<double> aR_;

	std::vector<double> aL_;

	std::vector<double> HR_;

	std::vector<double> HL_;

	std::vector<double> mach_P;

	std::vector<double> mach_M;

	std::vector<double> mach_L;

	std::vector<double> mach_R;

	std::vector<double> mach_interface;

	std::vector<double> fluxP_;

	std::vector<double> fluxM_;

	int num_grid_point_x;

	int num_grid_point_y;

	int num_total_grid_point;

	int num_center_x;

	int num_center_y;

	int num_total_center_point;

	int num_ghost_x;

	int num_ghost_y;

	FVS() {};
	~FVS() {};

	std::shared_ptr<Boundarycondition> BC_;

	double max_characteristic_speed;

	void initialize();

	void Cal_Flux(std::vector<double>& South_solution, std::vector<double>& East_solution, std::vector<double>& North_solution, std::vector<double>& West_solution, std::vector<double>& Flux);

	virtual void Cal_Flux_split(std::vector<double>& Left_solutuon, std::vector<double>& Right_solution, std::vector<double>& temp_flux, int total_size) = 0;


};

class VanLeer : public FVS {

public:

	VanLeer() {};
	~VanLeer() {};

	void Cal_Flux_split(std::vector<double>& Left_solutuon, std::vector<double>& Right_solution, std::vector<double>& temp_flux, int total_size);

};