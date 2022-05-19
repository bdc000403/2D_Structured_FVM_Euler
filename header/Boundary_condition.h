#pragma once

#include <vector>
#include <string>
#include "Config.h"
#include "Grid.h"

class Boundarycondition {

public:
	Boundarycondition() {};
	~Boundarycondition() {};

	int num_grid_point_x;

	int num_grid_point_y;

	int num_total_grid_point;

	int num_center_x;

	int num_center_y;

	int num_total_center_point;

	int num_ghost_x;

	int num_ghost_y;

	double Mach_number;

	double AOA;

	std::vector<int> BC;

	std::vector<double> nx;

	std::vector<double> ny;

	void Cell_boundary_condition(std::vector<double>& Ghost_solution);

	void Flux_boundary_condition(std::vector<double>& Solution, std::vector<double>& Flux);

	void(Boundarycondition::* BC_Imin) (std::vector<double>& Ghost_solution, int loc, int start_index, int end_index);

	void(Boundarycondition::* BC_Imax) (std::vector<double>& Ghost_solution, int loc, int start_index, int end_index);

	void(Boundarycondition::* BC_Jmin) (std::vector<double>& Ghost_solution, int loc, int start_index, int end_index);

	void(Boundarycondition::* BC_Jmax) (std::vector<double>& Ghost_solution, int loc, int start_index, int end_index);

	void(Boundarycondition::* Flux_BC_Imin) (std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

	void(Boundarycondition::* Flux_BC_Imax) (std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

	void(Boundarycondition::* Flux_BC_Jmin) (std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

	void(Boundarycondition::* Flux_BC_Jmax) (std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

	void BC_Imin_ptr(std::vector<double>& Ghost_solution, int start_index, int end_index) { (this->*BC_Imin)(Ghost_solution, 0, start_index, end_index); };

	void BC_Imax_ptr(std::vector<double>& Ghost_solution, int start_index, int end_index) { (this->*BC_Imax)(Ghost_solution, 1, start_index, end_index); };

	void BC_Jmin_ptr(std::vector<double>& Ghost_solution, int start_index, int end_index) { (this->*BC_Jmin)(Ghost_solution, 2, start_index, end_index); };

	void BC_Jmax_ptr(std::vector<double>& Ghost_solution, int start_index, int end_index) { (this->*BC_Jmax)(Ghost_solution, 3, start_index, end_index); };

	void Flux_BC_Imin_ptr(std::vector<double>& Solution, std::vector<double>& Flux, int start_index, int end_index) { (this->*Flux_BC_Imin)(Solution, Flux, 0, start_index, end_index); };

	void Flux_BC_Imax_ptr(std::vector<double>& Solution, std::vector<double>& Flux, int start_index, int end_index) { (this->*Flux_BC_Imax)(Solution, Flux, 1, start_index, end_index); };

	void Flux_BC_Jmin_ptr(std::vector<double>& Solution, std::vector<double>& Flux, int start_index, int end_index) { (this->*Flux_BC_Jmin)(Solution, Flux, 2, start_index, end_index); };

	void Flux_BC_Jmax_ptr(std::vector<double>& Solution, std::vector<double>& Flux, int start_index, int end_index) { (this->*Flux_BC_Jmax)(Solution, Flux, 3, start_index, end_index); };

	void BC_initialize();

	void Far_field(std::vector<double>& Ghost_olution, int loc, int start_index, int end_index);

	void Periodic(std::vector<double>& Ghost_olution, int loc, int start_index, int end_index);

	void Inviscid_wall(std::vector<double>& Ghost_olution, int loc, int start_index, int end_index);
	
	void Flux_Inviscid_wall(std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

	void Flux_NON(std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index);

};