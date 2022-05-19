#include "../header/Boundary_condition.h"



void Boundarycondition::BC_initialize() {

	Grid& grid = Grid::get_instance();
	BC = grid.get_BC();

	Config& config = Config::get_instance();
	Mach_number = config.get_Mach_number();
	AOA = config.get_AOA();
	AOA = AOA * 3.141592653589 / 180;

	num_grid_point_x = grid.get_num_grid_point_x();
	num_grid_point_y = grid.get_num_grid_point_y();
	num_total_grid_point = grid.get_num_total_grid_point();

	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center_point = grid.get_num_total_center_point();

	num_ghost_x = num_center_x + 2;
	num_ghost_y = num_center_y + 2;

	nx = grid.get_nx();
	ny = grid.get_ny();

	int  Imin_type = 0;
	int  Imax_type = 0;
	int  Jmin_type = 0;
	int  Jmax_type = 0;

	Imin_type = BC[0 * 4 + 2];
	Imax_type = BC[1 * 4 + 2];
	Jmin_type = BC[2 * 4 + 2];
	Jmax_type = BC[3 * 4 + 2];

	// Imin 
	if (Imin_type == 14) {
		BC_Imin = &Boundarycondition::Far_field;
		Flux_BC_Imin = &Boundarycondition::Flux_NON;
	}
	else if (Imin_type == 13) {
		BC_Imin = &Boundarycondition::Periodic;
		Flux_BC_Imin = &Boundarycondition::Flux_NON;
	}
	else if (Imin_type == 3) {
		BC_Imin = &Boundarycondition::Inviscid_wall;
		Flux_BC_Imin = &Boundarycondition::Flux_Inviscid_wall;
	}

	// Imax
	if (Imax_type == 14) {
		BC_Imax = &Boundarycondition::Far_field;
		Flux_BC_Imax = &Boundarycondition::Flux_NON;
	}
	else if (Imax_type == 13) {
		BC_Imax = &Boundarycondition::Periodic;
		Flux_BC_Imax = &Boundarycondition::Flux_NON;
	}
	else if (Imax_type == 3) {
		BC_Imax = &Boundarycondition::Inviscid_wall;
		Flux_BC_Imax = &Boundarycondition::Flux_Inviscid_wall;
	}

	// Jmin
	if (Jmin_type == 14) {
		BC_Jmin = &Boundarycondition::Far_field;
		Flux_BC_Jmin = &Boundarycondition::Flux_NON;
	}
	else if (Jmin_type == 13) {
		BC_Jmin = &Boundarycondition::Periodic;
		Flux_BC_Jmin = &Boundarycondition::Flux_NON;
	}
	else if (Jmin_type == 3) {
		BC_Jmin = &Boundarycondition::Inviscid_wall;
		Flux_BC_Jmin = &Boundarycondition::Flux_Inviscid_wall;
	}

	// Jmax
	if (Jmax_type == 14) {
		BC_Jmax = &Boundarycondition::Far_field;
		Flux_BC_Jmax = &Boundarycondition::Flux_NON;
	}
	else if (Jmax_type == 13) {
		BC_Jmax = &Boundarycondition::Periodic;
		Flux_BC_Jmax = &Boundarycondition::Flux_NON;
	}
	else if (Jmax_type == 3) {
		BC_Jmax = &Boundarycondition::Inviscid_wall;
		Flux_BC_Jmax = &Boundarycondition::Flux_Inviscid_wall;
	}

}

void Boundarycondition::Cell_boundary_condition(std::vector<double>& Ghost_solution) {

	int start_index = 0;
	int end_index = 0;

	for (int loc = 0; loc < 4; loc++) {

		// loc : 0 = imin, 1 = imax, 2 = jmin, 3 = jmax

		start_index = BC[4 * loc + 0];
		end_index = BC[4 * loc + 1];

		if (loc == 0) {
			BC_Imin_ptr(Ghost_solution, start_index, end_index);
		}
		else if (loc == 1) {
			BC_Imax_ptr(Ghost_solution, start_index, end_index);
		}
		else if (loc == 2) {
			BC_Jmin_ptr(Ghost_solution, start_index, end_index);
		}
		else if (loc == 3) {
			BC_Jmax_ptr(Ghost_solution, start_index, end_index);
		}

	}

}

void Boundarycondition::Flux_boundary_condition(std::vector<double>& Solution, std::vector<double>& Flux) {

	int ghost_index = 0;
	int cell_index = 0;

	int start_index = 0;
	int end_index = 0;
	int BC_type = 0;

	for (int loc = 0; loc < 4; loc++) {

		// loc : 0 = imin, 1 = imax, 2 = jmin, 3 = jmax

		start_index = BC[4 * loc + 0];
		end_index = BC[4 * loc + 1];

		if (loc == 0) {
			Flux_BC_Imin_ptr(Solution, Flux, start_index, end_index);
		}
		else if (loc == 1) {
			Flux_BC_Imax_ptr(Solution, Flux, start_index, end_index);
		}
		else if (loc == 2) {
			Flux_BC_Jmin_ptr(Solution, Flux, start_index, end_index);
		}
		else if (loc == 3) {
			Flux_BC_Jmax_ptr(Solution, Flux, start_index, end_index);
		}

	}

}


void Boundarycondition::Far_field(std::vector<double>& Ghost_solution, int loc, int start_index, int end_index) {

	int ghost_index = 0;

	double gamma = 1.4;
	double rho_inf = 1;
	double sos_inf = 1;
	double p_inf = rho_inf * pow(sos_inf, 2) / gamma;

	double x_mom = 1 * Mach_number * cos(AOA);
	double y_mom = 1 * Mach_number * sin(AOA);
	double totalE_per_volume = rho_inf * 0.5 * pow(Mach_number, 2) + p_inf / ((gamma - 1) * rho_inf);


	if (loc == 2) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * (icell + 1);

			Ghost_solution[ghost_index + 0] = rho_inf;
			Ghost_solution[ghost_index + 1] = x_mom;
			Ghost_solution[ghost_index + 2] = y_mom;
			Ghost_solution[ghost_index + 3] = totalE_per_volume;

		}
	}
	else if (loc == 3) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * ((icell + 1) + (num_center_y + 1) * num_ghost_x);

			Ghost_solution[ghost_index + 0] = rho_inf;
			Ghost_solution[ghost_index + 1] = x_mom;
			Ghost_solution[ghost_index + 2] = y_mom;
			Ghost_solution[ghost_index + 3] = totalE_per_volume;

		}
	}
	else if (loc == 0) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * (0 + num_ghost_x * (jcell + 1));

			Ghost_solution[ghost_index + 0] = rho_inf;
			Ghost_solution[ghost_index + 1] = x_mom;
			Ghost_solution[ghost_index + 2] = y_mom;
			Ghost_solution[ghost_index + 3] = totalE_per_volume;

		}
	}
	else if (loc == 1) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * ((num_center_x + 1) + num_ghost_x * (jcell + 1));

			Ghost_solution[ghost_index + 0] = rho_inf;
			Ghost_solution[ghost_index + 1] = x_mom;
			Ghost_solution[ghost_index + 2] = y_mom;
			Ghost_solution[ghost_index + 3] = totalE_per_volume;

		}
	}
}

void Boundarycondition::Periodic(std::vector<double>& Ghost_solution, int loc, int start_index, int end_index) {

	int ghost_index = 0;
	int period_index = 0;

	if (loc == 2) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * (icell + 1);
			period_index = 4 * ((icell + 1) + (num_center_y)*num_ghost_x);

			Ghost_solution[ghost_index + 0] = Ghost_solution[period_index + 0];
			Ghost_solution[ghost_index + 1] = Ghost_solution[period_index + 1];
			Ghost_solution[ghost_index + 2] = Ghost_solution[period_index + 2];
			Ghost_solution[ghost_index + 3] = Ghost_solution[period_index + 3];

		}
	}
	else if (loc == 3) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * ((icell + 1) + (num_center_y + 1) * num_ghost_x);
			period_index = 4 * ((icell + 1) + num_ghost_x);

			Ghost_solution[ghost_index + 0] = Ghost_solution[period_index + 0];
			Ghost_solution[ghost_index + 1] = Ghost_solution[period_index + 1];
			Ghost_solution[ghost_index + 2] = Ghost_solution[period_index + 2];
			Ghost_solution[ghost_index + 3] = Ghost_solution[period_index + 3];

		}
	}
	else if (loc == 0) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * (num_ghost_x * (jcell + 1));
			period_index = 4 * (num_ghost_x * (jcell + 1) + num_center_x);

			Ghost_solution[ghost_index + 0] = Ghost_solution[period_index + 0];
			Ghost_solution[ghost_index + 1] = Ghost_solution[period_index + 1];
			Ghost_solution[ghost_index + 2] = Ghost_solution[period_index + 2];
			Ghost_solution[ghost_index + 3] = Ghost_solution[period_index + 3];

		}
	}
	else if (loc == 1) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * ((num_center_x + 1) + num_ghost_x * (jcell + 1));
			period_index = 4 * (1 + num_ghost_x * (jcell + 1));

			Ghost_solution[ghost_index + 0] = Ghost_solution[period_index + 0];
			Ghost_solution[ghost_index + 1] = Ghost_solution[period_index + 1];
			Ghost_solution[ghost_index + 2] = Ghost_solution[period_index + 2];
			Ghost_solution[ghost_index + 3] = Ghost_solution[period_index + 3];

		}
	}




}

void Boundarycondition::Inviscid_wall(std::vector<double>& Ghost_solution, int loc, int start_index, int end_index) {
	int ghost_index = 0;
	int w2_index = 0;				//Blazek p.276 참고
	int w3_index = 0;

	if (loc == 2) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * (icell + 1);
			w2_index = 4 * (icell + 1 + (num_ghost_x));
			w3_index = 4 * (2*(num_ghost_x) + (icell + 1));

			Ghost_solution[ghost_index + 0] = 2 * Ghost_solution[w2_index + 0] - Ghost_solution[w3_index + 0];
			Ghost_solution[ghost_index + 1] = 2 * Ghost_solution[w2_index + 1] - Ghost_solution[w3_index + 1];
			Ghost_solution[ghost_index + 2] = 2 * Ghost_solution[w2_index + 2] - Ghost_solution[w3_index + 2];
			Ghost_solution[ghost_index + 3] = 2 * Ghost_solution[w2_index + 3] - Ghost_solution[w3_index + 3];

			//Ghost_solution[ghost_index + 0] = Ghost_solution[w2_index + 0];
			//Ghost_solution[ghost_index + 1] = Ghost_solution[w2_index + 1];
			//Ghost_solution[ghost_index + 2] = Ghost_solution[w2_index + 2];
			//Ghost_solution[ghost_index + 3] = Ghost_solution[w2_index + 3];

		}
	}
	else if (loc == 3) {
		for (int icell = 0; icell < num_center_x; icell++) {

			ghost_index = 4 * ((icell + 1) + (num_center_y + 1) * num_ghost_x);
			w2_index = ghost_index - 4 * num_ghost_x;
			w3_index = ghost_index - 4 * 2 * num_ghost_x;

			Ghost_solution[ghost_index + 0] = 2 * Ghost_solution[w2_index + 0] - Ghost_solution[w3_index + 0];
			Ghost_solution[ghost_index + 1] = 2 * Ghost_solution[w2_index + 1] - Ghost_solution[w3_index + 1];
			Ghost_solution[ghost_index + 2] = 2 * Ghost_solution[w2_index + 2] - Ghost_solution[w3_index + 2];
			Ghost_solution[ghost_index + 3] = 2 * Ghost_solution[w2_index + 3] - Ghost_solution[w3_index + 3];

		}
	}
	else if (loc == 0) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * (0 + num_ghost_x * (jcell + 1));
			w2_index = ghost_index + 4;
			w3_index = ghost_index + 4 * 2;

			Ghost_solution[ghost_index + 0] = 2 * Ghost_solution[w2_index + 0] - Ghost_solution[w3_index + 0];
			Ghost_solution[ghost_index + 1] = 2 * Ghost_solution[w2_index + 1] - Ghost_solution[w3_index + 1];
			Ghost_solution[ghost_index + 2] = 2 * Ghost_solution[w2_index + 2] - Ghost_solution[w3_index + 2];
			Ghost_solution[ghost_index + 3] = 2 * Ghost_solution[w2_index + 3] - Ghost_solution[w3_index + 3];

		}
	}
	else if (loc == 1) {
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			ghost_index = 4 * ((num_center_x + 1) + num_ghost_x * (jcell + 1));
			w2_index = ghost_index - 4;
			w3_index = ghost_index - 4 * 2;

			Ghost_solution[ghost_index + 0] = 2 * Ghost_solution[w2_index + 0] - Ghost_solution[w3_index + 0];
			Ghost_solution[ghost_index + 1] = 2 * Ghost_solution[w2_index + 1] - Ghost_solution[w3_index + 1];
			Ghost_solution[ghost_index + 2] = 2 * Ghost_solution[w2_index + 2] - Ghost_solution[w3_index + 2];
			Ghost_solution[ghost_index + 3] = 2 * Ghost_solution[w2_index + 3] - Ghost_solution[w3_index + 3];

		}
	}

}


void Boundarycondition::Flux_Inviscid_wall(std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index) {

	int cell_index = 0;
	int face_index = 0;
	int w2_index = 0;				//Blazek p.276 참고
	int w3_index = 0;
	double p2 = 0;
	double p3 = 0;
	double p_wall = 0;

	double r_ = 1.4;

	if (loc == 0) {												//Imin
		for (int jcell = 0; jcell < num_center_y; jcell++) {
	
				cell_index = 4 * (0 + num_center_x * jcell);
				w2_index = 4 * (1 + num_center_x * jcell);
				w3_index = 4 * (2 + num_center_x * jcell);
				p2 = (r_ - 1) * (Solution[w2_index + 3] - 0.5 * (pow(Solution[w2_index + 1], 2) + pow(Solution[w2_index + 2], 2)) / Solution[w2_index + 0]);
				p3 = (r_ - 1) * (Solution[w3_index + 3] - 0.5 * (pow(Solution[w3_index + 1], 2) + pow(Solution[w3_index + 2], 2)) / Solution[w3_index + 0]);
				p_wall = 0.5 * (3 * p2 - p3);
				face_index = 4 * (3 + 4 * (0 + num_center_x * jcell));			//west index
				Flux[face_index + 0] = 0;
				Flux[face_index + 1] = p_wall * nx[face_index];
				Flux[face_index + 2] = p_wall * ny[face_index];
				Flux[face_index + 3] = 0;

		}
	}
	else if (loc == 1) {										//Imax 
		for (int jcell = 0; jcell < num_center_y; jcell++) {

			cell_index = 4 * (num_center_x - 1 + num_center_x * jcell);
			w2_index = 4 * (num_center_x - 2 + num_center_x * jcell);
			w3_index = 4 * (num_center_x - 3 + num_center_x * jcell);
			p2 = (r_ - 1) * (Solution[w2_index + 3] - 0.5 * (pow(Solution[w2_index + 1], 2) + pow(Solution[w2_index + 2], 2)) / Solution[w2_index + 0]);
			p3 = (r_ - 1) * (Solution[w3_index + 3] - 0.5 * (pow(Solution[w3_index + 1], 2) + pow(Solution[w3_index + 2], 2)) / Solution[w3_index + 0]);
			p_wall = 0.5 * (3 * p2 - p3);
			face_index = 4 * (1 + 4 * (num_center_x - 1 + num_center_x * jcell));			//east index
			Flux[face_index + 0] = 0;
			Flux[face_index + 1] = p_wall * nx[face_index];
			Flux[face_index + 2] = p_wall * ny[face_index];
			Flux[face_index + 3] = 0;

		}
	}
	else if (loc == 2) {										//Jmin  face수정 필요
		for (int icell = 0; icell < num_center_x; icell++) {

			w2_index = 4 * (icell);
			w3_index = 4 * (icell + num_center_x);
			p2 = (r_ - 1) * (Solution[w2_index + 3] - 0.5 * (pow(Solution[w2_index + 1], 2) + pow(Solution[w2_index + 2], 2)) / Solution[w2_index + 0]);
			p3 = (r_ - 1) * (Solution[w3_index + 3] - 0.5 * (pow(Solution[w3_index + 1], 2) + pow(Solution[w3_index + 2], 2)) / Solution[w3_index + 0]);
			p_wall =  -0.5 * (3 * p2 - p3);
			face_index = 0 + 4 * (icell);			//south index
			Flux[face_index * 4 + 0] = 0;
			Flux[face_index * 4 + 1] = p_wall * nx[face_index];
			Flux[face_index * 4 + 2] = p_wall * ny[face_index];
			Flux[face_index * 4 + 3] = 0;

		}
	}
	else if (loc == 3) {										//Jmax
		for (int icell = 0; icell < num_center_x; icell++) {

			cell_index = 4 * (icell + num_center_x * (num_center_y - 1));
			w2_index = 4 * (icell + num_center_x * (num_center_y - 2));
			w3_index = 4 * (icell + num_center_x * (num_center_y - 3));
			p2 = (r_ - 1) * (Solution[w2_index + 3] - 0.5 * (pow(Solution[w2_index + 1], 2) + pow(Solution[w2_index + 2], 2)) / Solution[w2_index + 0]);
			p3 = (r_ - 1) * (Solution[w3_index + 3] - 0.5 * (pow(Solution[w3_index + 1], 2) + pow(Solution[w3_index + 2], 2)) / Solution[w3_index + 0]);
			p_wall = 0.5 * (3 * p2 - p3);
			face_index = 4 * (2 + 4 * (icell + num_center_x * (num_center_y - 1)));			//north index
			Flux[face_index + 0] = 0;
			Flux[face_index + 1] = p_wall * nx[face_index];
			Flux[face_index + 2] = p_wall * ny[face_index];
			Flux[face_index + 3] = 0;

		}
	}


}

void Boundarycondition::Flux_NON(std::vector<double>& Solution, std::vector<double>& Flux, int loc, int start_index, int end_index) {

	// NON

}