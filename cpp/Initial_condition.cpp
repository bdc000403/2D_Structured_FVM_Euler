#include "../header/Initial_condition.h"




void Initial::Free_stream(std::vector<double>& Solution) {


	Config& config = Config::get_instance();
	AOA = config.get_AOA();
	Mach_number = config.get_Mach_number();
	AOA = AOA * 3.141592653589 / 180;

	Grid& grid = Grid::get_instance();

	num_grid_point_x = grid.get_num_grid_point_x();
	num_grid_point_y = grid.get_num_grid_point_y();
	num_total_grid_point = grid.get_num_total_grid_point();
	std::vector<double> center_coord_x = grid.get_Cell_center_coordinates_x();
	std::vector<double> center_coord_y = grid.get_Cell_center_coordinates_y();
	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center_point = grid.get_num_total_center_point();
	double gamma = 1.4;
	double rho_inf = 1;
	double sos_inf = 1;
	double p_inf = rho_inf * pow(sos_inf, 2) / gamma;
	double x_mom = 1 * Mach_number * cos(AOA);
	double y_mom = 1 * Mach_number * sin(AOA);
	double totalE_per_volume = rho_inf * 0.5 * pow(Mach_number, 2) + p_inf / (gamma - 1);

	int cell_loc = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			cell_loc = 4 * (icell  + num_center_x * jcell);

			//rho_inf = - (pow(center_coord_x[cell_loc * 0.25] - 0.5, 2) + pow(center_coord_y[cell_loc * 0.25] - 0.5, 2)) + 5;
			//totalE_per_volume = rho_inf * 0.5 * pow(Mach_number, 2) + p_inf / (gamma - 1);

			Solution[cell_loc + 0] = rho_inf;
			Solution[cell_loc + 1] = x_mom;
			Solution[cell_loc + 2] = y_mom;
			Solution[cell_loc + 3] = totalE_per_volume;
			
		}
	}



}

void Initial::Test(std::vector<double>& Solution) {

	Grid& grid = Grid::get_instance();
	Config& config = Config::get_instance();

	num_grid_point_x = grid.get_num_grid_point_x();
	num_grid_point_y = grid.get_num_grid_point_y();
	num_total_grid_point = grid.get_num_total_grid_point();
	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center_point = grid.get_num_total_center_point();
	AOA = config.get_AOA();
	Mach_number = config.get_Mach_number();

	AOA = AOA * 3.141592653589 / 180;
	double gamma = 1.4;
	double rho_inf = 1;
	double sos_inf = 1;
	double p_inf = rho_inf * pow(sos_inf, 2) / gamma;
	double x_mom = 1 * Mach_number;
	double y_mom = 1 * Mach_number;
	double totalE_per_volume = rho_inf * 0.5 * pow(Mach_number, 2) + p_inf / (gamma - 1);

	int cell_loc = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			cell_loc = 4 * (icell + num_center_x * jcell);



			Solution[cell_loc + 0] = rho_inf;
			Solution[cell_loc + 1] = 1;
			Solution[cell_loc + 2] = 1;
			Solution[cell_loc + 3] = totalE_per_volume;

		}
	}



}