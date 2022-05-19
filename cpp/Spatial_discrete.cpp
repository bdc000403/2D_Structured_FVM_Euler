#include "../header/Spatial_discrete.h"

void Spatialdiscrete::initialize() {

	Grid& grid = Grid::get_instance();
	BC = grid.get_BC();

	Config& config = Config::get_instance();

	num_grid_point_x = grid.get_num_grid_point_x();
	num_grid_point_y = grid.get_num_grid_point_y();
	num_total_grid_point = grid.get_num_total_grid_point();

	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center_point = grid.get_num_total_center_point();

	num_ghost_x = num_center_x + 2;
	num_ghost_y = num_center_y + 2;

	dx = grid.get_dx();
	dy = grid.get_dy();

	CFL = config.get_CFL();
	Cell_area = grid.get_Cell_area();

	nx = grid.get_nx();
	ny = grid.get_ny();

	AOA = config.get_AOA();

	Mach_number = config.get_Mach_number();

	North_solution.resize(num_ghost_x * num_ghost_y * 4, 0.0);
	South_solution.resize(num_ghost_x * num_ghost_y * 4, 0.0);
	West_solution.resize(num_ghost_x * num_ghost_y * 4, 0.0);
	East_solution.resize(num_ghost_x * num_ghost_y * 4, 0.0);
	Ghost_solution.resize((num_ghost_x) * (num_ghost_y) * 4, 0.0);
	Flux.resize(num_center_x * num_center_y * 4 * 4, 0.0);

	x_phi.resize(num_center_x * num_center_y * 4, 0.0);
	y_phi.resize(num_center_x * num_center_y * 4, 0.0);
	x_theta.resize(num_center_x * num_center_y * 4, 0.0);
	y_theta.resize(num_center_x * num_center_y * 4, 0.0);

	Surface_length = grid.get_Surface_length();

	BC_ = std::make_shared<Boundarycondition>();

	BC_->BC_initialize();

	std::string FVS = config.get_Flux_splitting_type();

	if (FVS == "VanLeer") {
		FVS_ = std::make_shared<VanLeer>();
	}
	else if (FVS == "AUSM") {

	}

	FVS_->initialize();

}

void Spatialdiscrete::Cal_Residual(std::vector<double>& Solution, std::vector<double>& Residual) {


	Cal_interpolation(Solution);
	Cal_coordinate_transform();

	FVS_->Cal_Flux(South_solution, East_solution, North_solution, West_solution, Flux);

	Cal_inverse_coordinate_transform();							
	BC_->Flux_boundary_condition(Solution, Flux);

	int Flux_index = 0;
	int cell_index = 0;
	int face_index = 0;


	Residual.resize(num_total_center_point * 4, 0);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			
			cell_index = 4 * (icell + num_center_x * jcell);

			Residual[cell_index + 0] = 
				- Flux[0 + 4 * (0 + cell_index)] * Surface_length[0 + cell_index]	//south
				+ Flux[0 + 4 * (1 + cell_index)] * Surface_length[1 + cell_index]							//east
				+ Flux[0 + 4 * (2 + cell_index)] * Surface_length[2 + cell_index]							//north
				- Flux[0 + 4 * (3 + cell_index)] * Surface_length[3 + cell_index];					//west

			Residual[cell_index + 1] = 
				- Flux[1 + 4 * (0 + cell_index)] * Surface_length[0 + cell_index]	//south
				+ Flux[1 + 4 * (1 + cell_index)] * Surface_length[1 + cell_index]							//east
				+ Flux[1 + 4 * (2 + cell_index)] * Surface_length[2 + cell_index]							//north
				- Flux[1 + 4 * (3 + cell_index)] * Surface_length[3 + cell_index];						//west

			Residual[cell_index + 2] = 
				- Flux[2 + 4 * (0 + cell_index)] * Surface_length[0 + cell_index]	//south
				+ Flux[2 + 4 * (1 + cell_index)] * Surface_length[1 + cell_index]							//east
				+ Flux[2 + 4 * (2 + cell_index)] * Surface_length[2 + cell_index]							//north
				- Flux[2 + 4 * (3 + cell_index)] * Surface_length[3 + cell_index];					//west

			Residual[cell_index + 3] = 
				- Flux[3 + 4 * (0 + cell_index)] * Surface_length[0 + cell_index]	//south
				+ Flux[3 + 4 * (1 + cell_index)] * Surface_length[1 + cell_index]							//east
				+ Flux[3 + 4 * (2 + cell_index)] * Surface_length[2 + cell_index]							//north
				- Flux[3 + 4 * (3 + cell_index)] * Surface_length[3 + cell_index];					//west

	
		}
	}
}

double Spatialdiscrete::Cal_dt(std::vector<double>& Solution) {

	int cell_index = 0;
	double spectral_radii_I;
	double spectral_radii_J;
	double dt = 1000;
	double n_I_x = 0;
	double n_J_x = 0;
	double n_I_y = 0;
	double n_J_y = 0;
	double del_s_I = 0;
	double del_s_J = 0;
	double p = 0;
	double a = 0;
	double u = 0;
	double v = 0;
	double r_ = 1.4;
	double dt_temp = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			cell_index = 4 * (icell + num_center_x * jcell);

			u = Solution[cell_index + 1] / Solution[cell_index + 0];
			v = Solution[cell_index + 2] / Solution[cell_index + 0];
			p = (r_ - 1) * (Solution[cell_index + 3] - 0.5 * Solution[cell_index + 0] * (pow(u,2) + pow(v,2)));
			a = pow(r_ * p / Solution[cell_index + 0], 0.5);

			n_I_x = 0.5 * (nx[cell_index + 3] - nx[cell_index + 1]);
			n_I_y = 0.5 * (ny[cell_index + 3] - ny[cell_index + 1]);
			n_J_x = 0.5 * (nx[cell_index + 0] - nx[cell_index + 2]);
			n_J_y = 0.5 * (ny[cell_index + 0] - ny[cell_index + 2]);

			del_s_I = 0.5 * (Surface_length[cell_index + 3] + Surface_length[cell_index + 1]);
			del_s_J = 0.5 * (Surface_length[cell_index + 0] + Surface_length[cell_index + 2]);

			spectral_radii_I = (abs(u * n_I_x + v * n_I_y) + a) * del_s_I;
			spectral_radii_J = (abs(u * n_J_x + v * n_J_y) + a) * del_s_J;

			cell_index = 0.25 * cell_index;
			dt_temp = CFL * Cell_area[cell_index] / (spectral_radii_I + spectral_radii_J);

			if (dt_temp < dt) {
				dt = dt_temp;
			}

		}
	}

	return dt;

}

void Spatialdiscrete::Cal_interpolation(std::vector<double>& Solution) {

	int ghost_index = 0;
	int cell_index = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				ghost_index = var + 4 * ((icell + 1) + num_ghost_x * (jcell + 1));
				cell_index = var + 4 * (icell + num_center_x * jcell);

				Ghost_solution[ghost_index] = Solution[cell_index];

			}
		}
	}

	BC_->Cell_boundary_condition(Ghost_solution);

	Cal_smoothness();
	Cal_phi();

	double del_x = 0;
	double del_y = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			for (int var = 0; var < 4; var++) {

				cell_index = var + 4 * (icell + num_center_x * jcell);
				ghost_index = var + 4 * ((icell + 1) + num_ghost_x * (jcell + 1));

				del_x = Ghost_solution[ghost_index] - Ghost_solution[ghost_index - 4];
				del_y = Ghost_solution[ghost_index] - Ghost_solution[ghost_index - 4 * num_ghost_x];

				North_solution[ghost_index] = Ghost_solution[ghost_index] + 0.5 * y_phi[cell_index] * del_y;
				South_solution[ghost_index] = Ghost_solution[ghost_index] - 0.5 * y_phi[cell_index] * del_y;
				East_solution[ghost_index] = Ghost_solution[ghost_index] + 0.5 * x_phi[cell_index] * del_x;
				West_solution[ghost_index] = Ghost_solution[ghost_index] - 0.5 * x_phi[cell_index] * del_x;



			}
		}
	}
	//  ghost cell 

	BC_->Cell_boundary_condition(North_solution);
	BC_->Cell_boundary_condition(South_solution);
	BC_->Cell_boundary_condition(East_solution);
	BC_->Cell_boundary_condition(West_solution);

}

void Spatialdiscrete::Cal_coordinate_transform() {

	int ghost_index = 0;
	int cell_index = 0;

	double rho_u = 0;
	double rho_v = 0;
	double cos = 0;
	double sin = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

				cell_index = 4 * (icell + num_center_x * jcell);
				ghost_index = 4 * ((icell + 1) + num_ghost_x * (jcell + 1));


				rho_u = North_solution[ghost_index + 1];
				rho_v = North_solution[ghost_index + 2];
				cos = nx[2 + cell_index];
				sin = ny[2 + cell_index];
				North_solution[ghost_index + 1] = rho_u * cos + rho_v * sin;
				North_solution[ghost_index + 2] = -rho_u * sin + rho_v * cos;


				rho_u = South_solution[ghost_index + 1];
				rho_v = South_solution[ghost_index + 2];
				cos = -nx[0 + cell_index];
				sin = -ny[0 + cell_index];
				South_solution[ghost_index + 1] = rho_u * cos + rho_v * sin;
				South_solution[ghost_index + 2] = -rho_u * sin + rho_v * cos;


				rho_u = East_solution[ghost_index + 1];
				rho_v = East_solution[ghost_index + 2];
				cos = nx[1 + cell_index];
				sin = ny[1 + cell_index];
				East_solution[ghost_index + 1] = rho_u * cos + rho_v * sin;
				East_solution[ghost_index + 2] = -rho_u * sin + rho_v * cos;


				rho_u = West_solution[ghost_index + 1];
				rho_v = West_solution[ghost_index + 2];
				cos = -nx[3 + cell_index];
				sin = -ny[3 + cell_index];
				West_solution[ghost_index + 1] = rho_u * cos + rho_v * sin;
				West_solution[ghost_index + 2] = -rho_u * sin + rho_v * cos;

		}
	}



	int face_index = 0;

								
	for (int jcell = 0; jcell < num_center_y; jcell++) {		//Imin --> east

		cell_index = 4 * (0 + num_ghost_x * (jcell + 1));
		face_index = 4 * (0 + num_center_x * (jcell));
		rho_u = East_solution[cell_index + 1];
		rho_v = East_solution[cell_index + 2];
		cos = -nx[3 + face_index];
		sin = -ny[3 + face_index];
		East_solution[cell_index + 1] = rho_u * cos + rho_v * sin;
		East_solution[cell_index + 2] = -rho_u * sin + rho_v * cos;


	}
										
	for (int jcell = 0; jcell < num_center_y; jcell++) {		//Imax --> west

		cell_index = 4 * (num_center_x + 1 + num_ghost_x * (jcell + 1));
		face_index = 4 * (num_center_x - 1 + num_center_x * (jcell));
		rho_u = West_solution[cell_index + 1];
		rho_v = West_solution[cell_index + 2];
		cos = nx[1 + face_index];
		sin = ny[1 + face_index];
		West_solution[cell_index + 1] = rho_u * cos + rho_v * sin;
		West_solution[cell_index + 2] = -rho_u * sin + rho_v * cos;

	
	}

	for (int icell = 0; icell < num_center_x; icell++) {		//Jmin --> North

		cell_index = 4 * (icell + 1);
		face_index = 4 * (icell);
		rho_u = North_solution[cell_index + 1];
		rho_v = North_solution[cell_index + 2];
		cos = -nx[0 + face_index];
		sin = -ny[0 + face_index];
		North_solution[cell_index + 1] = rho_u * cos + rho_v * sin;
		North_solution[cell_index + 2] = -rho_u * sin + rho_v * cos;


	}
										//Jmax
	for (int icell = 0; icell < num_center_x; icell++) {

		cell_index = 4 * (icell + 1 + num_ghost_x * (num_center_y + 1));
		face_index = 4 * (icell + num_center_x * (num_center_y -1));
		rho_u = South_solution[cell_index + 1];			
		rho_v = South_solution[cell_index + 2];			
		cos = nx[2 + face_index];
		sin = ny[2 + face_index];
		South_solution[cell_index + 1] = rho_u * cos + rho_v * sin;		
		South_solution[cell_index + 2] = -rho_u * sin + rho_v * cos;
		
	}


}

void Spatialdiscrete::Cal_inverse_coordinate_transform() {


	int Flux_index = 0;
	int cell_index = 0;

	double rho_u = 0;
	double rho_v = 0;
	double cos = 0;
	double sin = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			cell_index = 4 * (icell + num_center_x * jcell);
			Flux_index = 4 * (0 + 4 * (icell  + num_center_x * jcell));			// south flux

			rho_u = Flux[Flux_index + 1];
			rho_v = Flux[Flux_index + 2];
			cos = -nx[0 + cell_index];
			sin = -ny[0 + cell_index];

			Flux[Flux_index + 1] = rho_u * cos - rho_v * sin;
			Flux[Flux_index + 2] = rho_u * sin + rho_v * cos;

			Flux_index = 4 * (2 + 4 * (icell + num_center_x * jcell));			// north flux

			rho_u = Flux[Flux_index + 1];
			rho_v = Flux[Flux_index + 2];
			cos = nx[2 + cell_index];
			sin = ny[2 + cell_index];
			Flux[Flux_index + 1] = rho_u * cos - rho_v * sin;
			Flux[Flux_index + 2] = rho_u * sin + rho_v * cos;


			Flux_index = 4 * (1 + 4 * (icell + num_center_x * jcell));			// east flux

			rho_u = Flux[Flux_index + 1];
			rho_v = Flux[Flux_index + 2];
			cos = nx[1 + cell_index];
			sin = ny[1 + cell_index];

			Flux[Flux_index + 1] = rho_u * cos - rho_v * sin;
			Flux[Flux_index + 2] = rho_u * sin + rho_v * cos;

			Flux_index = 4 * (3 + 4 * (icell + num_center_x * jcell));			// west flux

			rho_u = Flux[Flux_index + 1];
			rho_v = Flux[Flux_index + 2];
			cos = -nx[3 + cell_index];
			sin = -ny[3 + cell_index];

			Flux[Flux_index + 1] = rho_u * cos - rho_v * sin;
			Flux[Flux_index + 2] = rho_u * sin + rho_v * cos;

		}
	}



}

void Spatialdiscrete::Cal_smoothness() {

	int cell_index = 0;
	int ghost_index = 0;
	double x_up = 0;
	double x_down = 0;
	double y_up = 0;
	double y_down = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = var + 4 * (icell + num_center_x * jcell);
				ghost_index = var + 4 * ((icell + 1) + num_ghost_x * (jcell + 1));

				x_up = Ghost_solution[ghost_index + 4] - Ghost_solution[ghost_index];
				x_down = Ghost_solution[ghost_index] - Ghost_solution[ghost_index - 4];

				y_up = Ghost_solution[ghost_index + 4 * num_ghost_x] - Ghost_solution[ghost_index];
				y_down = Ghost_solution[ghost_index] - Ghost_solution[ghost_index - 4 * num_ghost_x];

				if (abs(x_up) < 1e-6) {
					if (x_up > 0) {
						x_up = 1e-6;
					}
					else if (x_up < 0){
						x_up = -1e-6;
					}
				}
				if (abs(x_down) < 1e-6) {
					if (x_down > 0) {
						x_down = 1e-6;
					}
					else if (x_down < 0) {
						x_down = -1e-6;
					}
				}
				if (abs(y_up) < 1e-6) {
					if (y_up > 0) {
						y_up = 1e-6;
					}
					else if (y_up < 0) {
						y_up = -1e-6;
					}
				}
				if (abs(y_down) < 1e-6) {
					if (y_down > 0) {
						y_down = 1e-6;
					}
					else if (y_down < 0) {
						y_down = -1e-6;
					}
				}

				x_theta[cell_index] = x_up / x_down;
				y_theta[cell_index] = y_up / y_down;
			}
		}
	}
}


void Minmod_limiter::Cal_phi() {

	int cell_index = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = var + 4 * (icell + num_center_x * jcell);

				if (x_theta[cell_index] <= 0) {
					x_phi[cell_index] = 0;
				}
				else {
					x_phi[cell_index] = std::min(1.0, x_theta[cell_index]);
				}

				if (y_theta[cell_index] <= 0) {
					y_phi[cell_index] = 0;
				}
				else {
					y_phi[cell_index] = std::min(1.0, y_theta[cell_index]);
				}
			}
		}
	}
}

void Off::Cal_phi() {



}