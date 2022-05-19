#include "../header/Time_integration.h"

void Timeintegration::initialize() {

	Config& config = Config::get_instance();
	Grid& grid = Grid::get_instance();
	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center = grid.get_num_total_center_point();

	num_total_var = 4 * num_total_center;

	test = 0;



	Cell_area = grid.get_Cell_area();

	std::string Limiter_type = config.get_Limiter_type();

	if (Limiter_type == "Minmod") {
		Spatialdiscrete_ = std::make_shared<Minmod_limiter>();
	}
	else if (Limiter_type == "Superbee") {
		Spatialdiscrete_ = std::make_shared<Minmod_limiter>();
	}
	else if (Limiter_type == "Off") {
		Spatialdiscrete_ = std::make_shared<Off>();
	}

	intermediate_1.resize(num_total_var, 0.0);
	intermediate_2.resize(num_total_var, 0.0);
	intermediate_3.resize(num_total_var, 0.0);
	intermediate_4.resize(num_total_var, 0.0);

	Residual.resize(num_total_center * 4, 0);



	Spatialdiscrete_->initialize();

}

double Timeintegration::Cal_dt(std::vector<double>& Solution) {

	double dt = Spatialdiscrete_->Cal_dt(Solution);

	return dt;

}

void RK3::Update(std::vector<double>& Solution, double dt) {

	Spatialdiscrete_->Cal_Residual(Solution, Residual);

	int var_index = 0;
	int cell_index = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {
				
				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_1[var_index] = Solution[var_index] - (dt / Cell_area[cell_index]) * Residual[var_index];
			}
		}
	}

	Spatialdiscrete_->Cal_Residual(intermediate_1, Residual);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_2[var_index] = 0.75 * Solution[var_index] + 0.25 * intermediate_1[var_index] - 0.25 *( dt / Cell_area[cell_index]) * Residual[var_index];
			}
		}
	}

	Spatialdiscrete_->Cal_Residual(intermediate_2, Residual);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				Solution[var_index] =  (1.0/ 3.0) * Solution[var_index] + (2.0/3.0 )* intermediate_2[var_index] -( 2.0/3.0) *( dt / Cell_area[cell_index] )* Residual[var_index];
			}
		}
	}


}

void SSPRK::Update(std::vector<double>& Solution, double dt) {

	Spatialdiscrete_->Cal_Residual(Solution, Residual);

	int var_index = 0;
	int cell_index = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_1[var_index] = Solution[var_index] - 0.391752226571890 * dt / Cell_area[cell_index] * Residual[var_index];
			}
		}
	}

	Spatialdiscrete_->Cal_Residual(intermediate_1, Residual);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_2[var_index] = 0.444370493651235 * Solution[var_index] + 0.555629506348765 * intermediate_1[var_index] - 0.368410593050371 * dt / Cell_area[cell_index] * Residual[var_index];
			}
		}
	}

	Spatialdiscrete_->Cal_Residual(intermediate_2, Residual);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_3[var_index] = 0.620101851488403 * Solution[var_index] + 0.379898148511597 * intermediate_2[var_index] - 0.251891774271694 * dt / Cell_area[cell_index] * Residual[var_index];
			}
		}
	}


	Spatialdiscrete_->Cal_Residual(intermediate_3, Residual);
	std::vector<double> temp_residual = Residual;
	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				intermediate_4[var_index] = 0.178079954393132 * Solution[var_index] + 0.821920045606868 * intermediate_3[var_index] - 0.544974750228521 * dt / Cell_area[cell_index] * Residual[var_index];
			}
		}
	}

	Spatialdiscrete_->Cal_Residual(intermediate_4, Residual);
	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				cell_index = icell + num_center_x * jcell;
				var_index = var + 4 * cell_index;

				Solution[var_index] = 0.517231671970585 * intermediate_2[var_index] + 0.096059710526147 * intermediate_3[var_index] - 0.063692468666290 * dt / Cell_area[cell_index] * temp_residual[var_index] + 0.386708617503269 * intermediate_4[var_index] - 0.226007483236906 * dt / Cell_area[cell_index] * Residual[var_index];

			}
		}
	}
}