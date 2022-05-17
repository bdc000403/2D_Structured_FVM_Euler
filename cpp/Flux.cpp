#include "../header/Flux.h"

void FVS::initialize() {

	Grid& grid = Grid::get_instance();
	Config& config = Config::get_instance();

	num_grid_point_x = grid.get_num_grid_point_x();
	num_grid_point_y = grid.get_num_grid_point_y();
	num_total_grid_point = grid.get_num_total_grid_point();

	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	num_total_center_point = grid.get_num_total_center_point();

	num_ghost_x = num_center_x + 2;
	num_ghost_y = num_center_y + 2;

	BC_ = std::make_shared<Boundarycondition>();

}

void FVS::Cal_Flux(std::vector<double>& South_solution, std::vector<double>& East_solution, std::vector<double>& North_solution, std::vector<double>& West_solution, std::vector<double>& Flux) {

	int West_index = 0;
	int East_index = 0;
	int South_index = 0;
	int North_index = 0;
	int x_index = 0;

	// West Flux 계산

	int i_size = num_center_x + 1;
	int j_size = num_center_y;
	int total_size = i_size * j_size;

	std::vector<double> Left_solution(4 * total_size);
	std::vector<double> Right_solution(4 * total_size);
	std::vector<double> temp_flux(4 * total_size);

	for (int jcell = 0; jcell < j_size; jcell++) {						// West Flux
		for (int icell = 0; icell < i_size; icell++) {
			for (int var = 0; var < 4; var++) {
				// Left --> West 
				// Right --> East 

				West_index = var + 4 * (icell + 1 + num_ghost_x * (jcell + 1));
				x_index = var + 4 * (icell + num_center_x * jcell);
				East_index = var + 4 * (icell + num_ghost_x * (jcell + 1));


				Right_solution[x_index] = West_solution[West_index];
				Left_solution[x_index] = East_solution[East_index];

			}
		}
	}

	Cal_Flux_split(Left_solution, Right_solution, temp_flux, total_size);

	int face_index = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				face_index = var + 4 * (3 + 4 * (icell + num_center_x * jcell));		// west face index 
				West_index = var + 4 * (icell + (num_center_x)*jcell);

				Flux[face_index] = temp_flux[West_index];								// west 저장

				face_index = var + 4 * (1 + 4 * (icell + num_center_x * jcell));		// east face index 
				East_index = var + 4 * (icell + 1 + (num_center_x)*jcell);

				Flux[face_index] = temp_flux[East_index];								// west 저장

			}
		}
	}

	i_size = num_center_x;
	j_size = num_center_y + 1;
	total_size = i_size * j_size;

	for (int jcell = 0; jcell < j_size; jcell++) {						// South Flux
		for (int icell = 0; icell < i_size; icell++) {
			for (int var = 0; var < 4; var++) {
				// Left --> South 
				// Right --> North 

				North_index = var + 4 * (icell + 1 + num_ghost_x * jcell);
				x_index = var + 4 * (icell + num_center_x * jcell);
				South_index = var + 4 * (icell + 1 + num_ghost_x * (jcell + 1));


				Left_solution[x_index] = North_solution[North_index];
				Right_solution[x_index] = South_solution[South_index];

			}
		}
	}

	Cal_Flux_split(Left_solution, Right_solution, temp_flux, total_size);

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < 4; var++) {

				face_index = var + 4 * (0 + 4 * (icell + num_center_x * jcell));		// South face index 
				South_index = var + 4 * (icell + num_center_x * jcell);

				Flux[face_index] = temp_flux[South_index];								// South 저장

				face_index = var + 4 * (2 + 4 * (icell + num_center_x * jcell));		// north face index
				North_index = var + 4 * (icell + num_center_x * (jcell + 1));

				Flux[face_index] = temp_flux[North_index];								// north 저장

			}

		}

	}

}

void VanLeer::Cal_Flux_split(std::vector<double>& Left_solution, std::vector<double>& Right_solution, std::vector<double>& temp_flux, int total_size) {

	uL_.resize(total_size, 0);
	uR_.resize(total_size, 0);
	vL_.resize(total_size, 0);
	vR_.resize(total_size, 0);
	pL_.resize(total_size, 0);
	pR_.resize(total_size, 0);
	aR_.resize(total_size, 0);
	aL_.resize(total_size, 0);
	HR_.resize(total_size, 0);
	HL_.resize(total_size, 0);

	double r_ = 1.4;


	for (int icell = 0; icell < total_size; icell++) {
		uL_[icell] = (Left_solution[icell * 4 + 1] / Left_solution[icell * 4 + 0]);
		uR_[icell] = (Right_solution[icell * 4 + 1] / Right_solution[icell * 4 + 0]);
		vL_[icell] = (Left_solution[icell * 4 + 2] / Left_solution[icell * 4 + 0]);
		vR_[icell] = (Right_solution[icell * 4 + 2] / Right_solution[icell * 4 + 0]);
		pR_[icell] = (r_ - 1) * (Right_solution[icell * 4 + 3] - 0.5 * ((pow(Right_solution[icell * 4 + 1], 2.0) + pow(Right_solution[icell * 4 + 2], 2.0)) / Right_solution[icell * 4 + 0]));
		pL_[icell] = (r_ - 1) * (Left_solution[icell * 4 + 3] - 0.5 * ((pow(Left_solution[icell * 4 + 1], 2.0) + pow(Left_solution[icell * 4 + 2], 2.0)) / Left_solution[icell * 4 + 0]));
		if (pR_[icell] < 0)
			pR_[icell] = (r_ - 1) * (Right_solution[icell * 4 + 3] - 0.5 * ((pow(Right_solution[icell * 4 + 1], 2.0) + pow(Right_solution[icell * 4 + 2], 2.0)) / Right_solution[icell * 4 + 0])) + pow(10.0, -10.0);
		else if (pL_[icell] < 0)
			pL_[icell] = (r_ - 1) * (Left_solution[icell * 4 + 3] - 0.5 * ((pow(Left_solution[icell * 4 + 1], 2.0) + pow(Left_solution[icell * 4 + 2], 2.0)) / Left_solution[icell * 4 + 0])) + pow(10.0, -10.0);
	}

	for (int icell = 0; icell < total_size; icell++) {
		aR_[icell] = sqrt(r_ * pR_[icell] / Right_solution[icell * 4 + 0]);
		aL_[icell] = sqrt(r_ * pL_[icell] / Left_solution[icell * 4 + 0]);
		HR_[icell] = ((Right_solution[icell * 4 + 3] + pR_[icell]) / Right_solution[icell * 4 + 0]);
		HL_[icell] = ((Left_solution[icell * 4 + 3] + pL_[icell]) / Left_solution[icell * 4 + 0]);
	}

	mach_P.resize(total_size, 0);
	mach_M.resize(total_size, 0);
	mach_L.resize(total_size, 0);
	mach_R.resize(total_size, 0);
	mach_interface.resize(total_size, 0);
	fluxP_.resize(total_size * 4, 0);
	fluxM_.resize(total_size * 4, 0);

	for (int icell = 0; icell < total_size; icell++) {
		mach_L[icell] = uL_[icell] / aL_[icell];		// i
		mach_R[icell] = uR_[icell] / aR_[icell];		// i 
		if (abs(mach_L[icell]) <= 1) {
			mach_P[icell] = 0.25 * (mach_L[icell] + 1) * (mach_L[icell] + 1);
		}
		else {
			mach_P[icell] = 0.5 * (mach_L[icell] + abs(mach_L[icell]));
		}
		if (abs(mach_R[icell]) <= 1) {
			mach_M[icell] = -0.25 * (mach_R[icell] - 1) * (mach_R[icell] - 1);
		}
		else {
			mach_M[icell] = 0.5 * (mach_R[icell] - abs(mach_R[icell]));
		}
		mach_interface[icell] =  mach_P[icell] + mach_M[icell];
	}

	for (int icell = 0; icell < total_size; icell++) {

		if (abs(mach_interface[icell]) < 1) {
			fluxP_[icell * 4 + 0] = 0.25 * Left_solution[icell * 4 + 0] * aL_[icell] * pow(mach_L[icell] + 1, 2.0);
			fluxP_[icell * 4 + 1] = fluxP_[icell * 4 + 0] * (2 * aL_[icell] / r_) * (0.5 * (r_ - 1.) * mach_L[icell] + 1.0);
			fluxP_[icell * 4 + 2] = fluxP_[icell * 4 + 0] * vL_[icell];
			fluxP_[icell * 4 + 3] = fluxP_[icell * 4 + 0] * ((1. / (2.0 * (r_ * r_ - 1.0))) * pow((r_ - 1) * mach_L[icell] * aL_[icell] + 2.0 * aL_[icell], 2.0) + 0.5 * pow(vL_[icell], 2));
			fluxM_[icell * 4 + 0] = - 0.25 * Right_solution[icell * 4 + 0] * aR_[icell] * pow(mach_R[icell] - 1, 2.0);
			fluxM_[icell * 4 + 1] = fluxM_[icell * 4 + 0] * (2 * aR_[icell] / r_) * (0.5 * (r_ - 1.) * mach_R[icell] - 1.0);
			fluxM_[icell * 4 + 2] = fluxM_[icell * 4 + 0] * vR_[icell];
			fluxM_[icell * 4 + 3] = fluxM_[icell * 4 + 0] * ((1. / (2.0 * (r_ * r_ - 1.0))) * pow((r_ - 1) * mach_R[icell] * aR_[icell] - 2.0 * aR_[icell], 2.0) + 0.5 * pow(vR_[icell], 2));
		}
		else if (mach_interface[icell] >= 1) {

			fluxP_[icell * 4 + 0] = Left_solution[icell * 4 + 1];
			fluxP_[icell * 4 + 1] = (pow(Left_solution[icell * 4 + 1], 2) / Left_solution[icell * 4 + 0]) + pL_[icell];
			fluxP_[icell * 4 + 2] = ((Left_solution[icell * 4 + 1] * Left_solution[icell * 4 + 2]) / Left_solution[icell * 4 + 0]);
			fluxP_[icell * 4 + 3] = Left_solution[icell * 4 + 1] * HL_[icell];
			fluxM_[icell * 4 + 0] = 0;
			fluxM_[icell * 4 + 1] = 0;
			fluxM_[icell * 4 + 2] = 0;
			fluxM_[icell * 4 + 3] = 0;
		}
		else {
			fluxP_[icell * 4 + 0] = 0;
			fluxP_[icell * 4 + 1] = 0;
			fluxP_[icell * 4 + 2] = 0;
			fluxP_[icell * 4 + 3] = 0;
			fluxM_[icell * 4 + 0] = Right_solution[icell * 4 + 1];
			fluxM_[icell * 4 + 1] = (pow(Right_solution[icell * 4 + 1], 2) / Right_solution[icell * 4 + 0]) + pR_[icell];
			fluxM_[icell * 4 + 2] = ((Right_solution[icell * 4 + 1] * Right_solution[icell * 4 + 2]) / Right_solution[icell * 4 + 0]);
			fluxM_[icell * 4 + 3] = Right_solution[icell * 4 + 1] * HR_[icell];
		}
	}

	for (int icell = 0; icell <  4 * (total_size); icell++) {
	
		temp_flux[icell] = fluxP_[icell] + fluxM_[icell];

	}




}

