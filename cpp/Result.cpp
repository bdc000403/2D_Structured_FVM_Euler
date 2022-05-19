#include "../header/Result.h"


void Result::Initialize() {

	Config& config = Config::get_instance();
	Grid& grid = Grid::get_instance();
	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();
	nx = grid.get_nx();
	ny = grid.get_ny();
	AOA = config.get_AOA();
	AOA = AOA * 3.141592653589 / 180;
	Free_mach = config.get_Mach_number();
	Surface_length = grid.get_Surface_length();
	num_var = 4;
	Prev_solution.resize(num_center_x * num_center_y * 4, 0.0);
	file_name = config.get_Result_file();
	conv_file_name = file_name + "_conv_histroy.plt";
	cp_file_name = file_name + "_cp.plt";
	Error = 0;

	conv_file.open(conv_file_name);
	conv_file.precision(10);
	conv_file << "variables =	step,time,cl,cd,error" << std::endl;
	conv_file << "zone t =	'convhist'" << std::endl;

}

void Result::Post(std::vector<double>& Solution, double time, std::string file_name) {

	std::ofstream outfile(file_name);

	if (!outfile.is_open()) {
		std::cout << "Fail to open post file : " << std::endl;
		return;
	}

	int total = num_center_x * num_center_y;
	int index = 0;
	int sol_index = 0;
	double r_ = 1.4;

	Grid& grid = Grid::get_instance();

	std::vector<double> rho(total);
	std::vector<double> p(total);
	std::vector<double> mach(total);
	std::vector<double> u(total);
	std::vector<double> v(total);
	std::vector<double> a(total);
	x = grid.get_Cell_center_coordinates_x();
	y = grid.get_Cell_center_coordinates_y();
	std::vector<double> grid_coordinate = grid.get_Grid_point_coordinates();

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			index = icell + num_center_x * jcell;
			sol_index = 4 * index;

			rho[index] = Solution[sol_index + 0];
			u[index] = Solution[sol_index + 1] / rho[index];
			v[index] = Solution[sol_index + 2] / rho[index];
			p[index] = (r_ - 1) * (Solution[sol_index + 3] - 0.5 * rho[index] * (pow(u[index], 2) + pow(v[index], 2)));
			a[index] = pow(r_ * p[index] / rho[index], 0.5);
			mach[index] = u[index] / a[index];

		}
	}

	outfile.precision(10);
	outfile << "variables =	x,y,rho,p,mach,u,v" << std::endl;
	outfile << "zone t = '2D_structured'" << "," << "i= " << "\t" << (num_center_x + 1) << " ,j= " << "\t" << (num_center_y) << " ,f=point" << std::endl;
	outfile << "Solutiontime =" << time << std::endl;
	outfile.precision(10);
	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			index = icell + num_center_x * jcell;
			outfile << "\t" << x[index] << "\t" << y[index] << "\t" << rho[index] << "\t" << p[index] << "\t" << mach[index] << "\t" << u[index] << "\t" << v[index] << std::endl;
		}
		index = 0 + num_center_x * jcell;
		outfile << "\t" << x[index] << "\t" << y[index] << "\t" << rho[index] << "\t" << p[index] << "\t" << mach[index] << "\t" << u[index] << "\t" << v[index] << std::endl;


	}

	outfile.close();
}

void Result::Save(std::vector<double>& Solution, double time, int file_index) {

	std::string file;
	file = file_name + "_" + std::to_string(file_index) + ".plt";

	Post(Solution, time, file);

}

void Result::Save(std::vector<double>& Solution, double time, std::string file_index) {

	std::string file;
	file = file_name + "_" + file_index + ".plt";

	Post(Solution, time, file);

}

void Result::Save_conv_hist(std::vector<double>& Solution, int step, double time, double Error) {

	Cal_coefficient(Solution);

	conv_file << "\t" << step << "\t" << time << "\t" << Cl << "\t" << Cd << "\t" << Error << std::endl;

}

void Result::Cal_coefficient(std::vector<double>& Solution) {

	double p = 0;
	double rho_inf = 1;
	int index = 0;
	int face_index = 0;
	double r_ = 1.4;

	Fx = 0;
	Fy = 0;

	for (int icell = 0; icell < num_center_x; icell ++ ) {

		index = num_var * icell;
		face_index = 0 + 4 * (icell);
		p = (r_ - 1) * (Solution[index + 3] - 0.5 * (pow(Solution[index + 1], 2) + pow(Solution[index + 2], 2)) / Solution[index + 0]);

		Fx += p * nx[face_index] * Surface_length[face_index];
		Fy += p * ny[face_index] * Surface_length[face_index];

	}

	Drag = Fy * sin(AOA) + Fx * cos(AOA);
	Lift = Fy * cos(AOA) - Fx * sin(AOA);
	Cl = Lift / (0.5 * rho_inf * pow(Free_mach, 2));
	Cd = Drag / (0.5 * rho_inf * pow(Free_mach, 2));

}

void Result::Save_cp(std::vector<double>& Solution) {

	std::ofstream cp_file(cp_file_name);
	cp_file.precision(10);
	cp_file << "variables =	x,y," << '"' <<  "- cp" << '"' << std::endl;
	cp_file << "zone t =	'-cp'" << std::endl;

	double Cp;

	double p = 0;
	double r_ = 1.4;
	double rho_inf = 1;
	double sos_inf = 1;
	double p_inf = rho_inf * pow(sos_inf, 2) / r_;
	int index = 0;
	int face_index = 0;


	for (int icell = 0; icell < num_center_x; icell++) {

		index = num_var * icell;
		face_index = 0 + 4 * (icell);
		p = (r_ - 1) * (Solution[index + 3] - 0.5 * (pow(Solution[index + 1], 2) + pow(Solution[index + 2], 2)) / Solution[index + 0]);

		Cp = -(p - p_inf)/ (0.5 * rho_inf * pow(Free_mach, 2));

		index = 0.25 * index;

		cp_file << "\t" << x[index] << "\t" << y[index] << "\t" << Cp << std::endl;

	}



}

double Result::Cal_Error(std::vector<double>& Solution) {
	
	int index = 0;

	Error = 0;

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int var = 0; var < num_var; var++) {
				index = var + num_var * (icell + num_center_x * jcell);

				Error += pow(Solution[index] - Prev_solution[index] , 2);

			}
		}
	}

	Prev_solution = Solution;

	Error = Error / (static_cast<double>(num_center_x) * static_cast<double>(num_center_y) * static_cast<double>(num_var));
	Error = sqrt(Error);

	return Error;
}

void Result::Close_file() {

	conv_file.close();

}

void Result::Test_save(std::vector<double>& Solution, double time, int file_index) {

	std::string file;
	file = file_name + "_" + std::to_string(file_index) + ".plt";

	std::ofstream outfile(file);

	if (!outfile.is_open()) {
		std::cout << "Fail to open post file : " << std::endl;
		return;
	}

	int total = num_center_x * num_center_y;
	int index = 0;
	int sol_index = 0;
	double r_ = 1.4;

	Grid& grid = Grid::get_instance();

	std::vector<double> rho(total);
	std::vector<double> p(total);
	std::vector<double> mach(total);
	std::vector<double> u(total);
	std::vector<double> v(total);
	std::vector<double> a(total);
	std::vector<double> x = grid.get_Cell_center_coordinates_x();
	std::vector<double> y = grid.get_Cell_center_coordinates_y();
	std::vector<double> grid_coordinate = grid.get_Grid_point_coordinates();

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			index = icell + num_center_x * jcell;
			sol_index = 4 * index;

			rho[index] = Solution[sol_index + 0];
			u[index] = Solution[sol_index + 1] / rho[index];
			v[index] = Solution[sol_index + 2] / rho[index];
			p[index] = (r_ - 1) * (Solution[sol_index + 3] - 0.5 * rho[index] * (pow(u[index], 2) + pow(v[index], 2)));
			a[index] = pow(r_ * p[index] / rho[index], 0.5);
			mach[index] = u[index] / a[index];

		}
	}

	outfile.precision(10);
	outfile << "variables =	x,y,rho,p,mach,u,v" << std::endl;
	outfile << "zone t = '2D_structured'" << "," << "i= " << "\t" << (num_center_x + 1) << " ,j= " << "\t" << (num_center_y) << " ,f=point" << std::endl;

	outfile.precision(10);
	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			index = icell + num_center_x * jcell;
			outfile << "\t" << x[index] << "\t" << y[index] << "\t" << rho[index] << "\t" << p[index] << "\t" << mach[index] << "\t" << u[index] << "\t" << v[index] << std::endl;
		}
		index = 0 + num_center_x * jcell;
		outfile << "\t" << x[index] << "\t" << y[index] << "\t" << rho[index] << "\t" << p[index] << "\t" << mach[index] << "\t" << u[index] << "\t" << v[index] << std::endl;


	}

	outfile.close();
}
