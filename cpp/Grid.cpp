#include "../header/Grid.h"


bool Grid::initialize() {

	if (open_grid_file() != 1) {
		return 0;
	}

	read_grid_point();
	cal_grid_info();


	return 1;
}

bool Grid::open_grid_file() {

	Config& config = Config::get_instance();
	Grid_file = config.get_Grid_file();

	std::string Grid_file_loc = "grid/" + Grid_file;


	inFile.open(Grid_file_loc.c_str());

	if (inFile.is_open()) {
		std::cout << "		*Grid file was opend successfully" << std::endl << std::endl;
	}
	else {
		std::cout << "		*Can't open Grid file" << std::endl << std::endl;
	}

	double temp = 0;

	inFile >> temp;
	if (temp != 1) {
		std::cout << "*Inappropriate Grid File (Number of Block must be 1)" << std::endl << std::endl;
		return 0;
	}

	inFile >> temp;
	num_grid_point_x = static_cast<int>(temp);
	num_center_x = num_grid_point_x - 1;

	inFile >> temp;
	num_grid_point_y = static_cast<int>(temp);
	num_center_y = num_grid_point_y - 1;

	num_total_grid_point = num_grid_point_x * num_grid_point_y;
	num_total_center_point = num_center_x * num_center_y;


	inFile >> temp;
	if (temp != 1) {
		std::cout << "		*Inappropriate Grid File (Grid must be in 2D format)" << std::endl << std::endl;
		return 0;
	}

	return 1;
}

void Grid::read_grid_point() {

	Grid_point_coordinates.resize(num_total_grid_point * 2);
	Cell_center_coordinates_x.resize(num_total_center_point);
	Cell_center_coordinates_y.resize(num_total_center_point);
	Cell_area.resize(num_total_center_point);
	nx.resize(num_total_center_point * 4);
	ny.resize(num_total_center_point * 4);
	Surface_length.resize(num_total_center_point * 4);
	dx.resize(num_total_center_point, 0);
	dy.resize(num_total_center_point, 0);


	double temp = 0;

	// Grid 저장

	for (int jnode = 0; jnode < num_grid_point_y; jnode++) {
		for (int inode = 0; inode < num_grid_point_x; inode++) {
			inFile >> temp;
			Grid_point_coordinates[inode + jnode * num_grid_point_x] = temp;
		}
	}

	for (int jnode = 0; jnode < num_grid_point_y; jnode++) {
		for (int inode = 0; inode < num_grid_point_x; inode++) {
			inFile >> temp;
			Grid_point_coordinates[num_total_grid_point + inode + jnode * num_grid_point_x] = temp;
		}
	}

	for (int i = 0; i < num_total_grid_point; i++) {
		inFile >> temp;			// Z 좌표 뛰어 넘기
		//std::cout << temp << i<<"  ";
	}

	std::string rank;		// "rank"문자 --> string에 저장해서 처리 (double temp에 저장하면 오류)
	inFile >> rank;
	inFile >> rank;

	// BC 읽어오기

	BC.resize(4 * 4);

	for (int BC_index = 0; BC_index < 4; BC_index++) {

		inFile >> temp;
		if (temp != 1.0) {
			std::cout << "Too many BC in BC_index = " << BC_index << "  " << temp << std::endl << std::endl;
		}
		

		for (int BC_info = 0; BC_info < 4; BC_info++) {		// BC_info, 1: start_index, 2: end_indx, 3: type, 4: additional

			inFile >> temp;
			BC[BC_info + BC_index * 4] = temp;

		}
	}

	inFile.close();

	for (int BC_index = 0; BC_index < 4; BC_index++) {
		for (int BC_info = 0; BC_info < 4; BC_info++) {		// BC_info, 1: start_index, 2: end_indx, 3: type, 4: additional

			//std::cout << BC[BC_info + BC_index * 4] << "   ";

		}
	}

}

void Grid::cal_grid_info(){

	int cell_loc = 0;
	int grid_loc = 0;

	double x_1 = 0;
	double x_2 = 0;
	double x_3 = 0;
	double x_4 = 0;
	double y_1 = 0;
	double y_2 = 0;
	double y_3 = 0;
	double y_4 = 0;


	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {

			cell_loc = icell + jcell * num_center_x;
			grid_loc = icell + jcell * num_grid_point_x;

			x_1 = Grid_point_coordinates[grid_loc];														// (i,j)
			x_2 = Grid_point_coordinates[grid_loc + 1];													// (i + 1,j)
			x_3 = Grid_point_coordinates[grid_loc + num_grid_point_x + 1];								// (i + 1,j + 1)
			x_4 = Grid_point_coordinates[grid_loc + num_grid_point_x];									// (i,j + 1)
			y_1 = Grid_point_coordinates[num_total_grid_point + grid_loc];								// (i,j)
			y_2 = Grid_point_coordinates[num_total_grid_point + grid_loc + 1];							// (i + 1,j)
			y_3 = Grid_point_coordinates[num_total_grid_point + grid_loc + num_grid_point_x + 1];		// (i + 1,j + 1)
			y_4 = Grid_point_coordinates[num_total_grid_point + grid_loc + num_grid_point_x];			// (i,j + 1)

			Cell_center_coordinates_x[cell_loc] = 0.25 * (x_1 + x_2 + x_3 + x_4);
			Cell_center_coordinates_y[cell_loc] = 0.25 * (y_1 + y_2 + y_3 + y_4);
			Cell_area[cell_loc] = 0.5 * ((x_1 - x_3) * (y_2 - y_4) + (x_4 - x_2) * (y_1 - y_3));

			// surf_index --> south부터 반시계방향
			
			Surface_length[0 + 4 * cell_loc] = pow(pow(y_2 - y_1, 2) + pow(x_1 - x_2, 2), 0.5);			//South
			Surface_length[1 + 4 * cell_loc] = pow(pow(y_3 - y_2, 2) + pow(x_2 - x_3, 2), 0.5);			//East
			Surface_length[2 + 4 * cell_loc] = pow(pow(y_4 - y_3, 2) + pow(x_3 - x_4, 2), 0.5);			//North
			Surface_length[3 + 4 * cell_loc] = pow(pow(y_1 - y_4, 2) + pow(x_4 - x_1, 2), 0.5);			//West

			nx[0 + 4 * cell_loc] = (y_2 - y_1) / Surface_length[0 + 4 * cell_loc];				//South
			ny[0 + 4 * cell_loc] = (x_1 - x_2) / Surface_length[0 + 4 * cell_loc];				//South
			nx[1 + 4 * cell_loc] = (y_3 - y_2) / Surface_length[1 + 4 * cell_loc];				//East
			ny[1 + 4 * cell_loc] = (x_2 - x_3) / Surface_length[1 + 4 * cell_loc];				//East
			nx[2 + 4 * cell_loc] = (y_4 - y_3) / Surface_length[2 + 4 * cell_loc];				//North
			ny[2 + 4 * cell_loc] = (x_3 - x_4) / Surface_length[2 + 4 * cell_loc];				//North
			nx[3 + 4 * cell_loc] = (y_1 - y_4) / Surface_length[3 + 4 * cell_loc];				//West
			ny[3 + 4 * cell_loc] = (x_4 - x_1) / Surface_length[3 + 4 * cell_loc];				//West




		}
	}

	for (int jcell = 0; jcell < num_center_y; jcell++) {
		for (int icell = 0; icell < num_center_x; icell++) {
			for (int face = 0; face < 4; face++) {

				cell_loc = face + 4 * (icell + num_center_x * jcell);

				if (abs(ny[cell_loc]) < pow(10, -14)) {
					ny[cell_loc] = 0;
				}

				if (abs(nx[cell_loc]) < pow(10, -14)) {
					nx[cell_loc] = 0;
				}


			}
		}
	}



}