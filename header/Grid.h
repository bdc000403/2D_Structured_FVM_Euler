#pragma once

#include "Config.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>


class Grid {

private:

	static Grid instance;

	std::string Grid_file;

	int num_grid_point_x;

	int num_grid_point_y;

	int num_total_grid_point;

	int num_center_x;

	int num_center_y;

	int num_total_center_point;
	
	std::vector<double> Grid_point_coordinates;

	std::vector<double> Cell_center_coordinates_x;

	std::vector<double> Cell_center_coordinates_y;

	std::vector<double> Cell_area;

	std::vector<double> Surface_length; // 4 * totlal_center - south 부터 반시계방향

	std::vector<double> nx; // 4 * total_center

	std::vector<double> ny; // 4* total_center

	std::vector<int> BC; // imin imax jmin jmax 순서 (start_index, end_index, type)

	std::ifstream inFile;

	std::vector<double> dx;

	std::vector<double> dy;


public:

	static Grid& get_instance() {		// static 함수 -> Class의 소유, static 멤버 instance 반환

		static Grid instance;

		return instance;
	}

	std::vector<double>& get_Grid_point_coordinates() {
		return Grid_point_coordinates;
	}
	
	int get_num_grid_point_x() {
		return num_grid_point_x;
	}

	int get_num_grid_point_y() {
		return num_grid_point_y;
	}

	int get_num_total_grid_point() {
		return num_total_grid_point;
	}

	int get_num_center_x() {
		return num_center_x;
	}

	int get_num_center_y() {
		return num_center_y;
	}

	int get_num_total_center_point() {
		return num_total_center_point;
	}

	std::vector<double>& get_Cell_center_coordinates_x() {
		return Cell_center_coordinates_x;
	}

	std::vector<double>& get_Cell_center_coordinates_y() {
		return Cell_center_coordinates_y;
	}

	std::vector<double>& get_Cell_area() {
		return Cell_area;
	}

	std::vector<double>& get_Surface_length() {
		return Surface_length;
	}

	std::vector<double>& get_nx() {
		return nx;
	}

	std::vector<double>& get_ny() {
		return ny;
	}

	std::vector<int>& get_BC() {
		return BC;
	}

	std::vector<double>& get_dx() {
		return dx;
	}

	std::vector<double>& get_dy() {
		return dy;
	}

	bool open_grid_file();

	void read_grid_point();

	void cal_grid_info();

	bool initialize();

protected:

	Grid() {};
	~Grid() {};


};