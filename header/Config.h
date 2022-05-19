#pragma once

#include <string>
#include <iostream>
#include <fstream>


class Config {
private:

	std::string Grid_file;

	std::string Flow_type;

	double Mach_number;

	double Reynolds_number;

	double AOA;

	double End_time;

	int num_result_file;

	int Max_iter;

	int Save_interval;

	double CFL_number;

	std::string Flux_type;
	
	std::string Flux_splitting_type;

	std::string Flux_difference_type;

	std::string Limiter_type;

	std::string Result_file;

	std::string Time_integration_method;

	std::string Solution_type;

	double Tolerance;

	double Time_interval;

	static Config instance;

public:

	static Config& get_instance() {		// static 함수 -> Class의 소유, static 멤버 instance 반환

		static Config instance;

		return instance;
	}

	std::string& get_Grid_file() {
		return Grid_file;
	}

	std::string& get_Solution_type() {
		return Solution_type;
	}

	std::string& get_Result_file() {
		return Result_file;
	}

	std::string& get_Flow_type() {
		return Flow_type;
	}

	int& get_Max_iter() {
		return Max_iter;
	}

	int& get_Save_interval() {
		return Save_interval;
	}

	double& get_Mach_number() {
		return Mach_number;
	}

	double& get_Tolerance() {
		return Tolerance;
	}

	double& get_Time_interval() {
		return Time_interval;
	}

	double& get_Reynolds_number() {
		return Reynolds_number;
	}

	double& get_End_tume() {
		return End_time;
	}

	int& get_Num_result_file() {
		return num_result_file;
	}

	std::string& get_Flux_type() {
		return Flux_type;
	}

	std::string& get_Flux_splitting_type() {
		return Flux_splitting_type;
	}

	std::string& get_Flux_difference_type() {
		return Flux_difference_type;
	}

	std::string& get_Limiter_type() {
		return Limiter_type;
	}

	std::string& get_Time_integration_method() {
		return Time_integration_method;
	}

	double& get_AOA() {
		return AOA;
	}

	double& get_CFL() {
		return CFL_number;
	}

	bool read_config();

	bool search_and_save(std::ifstream& config_file, const std::string& text, const std::string& data_type, std::string& data_variable);

	bool search_and_save(std::ifstream& config_file, const std::string& text, const std::string& data_type, int* data_variable);

	bool search_and_save(std::ifstream& config_file, const std::string& text, const std::string& data_type, double* data_variable);

protected:

	Config() {};		// 생성자 protected --> 멤버 함수를 통해서 접근 (getinstance)

	~Config() {};

};