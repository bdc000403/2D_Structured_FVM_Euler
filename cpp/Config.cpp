#include "../header/Config.h"

bool Config::read_config() {

	std::string config_file_name = "config/config.dat";
	std::ifstream inFile;

	inFile.open(config_file_name.c_str());

	if (inFile.is_open()) {
		std::cout << "*Configuration file was opend successfully" << std::endl << std::endl;
	}
	else {
		std::cout << "*Can't open Configuratio file" << std::endl << std::endl;
	}

	int num_config = 0;
	int total_config = 18;

	std::string text;
	while (std::getline(inFile, text)) { // std::getline --> inFile���� ���� �о���� text�� ���� +  ture ��ȯ + �����ٷ� Ŀ��  (������ ������ false ��ȯ)
		num_config += search_and_save(inFile, text, "Grid file=", Grid_file);
		num_config += search_and_save(inFile, text, "Max iteration=", &Max_iter);
		num_config += search_and_save(inFile, text, "Result file=", Result_file);
		num_config += search_and_save(inFile, text, "Save interval=", &Save_interval);
		num_config += search_and_save(inFile, text, "Flow type=", Flow_type);
		num_config += search_and_save(inFile, text, "Solution type=", Solution_type);
		num_config += search_and_save(inFile, text, "Tolerance=", &Tolerance);
		num_config += search_and_save(inFile, text, "Time interval=", &Time_interval);
		num_config += search_and_save(inFile, text, "Mach number=", &Mach_number);
		num_config += search_and_save(inFile, text, "Reynolds number=", &Reynolds_number);	//	double �Ǵ� int --> & �ٿ��� �ּҿ� ����
		num_config += search_and_save(inFile, text, "End time=", &End_time);
		num_config += search_and_save(inFile, text, "Flux type=", Flux_type);
		num_config += search_and_save(inFile, text, "Flux splitting type=", Flux_splitting_type);
		num_config += search_and_save(inFile, text, "Flux difference type=", Flux_difference_type);
		num_config += search_and_save(inFile, text, "Limiter type=", Limiter_type);
		num_config += search_and_save(inFile, text, "Time integration method=", Time_integration_method);
		num_config += search_and_save(inFile, text, "CFL number=", &CFL_number);
		num_config += search_and_save(inFile, text, "AOA=", &AOA);
	}	
	if (num_config == total_config) {
		std::cout << std::endl;
		return 1;
	}
	else {
		std::cout << "		*Insufficient configuration information" << std::endl << std::endl;
		std::cout << num_config;
		return 0;
	}

	inFile.close();

}



bool Config::search_and_save(std::ifstream& inFile, const std::string& text, const std::string& data_type, std::string& data_variable) {
	if (text.find(data_type, 0) != std::string::npos) {		// ã���� ù ��ġ ��ȯ
		getline(inFile, data_variable);						// ���� �о���� ����
		std::cout << data_type << " " << data_variable << std::endl;
		return 1;
	}
	else {
		return 0;
	}
}

bool Config::search_and_save(std::ifstream& inFile, const std::string& text, const std::string& data_type, int* data_variable) {
	if (text.find(data_type, 0) != std::string::npos) {		// ã���� ù ��ġ ��ȯ
		inFile >> *data_variable;					// ���� �о���� ����
		std::cout << data_type << " " << std::to_string(*data_variable) << std::endl;
		return 1;
	}
	else {
		return 0;
	}
}

bool Config::search_and_save(std::ifstream& inFile, const std::string& text, const std::string& data_type, double* data_variable) {
	if (text.find(data_type, 0) != std::string::npos) {		// ã���� ù ��ġ ��ȯ
		inFile >> *data_variable;						// ���� �о���� ����
		std::cout << data_type << " " << std::to_string(*data_variable) << std::endl;
		return 1;
	}
	else {
		return 0;
	}
}