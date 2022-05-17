#include "../header/Config.h"
#include "../header/grid.h"
#include "../header/Solver.h"
#include <iostream>
#include <string>

int main() {

	std::cout << std::endl;
	std::cout << "	2D Euler Equation Solver for Structured Grid" << std::endl;

	Config& config = Config::get_instance();		
	if (config.read_config()) {}
	else { return 0; }

	Grid& grid = Grid::get_instance();
	if (grid.initialize()){}
	else { return 0; }

	std::string yesorno;
	std::cout << "	Continue to Solve? (y/n) : ";
	std::cin >> yesorno;
	if (yesorno != "y") {
		return 0;
	}

	std::shared_ptr<Solver> solver = std::make_shared<Solver>();
	solver->initialize();

	std::string Solution_type = config.get_Solution_type();

	if (Solution_type == "Unsteady") {
		solver->Solve_Unsteady();
	}
	else if (Solution_type == "Steady"){
		solver->Solve_Steady();
	}

	return 0;

}


