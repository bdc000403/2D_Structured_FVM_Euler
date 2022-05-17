#include "../header/Solver.h"

void Solver::initialize() {

	Config& config = Config::get_instance();
	Grid& grid = Grid::get_instance();
	num_center_x = grid.get_num_center_x();
	num_center_y = grid.get_num_center_y();

	std::string Limiter_type = config.get_Limiter_type();
	std::string Time_integration = config.get_Time_integration_method();
	std::string flux_splitting_type = config.get_Flux_splitting_type();
	num_result_file = config.get_Num_result_file();
	Max_iter = config.get_Max_iter();
	End_time = config.get_End_tume();
	Save_interval = config.get_Save_interval();
	Tolerance = config.get_Tolerance();
	Time_interval = config.get_Time_interval();
	CFL = config.get_CFL();

	Solution.resize(num_center_x * num_center_y * 4, 0.0);
	Residual.resize(num_center_x * num_center_y * 4, 0);

	Initial_ = std::make_shared<Initial>();

	Result_ = std::make_shared<Result>();

	if (Time_integration == "SSPRK") {
		Timeintegration_ = std::make_shared<SSPRK>();
	}
	if (Time_integration == "RK3") {
		Timeintegration_ = std::make_shared<RK3>();
	}

	Initial_->Free_stream(Solution);

	//Initial_->Test(Solution);

	Timeintegration_->initialize();

	Result_->Initialize();
}

void Solver::Solve_Steady() {


	double time = 0.0;
	double dt = 0.0;
	double Error = 100;
	bool End_flag = 0;
	double Save_step = 0;
	int step = 0;
	int file_index = 1;

	Result_->Save(Solution, 0, 0);		// initial solution

	while (!End_flag) {

		dt = Timeintegration_->Cal_dt(Solution);

		if (isnan(dt)) {
			std::cout << "			\n\nInfinity loop\n" << std::endl;
			return;
		}

		time += dt;
		Save_step++;
		step++;
		Error = Result_->Cal_Error(Solution);

		std::cout << "Step: " << step << ",   dt: " << dt << ",   Error: " << Error << std::endl << std::endl;
		Timeintegration_->Update(Solution, dt);

		if (Save_step >= Save_interval) {
			Result_->Save(Solution, time, file_index);
			Result_->Save_conv_hist(Solution, step, time, Error);
			Save_step -= Save_interval;
			file_index++;
		}

		if ((Error <= Tolerance) || (step >= Max_iter))
		{
			End_flag = true;
		}

	}

	Result_->Save(Solution, time, "result");
	Result_->Save_conv_hist(Solution, step, time, Error);
	Result_->Save_cp(Solution);
	Result_->Close_file();

	return;
}

void Solver::Solve_Unsteady() {

	double time = 0.0;
	double dt = 0.0;
	bool End_flag = 0;
	double Save_time = 0;
	int step = 0;
	int file_index = 1;

	Result_->Save(Solution, 0, 0);		// initial solution

	while (!End_flag) {

		dt = Timeintegration_->Cal_dt(Solution);

		if (isnan(dt)) {
			std::cout << "			\n\nInfinity loop\n" << std::endl;
			return;
		}

		if (time + dt >= End_time)
		{
			dt = End_time - time;
			time = End_time;
			End_flag = true;

		}
		else {
			time += dt;
		}
		Save_time += dt;
		step++;

		std::cout << "Step: " << step << ",   dt: " << dt << ",   time: " << time << std::endl << std::endl;
		Timeintegration_->Update(Solution, dt);

		if (Save_time >= Time_interval) {
			Result_->Save(Solution, time, file_index);
			//Result_->Save_conv_hist(Solution, time);
			Save_time -= Time_interval;
			file_index++;

		}


	}

	Result_->Save(Solution, time, "result");
	//Result_->Save_conv_hist(Solution, time);
	Result_->Save_cp(Solution);
	Result_->Close_file();
	
	return;
}
