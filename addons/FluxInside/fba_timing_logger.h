#pragma once

#include <chrono>
#include <fstream>
#include <mutex>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

struct FbaTimingLogger
{
	std::mutex mutex;
	std::ofstream stream;
	std::chrono::steady_clock::time_point start_wall_clock;
	std::chrono::steady_clock::time_point previous_wall_clock;
	double previous_process_cpu_s = 0.0;
	bool has_previous_process_sample = false;
	bool initialized = false;
	bool has_current_step = false;
	double current_step_time = 0.0;
	double current_step_wall_time_s = 0.0;
	double previous_step_wall_time_s = 0.0;
	double current_step_cpu_util_percent = 0.0;
	double current_step_working_set_mb = 0.0;
	std::vector<double> samples_ms;
	std::vector<double> samples_working_set_delta_mb;
	std::vector<double> samples_private_bytes_delta_mb;
	std::vector<double> samples_dead_cell_count;
	#ifdef _WIN32
	static double filetime_to_seconds(const FILETIME &ft);
	#endif
	static double get_process_cpu_seconds();
	static double get_process_working_set_mb();
	static double get_process_private_bytes_mb();
	static double get_cpu_cores();
	static double normalize_step_time(double time);

	void open_if_needed();
	void flush_locked();
	void record(double time, double elapsed_ms, double working_set_delta_mb, double private_bytes_delta_mb, int dead_cell_count);
	void finalize();
	~FbaTimingLogger();
};

extern FbaTimingLogger fba_timing_logger;