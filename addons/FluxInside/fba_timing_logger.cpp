#include "fba_timing_logger.h"
#include "../../modules/PhysiCell_settings.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

double FbaTimingLogger::filetime_to_seconds(const FILETIME &ft)
{
	ULARGE_INTEGER value;
	value.LowPart = ft.dwLowDateTime;
	value.HighPart = ft.dwHighDateTime;
	return static_cast<double>(value.QuadPart) * 1e-7;
}

double FbaTimingLogger::get_process_cpu_seconds()
{
#ifdef _WIN32
	FILETIME creation_time;
	FILETIME exit_time;
	FILETIME kernel_time;
	FILETIME user_time;
	if (!GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time, &kernel_time, &user_time))
		throw std::runtime_error("GetProcessTimes failed while collecting process CPU usage");

	return filetime_to_seconds(kernel_time) + filetime_to_seconds(user_time);
#else
	throw std::runtime_error("Process CPU usage is only implemented on Windows");
#endif
}

double FbaTimingLogger::get_process_working_set_mb()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS_EX counters;
	if (!GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS *>(&counters), sizeof(counters)))
		throw std::runtime_error("GetProcessMemoryInfo failed while collecting process memory usage");

	return static_cast<double>(counters.WorkingSetSize) / (1024.0 * 1024.0);
#else
	throw std::runtime_error("Process memory usage is only implemented on Windows");
#endif
}

double FbaTimingLogger::get_process_private_bytes_mb()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS_EX counters;
	if (!GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS *>(&counters), sizeof(counters)))
		throw std::runtime_error("GetProcessMemoryInfo failed while collecting process private bytes");

	return static_cast<double>(counters.PrivateUsage) / (1024.0 * 1024.0);
#else
	throw std::runtime_error("Process private bytes are only implemented on Windows");
#endif
}

double FbaTimingLogger::get_cpu_cores()
{
#ifdef _WIN32
	SYSTEM_INFO system_info;
	GetSystemInfo(&system_info);
	return static_cast<double>(std::max<DWORD>(1, system_info.dwNumberOfProcessors));
#else
	return 1.0;
#endif
}

double FbaTimingLogger::normalize_step_time(double time)
{
	return std::round(time * 1e9) / 1e9;
}

void FbaTimingLogger::open_if_needed()
{
	if (initialized)
		return;

	//get the output directory from physicell configuration or use default
	stream.open(PhysiCell::PhysiCell_settings.folder + "/fba_linear_optimization_times.csv", std::ios::out | std::ios::trunc);
	if (!stream)
		throw std::runtime_error("Failed to open output/fba_linear_optimization_times.csv for writing");

	start_wall_clock = std::chrono::steady_clock::now();
	previous_wall_clock = start_wall_clock;
	stream << "time_step,physicell_time,current_simulated_time,wall_time_s,interval_wall_time_s,process_cpu_util_percent,process_working_set_mb,working_set_delta_mb_avg,private_bytes_delta_mb_avg,dead_cell_count,total_ms,min_ms,max_ms,avg_ms,median_ms,std_err_ms,n_samples\n";
	initialized = true;
}

void FbaTimingLogger::flush_locked()
{
	if (!has_current_step || samples_ms.empty())
		return;

	std::vector<double> sorted = samples_ms;
	std::sort(sorted.begin(), sorted.end());

	const auto minmax = std::minmax_element(samples_ms.begin(), samples_ms.end());
	const double sum = std::accumulate(samples_ms.begin(), samples_ms.end(), 0.0);
	const double avg = sum / static_cast<double>(samples_ms.size());

	double median = 0.0;
	const size_t mid = sorted.size() / 2;
	if (sorted.size() % 2 == 0)
		median = 0.5 * (sorted[mid - 1] + sorted[mid]);
	else
		median = sorted[mid];

	double std_err = 0.0;
	if (samples_ms.size() > 1)
	{
		double variance = 0.0;
		for (double sample : samples_ms)
		{
			const double delta = sample - avg;
			variance += delta * delta;
		}
		variance /= static_cast<double>(samples_ms.size() - 1);
		std_err = std::sqrt(variance) / std::sqrt(static_cast<double>(samples_ms.size()));
	}

	double process_cpu_util_percent = current_step_cpu_util_percent;
	double process_working_set_mb = current_step_working_set_mb;
	auto average_or_zero = [](const std::vector<double> &values) {
		if (values.empty())
			return 0.0;
		return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
	};


	

	stream << current_step_time << ','
		   << current_step_time << ','
		   << current_step_time << ','
		   << current_step_wall_time_s << ','
		   << (current_step_wall_time_s - previous_step_wall_time_s) << ','
		   << process_cpu_util_percent << ','
		   << process_working_set_mb << ','
		   << average_or_zero(samples_working_set_delta_mb) << ','
		   << average_or_zero(samples_private_bytes_delta_mb) << ','
		   << static_cast<int>(samples_dead_cell_count.back()) << ','
		   << sum << ','
		   << *minmax.first << ','
		   << *minmax.second << ','
		   << avg << ','
		   << median << ','
		   << std_err << ','
		   << samples_ms.size() << ','
		   <<'\n';

	stream.flush();
	samples_ms.clear();
	samples_working_set_delta_mb.clear();
	samples_private_bytes_delta_mb.clear();
	samples_dead_cell_count.clear();
}

void FbaTimingLogger::record(double time, double elapsed_ms, double working_set_delta_mb, double private_bytes_delta_mb, int dead_cell_count)
{
	std::lock_guard<std::mutex> lock(mutex);
	open_if_needed();

	const double step_time = normalize_step_time(time);
	const double wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_wall_clock).count();
	const double process_working_set_mb = get_process_working_set_mb();
	const double process_private_bytes_mb = get_process_private_bytes_mb();
	double process_cpu_seconds = get_process_cpu_seconds();
	double wall_interval_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - previous_wall_clock).count();
	double cpu_interval_s = 0.0;
	if (has_previous_process_sample)
		cpu_interval_s = process_cpu_seconds - previous_process_cpu_s;
	double cpu_util_percent = 0.0;
	if (wall_interval_s > 0.0)
		cpu_util_percent = (cpu_interval_s / (wall_interval_s * get_cpu_cores())) * 100.0;
	if (!has_current_step)
	{
		current_step_time = step_time;
		current_step_wall_time_s = wall_time_s;
		current_step_cpu_util_percent = cpu_util_percent;
		current_step_working_set_mb = process_working_set_mb;
		has_current_step = true;
	}
	else if (std::abs(step_time - current_step_time) > 1e-9)
	{
		flush_locked();
		previous_step_wall_time_s = current_step_wall_time_s;
		current_step_time = step_time;
		current_step_wall_time_s = wall_time_s;
		current_step_cpu_util_percent = cpu_util_percent;
		current_step_working_set_mb = process_working_set_mb;
	}
	else
	{
		current_step_wall_time_s = wall_time_s;
		current_step_cpu_util_percent = cpu_util_percent;
		current_step_working_set_mb = process_working_set_mb;
	}
	previous_wall_clock = std::chrono::steady_clock::now();
	previous_process_cpu_s = process_cpu_seconds;
	has_previous_process_sample = true;

	samples_ms.push_back(elapsed_ms);
	samples_working_set_delta_mb.push_back(working_set_delta_mb);
	samples_private_bytes_delta_mb.push_back(private_bytes_delta_mb);
	samples_dead_cell_count.push_back(dead_cell_count);
}

void FbaTimingLogger::finalize()
{
	std::lock_guard<std::mutex> lock(mutex);
	if (initialized)
	{
		flush_locked();
		stream.flush();
	}
}

FbaTimingLogger::~FbaTimingLogger()
{
	try
	{
		finalize();
	}
	catch (...)
	{
	}
}

FbaTimingLogger fba_timing_logger;