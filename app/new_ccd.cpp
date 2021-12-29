#include "timer.hpp"
#include <algorithm>
#include <cputi/root_finder.h>
#include <float.h>
#include <iostream>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>

#include <vector>
using namespace std;
void get_numerical_error_vf_memory_pool(CCDdata &data_in)
{
	Scalar vffilter;

#ifdef GPUTI_USE_DOUBLE_PRECISION
	vffilter = 6.661338147750939e-15;
#else
	vffilter = 3.576279e-06;
#endif
	Scalar xmax = fabs(data_in.v0s[0]);
	Scalar ymax = fabs(data_in.v0s[1]);
	Scalar zmax = fabs(data_in.v0s[2]);

	xmax = max(xmax, fabs(data_in.v1s[0]));
	ymax = max(ymax, fabs(data_in.v1s[1]));
	zmax = max(zmax, fabs(data_in.v1s[2]));

	xmax = max(xmax, fabs(data_in.v2s[0]));
	ymax = max(ymax, fabs(data_in.v2s[1]));
	zmax = max(zmax, fabs(data_in.v2s[2]));

	xmax = max(xmax, fabs(data_in.v3s[0]));
	ymax = max(ymax, fabs(data_in.v3s[1]));
	zmax = max(zmax, fabs(data_in.v3s[2]));

	xmax = max(xmax, fabs(data_in.v0e[0]));
	ymax = max(ymax, fabs(data_in.v0e[1]));
	zmax = max(zmax, fabs(data_in.v0e[2]));

	xmax = max(xmax, fabs(data_in.v1e[0]));
	ymax = max(ymax, fabs(data_in.v1e[1]));
	zmax = max(zmax, fabs(data_in.v1e[2]));

	xmax = max(xmax, fabs(data_in.v2e[0]));
	ymax = max(ymax, fabs(data_in.v2e[1]));
	zmax = max(zmax, fabs(data_in.v2e[2]));

	xmax = max(xmax, fabs(data_in.v3e[0]));
	ymax = max(ymax, fabs(data_in.v3e[1]));
	zmax = max(zmax, fabs(data_in.v3e[2]));

	xmax = max(xmax, Scalar(1));
	ymax = max(ymax, Scalar(1));
	zmax = max(zmax, Scalar(1));

	data_in.err[0] = xmax * xmax * xmax * vffilter;
	data_in.err[1] = ymax * ymax * ymax * vffilter;
	data_in.err[2] = zmax * zmax * zmax * vffilter;
	return;
}
void compute_face_vertex_tolerance_memory_pool(CCDdata &data_in,
											   const CCDConfig &config)
{
	Scalar p000[3], p001[3], p011[3], p010[3], p100[3], p101[3], p111[3], p110[3];
	for (int i = 0; i < 3; i++)
	{
		p000[i] = data_in.v0s[i] - data_in.v1s[i];
		p001[i] = data_in.v0s[i] - data_in.v3s[i];
		p011[i] =
			data_in.v0s[i] - (data_in.v2s[i] + data_in.v3s[i] - data_in.v1s[i]);
		p010[i] = data_in.v0s[i] - data_in.v2s[i];
		p100[i] = data_in.v0e[i] - data_in.v1e[i];
		p101[i] = data_in.v0e[i] - data_in.v3e[i];
		p111[i] =
			data_in.v0e[i] - (data_in.v2e[i] + data_in.v3e[i] - data_in.v1e[i]);
		p110[i] = data_in.v0e[i] - data_in.v2e[i];
	}
	Scalar dl = 0;
	for (int i = 0; i < 3; i++)
	{
		dl = max(dl, fabs(p100[i] - p000[i]));
		dl = max(dl, fabs(p101[i] - p001[i]));
		dl = max(dl, fabs(p111[i] - p011[i]));
		dl = max(dl, fabs(p110[i] - p010[i]));
	}
	dl *= 3;
	data_in.tol[0] = config.co_domain_tolerance / dl;

	dl = 0;
	for (int i = 0; i < 3; i++)
	{
		dl = max(dl, fabs(p010[i] - p000[i]));
		dl = max(dl, fabs(p110[i] - p100[i]));
		dl = max(dl, fabs(p111[i] - p101[i]));
		dl = max(dl, fabs(p011[i] - p001[i]));
	}
	dl *= 3;
	data_in.tol[1] = config.co_domain_tolerance / dl;

	dl = 0;
	for (int i = 0; i < 3; i++)
	{
		dl = max(dl, fabs(p001[i] - p000[i]));
		dl = max(dl, fabs(p101[i] - p100[i]));
		dl = max(dl, fabs(p111[i] - p110[i]));
		dl = max(dl, fabs(p011[i] - p010[i]));
	}
	dl *= 3;
	data_in.tol[2] = config.co_domain_tolerance / dl;
}

void compute_vf_tolerance_and_error_bound_memory_pool(CCDdata &data,
													  const CCDConfig &config)
{
	compute_face_vertex_tolerance_memory_pool(data, config);

	data.last_round_has_root = false;
	data.last_round_has_root_record = false;
	data.sure_have_root = false;
	data.nbr_pushed = 1; // initially the number of pushed element is 1
#ifdef CALCULATE_ERROR_BOUND
	get_numerical_error_vf_memory_pool(data);
#endif
}
// void BoxPrimatives::calculate_tuv(const MP_unit &unit)
// {
// 	if (b[0] == 0)
// 	{ // t0
// 		t = unit.itv[0].first;
// 	}
// 	else
// 	{ // t1
// 		t = unit.itv[0].second;
// 	}

// 	if (b[1] == 0)
// 	{ // u0
// 		u = unit.itv[1].first;
// 	}
// 	else
// 	{ // u1
// 		u = unit.itv[1].second;
// 	}

// 	if (b[2] == 0)
// 	{ // v0
// 		v = unit.itv[2].first;
// 	}
// 	else
// 	{ // v1
// 		v = unit.itv[2].second;
// 	}
// }
inline bool Origin_in_vf_inclusion_function_memory_pool(const CCDdata &data_in,
														MP_unit &unit)
{
	BoxPrimatives bp;
	Scalar vmin = SCALAR_LIMIT;
	Scalar vmax = -SCALAR_LIMIT;
	Scalar value;

	for (bp.dim = 0; bp.dim < 3; bp.dim++)
	{
		vmin = SCALAR_LIMIT;
		vmax = -SCALAR_LIMIT;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					bp.b[0] = i;
					bp.b[1] = j;
					bp.b[2] = k; // 100
					bp.calculate_tuv(unit);
					value = calculate_vf(data_in, bp);
					vmin = min(vmin, value);
					vmax = max(vmax, value);
				}
			}
		}

		// get the min and max in one dimension
		unit.true_tol =
			max(unit.true_tol, vmax - vmin); // this is the real tolerance

		if (vmin > data_in.err[bp.dim] || vmax < -data_in.err[bp.dim])
		{
			return false;
		}

		if (vmin < -data_in.err[bp.dim] || vmax > data_in.err[bp.dim])
		{
			unit.box_in = false;
		}
	}
	return true;
}
int mutex_add_1(tbb::mutex &mutex, int &a)
{
	int result;
	mutex.lock();
	a++;
	result = a;
	mutex.unlock();
	return result;
}

void mutex_equal(tbb::mutex &mutex, int &a, const int &b)
{
	mutex.lock();
	a = b;
	mutex.unlock();
}
int mutex_add(tbb::mutex &mutex, int &a, const int &value)
{
	int result;
	mutex.lock();
	a += value;
	result = a;
	mutex.unlock();
	return result;
}

int mutex_increase(tbb::mutex &mutex, int &a, const int size)
{
	int result;
	mutex.lock();
	a = a + 1 >= size ? 0 : a + 1;
	result = a;
	mutex.unlock();
	return result;
}

inline void split_dimension_memory_pool(const CCDdata &data, Scalar width[3],
										int &split)
{ // clarified in queue.h
	Scalar res[3];
	res[0] = width[0] / data.tol[0];
	res[1] = width[1] / data.tol[1];
	res[2] = width[2] / data.tol[2];
	if (res[0] >= res[1] && res[0] >= res[2])
	{
		split = 0;
	}
	if (res[1] >= res[0] && res[1] >= res[2])
	{
		split = 1;
	}
	if (res[2] >= res[1] && res[2] >= res[0])
	{
		split = 2;
	}
}

inline void bisect_vf_memory_pool(MP_unit &unit, int split,
								  const CCDConfig &config, MP_unit bisected[2],
								  int &valid_nbr)
{

	interval_pair halves(unit.itv[split]); // bisected

	if (halves.first.first >= halves.first.second)
	{
		valid_nbr = 0;
		return;
	}
	if (halves.second.first >= halves.second.second)
	{
		valid_nbr = 0;
		return;
	}

	bisected[0] = unit;
	bisected[1] = unit;
	valid_nbr = 1;
	bisected[0].itv[split] = halves.first;

	if (split == 0)
	{
		if (config.max_t != 1)
		{
			if (halves.second.first <= config.max_t)
			{
				bisected[1].itv[split] = halves.second;
				valid_nbr = 2;
			}
		}
		else
		{
			bisected[1].itv[split] = halves.second;
			valid_nbr = 2;
		}
	}

	if (split == 1)
	{
		if (sum_no_larger_1(halves.second.first,
							bisected[1].itv[2].first)) // check if u+v<=1
		{

			bisected[1].itv[1] = halves.second;
			valid_nbr = 2;
		}
	}
	if (split == 2)
	{
		if (sum_no_larger_1(halves.second.first,
							bisected[1].itv[1].first)) // check if u+v<=1
		{

			bisected[1].itv[2] = halves.second;
			valid_nbr = 2;
		}
	}
}
// input: "refine" is the number of iterations out side this function
//        "qmutex" is the mutex for the queue
void vf_ccd_memory_pool_parallel( // parallel with different unit_id
	const std::vector<MP_unit> &vec_in, std::vector<MP_unit> &vec_out,
	const std::vector<CCDdata> &data_in, CCDConfig &config,
	std::vector<int> &sure_have_root, std::vector<int> &results, int unit_id, std::vector<tbb::mutex> &mutex,
	tbb::mutex &qmutex)
{
	const int box_id = vec_in[unit_id].query_id;
	const int query_size = data_in.size();
  CCDdata data=data_in[box_id];
	if (sure_have_root[box_id] > 0)
	{
		return;
	}

	Scalar widths[3];
	bool condition;
	int split;

	MP_unit temp_unit = vec_in[unit_id];
	const bool zero_in = Origin_in_vf_inclusion_function_memory_pool(data, temp_unit);
	// mutex_add(mutex[box_id + query_size], data[box_id].nbr_pushed, -1); // queue size-=1
  
	if (zero_in)
	{
		// std::cout<<"###have root"<<std::endl;
		widths[0] = temp_unit.itv[0].second - temp_unit.itv[0].first;
		widths[1] = temp_unit.itv[1].second - temp_unit.itv[1].first;
		widths[2] = temp_unit.itv[2].second - temp_unit.itv[2].first;

		// Condition 1
		condition = widths[0] <= data.tol[0] && widths[1] <= data.tol[1] && widths[2] <= data.tol[2];
		if (condition)
		{
			sure_have_root[box_id] = 1;
			return;
		}
		// Condition 2, the box is inside the epsilon box, have a root, return true;
		condition = temp_unit.box_in;
		if (condition)
		{
			sure_have_root[box_id] = 1;
			return;
		}
		// Condition 3, real tolerance is smaller than the input tolerance, return true
		condition = temp_unit.true_tol <= config.co_domain_tolerance;
		if (condition)
		{
			sure_have_root[box_id] = 1;
			return;
		}

		split_dimension_memory_pool(data, widths, split);

		MP_unit bisected[2];
		int valid_nbr;

		bisect_vf_memory_pool(temp_unit, split, config, bisected, valid_nbr);

		bisected[0].query_id = box_id;
		bisected[1].query_id = box_id;
		if (valid_nbr == 0) // in this case, the interval is too small that overflow happens. it should be rare to happen
		{
			sure_have_root[box_id] = 1;
			return;
		}
		else if (valid_nbr == 1)
		{
			vec_out.push_back(bisected[0]);
		}
		else if (valid_nbr == 2)
		{
			vec_out.push_back(bisected[0]);
			vec_out.push_back(bisected[1]);
		}

		// if (data[box_id].nbr_pushed > UNIT_SIZE)
		// { // if heap overflow happens, we
		//   // regard it as having root
		// 	sure_have_root = 1;
		// 	// mutex_equal(mutex[box_id], data[box_id].sure_have_root,
		// 	//             1); // TODO remove this to make sure the queue is usable
		// 	// //                 // for next stage
		// }
	}
	// if (!no_need_check)
	// if (sure_have_root)
	// {
	// 	// data[box_id].sure_have_root=1;
	// }
}

void memory_pool_ccd_run(
	const std::vector<std::array<std::array<Scalar, 3>, 8>> &V, bool is_edge,
	std::vector<int> &result_list, double &run_time,
	std::vector<Scalar> &time_impact, int parallel_nbr)
{
	std::cout << "runnin CPU parallization" << std::endl;

	ccd::Timer timer;

	int query_size = V.size();
	result_list.resize(V.size());

	std::array<std::vector<MP_unit>, 2> units; // the input units and output units
	std::vector<CCDdata> data;                 // input data list
	std::vector<int> results;

	data.resize(query_size);
	results.resize(query_size);
	int vec_in = 0; // the id of the input vec, switch between 0 and 1
	int vec_out = 1;
	units[vec_in].resize(query_size * 2);
	units[vec_out].resize(
		query_size * 2); // set the capicity of the output queue as twice large as the input
	CCDConfig config;
	config.err_in[0] =
		-1;                            // the input error bound calculate from the AABB of the whole mesh
	config.co_domain_tolerance = 1e-6; // tolerance of the co-domain
	config.max_t = 1;                  // the upper bound of the time interval
	config.max_itr = 1e6;              // the maximal nbr of iterations
	config.mp_end = query_size - 1;    // the initialized trunk is from 0 to nbr-1;
	config.mp_status = true;           // in the begining, start < end
	config.not_empty = 0;
	std::vector<tbb::mutex> mutexes(query_size * 2);
	tbb::mutex qmutex;
	int nbr_itr = 0;
	tbb::enumerable_thread_specific<std::vector<MP_unit>> storage;

	// now initialize the CCDData and the units.
	tbb::parallel_for(tbb::blocked_range<int>(0, query_size),
					  [&](tbb::blocked_range<int> r) // initialize the units
					  {
						  for (int i = r.begin(); i < r.end(); ++i)
						  {
							  units[vec_in][i].init(i);
							  data[i] = array_to_ccd(V[i], is_edge);
							  compute_vf_tolerance_and_error_bound_memory_pool(
								  data[i], config);
						  }
					  });
	// std::cout<<"initialized"<<std::endl;
	int max_queue_size = 0;
	std::vector<int> must_skip(data.size(), 0);
	timer.start();
	while (1)
	{
		// check if the queue is empty
		int remain_unit_size = config.mp_end + 1; // the size of the input queue
		max_queue_size = std::max(max_queue_size, remain_unit_size);
		// std::cout<<"remain_unit_size "<<remain_unit_size<<std::endl;
		// std::cout<<"## remain_unit_size "<<remain_unit_size<<" itr "<<nbr_itr<<"
		// while checking vec_in "<<vec_in<<std::endl;
		if (remain_unit_size == 0)
		{
			break;
		}
		// config.mp_end = -1; // clear the output queue

		tbb::parallel_for(tbb::blocked_range<int>(0, remain_unit_size),
						  [&](tbb::blocked_range<int> r) {
							  std::vector<MP_unit> &local_storage = storage.local();
							  local_storage.reserve(remain_unit_size);
							  for (int i = r.begin(); i < r.end(); ++i)
							  {
								  vf_ccd_memory_pool_parallel(units[vec_in], local_storage, data, config,
															  must_skip, result_list, i,
															  mutexes, qmutex);
							  }
						  });
		// read from in, write to out
		// std::cout<<"itr "<<nbr_itr<<" finished"<<std::endl;
		timer.stop();
		run_time = timer.getElapsedTimeInMicroSec();
		std::cout << "timing temp " << run_time / 1000 << std::endl;

		timer.start();
		int size = 0;
		units[vec_out].clear();
		for (const auto &local_storage : storage)
		{
			size += local_storage.size();
			units[vec_out].insert(units[vec_out].end(), local_storage.begin(), local_storage.end());
		}

		timer.stop();
		run_time = timer.getElapsedTimeInMicroSec();
		std::cout << "timing merge " << run_time / 1000 << " size=" << size << "/" << units[vec_out].size() << std::endl;

		timer.start();

		// tbb::parallel_sort
		// std::sort
		tbb::parallel_sort(units[vec_out].begin(), units[vec_out].end(), [](const MP_unit &a, const MP_unit &b) {
			if (a.query_id == b.query_id)
				return a.itv[2].first < b.itv[2].first;
			return a.query_id < b.query_id;
		});
		timer.stop();
		run_time = timer.getElapsedTimeInMicroSec();
		std::cout << "timing sort " << run_time / 1000 << std::endl;
		exit(0);

		vec_out = vec_in;
		vec_in = !bool(vec_out); // switch in and out, now the out contains the old queue

		nbr_itr++;
		if (nbr_itr > UNIT_SIZE + 10)
		{
			// std::cout<<"WARNING the queue should be empty, this should not
			// happen!!!!"<<std::endl;
			break;
		}
	}

	// get the results
	tbb::parallel_for(
		tbb::blocked_range<int>(0, query_size), [&](tbb::blocked_range<int> r) {
			for (int i = r.begin(); i < r.end(); ++i)
			{
				result_list[i] =
					data[i].sure_have_root; // TODO deal with the remained blocks in
											// the current unit list
			}
		});
	// if (units[vec_out].size() > 0) {
	//   tbb::mutex final_mutex;
	//   tbb::parallel_for(
	//       tbb::blocked_range<int>(0, units[vec_out].size()),
	//       [&](tbb::blocked_range<int> r) {
	//         for (int i = r.begin(); i < r.end(); ++i) {
	//           if()
	//           result_list[i] =
	//               data[i].sure_have_root; // TODO deal with the remained blocks
	//               in
	//                                       // the current unit list
	//         }
	//       });
	// }

	timer.stop();
	run_time = timer.getElapsedTimeInMicroSec();
	int trues = 0;
	for (int i = 0; i < result_list.size(); i++)
	{
		if (result_list[i])
		{
			trues++;
		}
	}
	std::cout << "maximal queue size " << max_queue_size << std::endl;
	std::cout << "THE number of returned trues " << trues << std::endl
			  << std::endl;
	;
	return;
}
// input: "refine" is the number of iterations out side this function
void vf_ccd_memory_pool_parallel_narrow_phase( // parallel with different
											   // unit_id
	std::array<tbb::concurrent_vector<MP_unit>, 2> &units,
	std::vector<CCDdata> &data, CCDConfig &config, std::vector<int> &results,
	int vec_in, int unit_id, std::vector<tbb::mutex> &mutex, int refine)
{
	// bool no_need_check = false;
	// int vec_out = !bool(vec_in);
	// int box_id = units[vec_in][unit_id].query_id;
	// int query_size = data.size();
	// if (data[box_id].sure_have_root >
	//     0) { // if it is sure that have root, then no need to check
	//   no_need_check = true;
	// }
	// Scalar widths[3];
	// bool condition;
	// int split;
	// // mind here we do not need to set last_round_has_root status
	// // because we can get the info by checking the units

	// if (!no_need_check) { // if need check, and the for loop
	//                       // is not broken, do the check
	//   bool zero_in = Origin_in_vf_inclusion_function_memory_pool(
	//       data[box_id], units[vec_in][unit_id]);
	//   mutex_add(mutex[box_id + query_size], data[box_id].nbr_pushed,
	//             -1); // queue size-=1
	//   if (zero_in) {
	//     // std::cout<<"###have root"<<std::endl;
	//     widths[0] = units[vec_in][unit_id].itv[0].second -
	//                 units[vec_in][unit_id].itv[0].first;
	//     widths[1] = units[vec_in][unit_id].itv[1].second -
	//                 units[vec_in][unit_id].itv[1].first;
	//     widths[2] = units[vec_in][unit_id].itv[2].second -
	//                 units[vec_in][unit_id].itv[2].first;

	//     // Condition 1
	//     condition = widths[0] <= data[box_id].tol[0] &&
	//                 widths[1] <= data[box_id].tol[1] &&
	//                 widths[2] <= data[box_id].tol[2];
	//     if (condition) {

	//       mutex_equal(mutex[box_id], data[box_id].sure_have_root, 1);
	//     }
	//     // Condition 2, the box is inside the epsilon box, have a root, return
	//     // true;
	//     condition = units[vec_in][unit_id].box_in;
	//     if (condition) {
	//       mutex_equal(mutex[box_id], data[box_id].sure_have_root, 1);
	//     }

	//     // Condition 3, real tolerance is smaller than the input tolerance,
	//     // return true
	//     condition = units[vec_in][unit_id].true_tol <=
	//     config.co_domain_tolerance; if (condition) {
	//       mutex_equal(mutex[box_id], data[box_id].sure_have_root, 1);
	//     }
	//     split_dimension_memory_pool(data[box_id], widths, split);
	//     MP_unit bisected[2];
	//     int valid_nbr;
	//     bisect_vf_memory_pool(units[vec_in][unit_id], split, config, bisected,
	//                           valid_nbr);

	//     bisected[0].query_id = box_id;
	//     bisected[1].query_id = box_id;
	//     if (valid_nbr == 0) { // in this case, the interval is too small that
	//                           // overflow happens. it should be rare to happen
	//       mutex_equal(mutex[box_id], data[box_id].sure_have_root, 1);
	//     }
	//     if (valid_nbr == 1) {
	//       units[vec_out].push_back(bisected[0]);
	//       mutex_add(mutex[box_id + query_size], data[box_id].nbr_pushed,
	//                 1); // substract add 1
	//       // std::cout<<"###pushed 1, vec_out "<<vec_out<<" size
	//       // "<<units[vec_out].size()<<std::endl;
	//     }
	//     if (valid_nbr == 2) {
	//       units[vec_out].push_back(bisected[0]);
	//       units[vec_out].push_back(bisected[1]);
	//       mutex_add(mutex[box_id + query_size], data[box_id].nbr_pushed,
	//                 2); // substract add 1
	//       // std::cout<<"###pushed 2, vec_out "<<vec_out<<" size
	//       // "<<units[vec_out].size()<<std::endl;
	//     }
	//     if (data[box_id].nbr_pushed > UNIT_SIZE) { // if heap overflow happens,
	//     we
	//                                                // regard it as having root
	//       mutex_equal(mutex[box_id], data[box_id].sure_have_root,
	//                   1); // TODO remove this to make sure the queue is usable
	//                       // for next stage
	//     }
	//   }
	// }
}