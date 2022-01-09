#include "timer.hpp"

#include "type.hpp"

#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/mutex.h>

#include <tbb/parallel_for.h>
// #include <tbb/task_scheduler_init.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <algorithm>
#include <float.h>
#include <iostream>
#include <vector>

using namespace std;

namespace cpu_ti
{

	inline void get_numerical_error_vf_memory_pool(CCDdata &data_in)
	{
#ifdef GPUTI_USE_DOUBLE_PRECISION
		static const Scalar vffilter = 6.661338147750939e-15;
#else
		static const Scalar vffilter = 3.576279e-06;
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

	inline void get_numerical_error_ee_memory_pool(CCDdata &data_in)
	{
#ifdef GPUTI_USE_DOUBLE_PRECISION
		static const Scalar vffilter = 6.217248937900877e-15;
#else
		static const Scalar vffilter = 3.337861e-06;
#endif
		Scalar xmax = fabs(data_in.v0s[0]);
		Scalar ymax = fabs(data_in.v0s[1]);
		Scalar zmax = fabs(data_in.v0s[2]);

		xmax = std::max(xmax, fabs(data_in.v1s[0]));
		ymax = std::max(ymax, fabs(data_in.v1s[1]));
		zmax = std::max(zmax, fabs(data_in.v1s[2]));

		xmax = std::max(xmax, fabs(data_in.v2s[0]));
		ymax = std::max(ymax, fabs(data_in.v2s[1]));
		zmax = std::max(zmax, fabs(data_in.v2s[2]));

		xmax = std::max(xmax, fabs(data_in.v3s[0]));
		ymax = std::max(ymax, fabs(data_in.v3s[1]));
		zmax = std::max(zmax, fabs(data_in.v3s[2]));

		xmax = std::max(xmax, fabs(data_in.v0e[0]));
		ymax = std::max(ymax, fabs(data_in.v0e[1]));
		zmax = std::max(zmax, fabs(data_in.v0e[2]));

		xmax = std::max(xmax, fabs(data_in.v1e[0]));
		ymax = std::max(ymax, fabs(data_in.v1e[1]));
		zmax = std::max(zmax, fabs(data_in.v1e[2]));

		xmax = std::max(xmax, fabs(data_in.v2e[0]));
		ymax = std::max(ymax, fabs(data_in.v2e[1]));
		zmax = std::max(zmax, fabs(data_in.v2e[2]));

		xmax = std::max(xmax, fabs(data_in.v3e[0]));
		ymax = std::max(ymax, fabs(data_in.v3e[1]));
		zmax = std::max(zmax, fabs(data_in.v3e[2]));

		xmax = std::max(xmax, Scalar(1));
		ymax = std::max(ymax, Scalar(1));
		zmax = std::max(zmax, Scalar(1));

		data_in.err[0] = xmax * xmax * xmax * vffilter;
		data_in.err[1] = ymax * ymax * ymax * vffilter;
		data_in.err[2] = zmax * zmax * zmax * vffilter;
		return;
	}

	inline void compute_face_vertex_tolerance_memory_pool(CCDdata &data_in, const CCDConfig &config)
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

	inline void compute_edge_edge_tolerance_memory_pool(CCDdata &data_in, const CCDConfig &config)
	{
		Scalar p000[3], p001[3], p011[3], p010[3], p100[3], p101[3], p111[3], p110[3];
		for (int i = 0; i < 3; i++)
		{
			p000[i] = data_in.v0s[i] - data_in.v2s[i];
			p001[i] = data_in.v0s[i] - data_in.v3s[i];
			p011[i] = data_in.v1s[i] - data_in.v3s[i];
			p010[i] = data_in.v1s[i] - data_in.v2s[i];
			p100[i] = data_in.v0e[i] - data_in.v2e[i];
			p101[i] = data_in.v0e[i] - data_in.v3e[i];
			p111[i] = data_in.v1e[i] - data_in.v3e[i];
			p110[i] = data_in.v1e[i] - data_in.v2e[i];
		}
		Scalar dl = 0;
		for (int i = 0; i < 3; i++)
		{
			dl = std::max(dl, fabs(p100[i] - p000[i]));
			dl = std::max(dl, fabs(p101[i] - p001[i]));
			dl = std::max(dl, fabs(p111[i] - p011[i]));
			dl = std::max(dl, fabs(p110[i] - p010[i]));
		}
		dl *= 3;
		data_in.tol[0] = config.co_domain_tolerance / dl;

		dl = 0;
		for (int i = 0; i < 3; i++)
		{
			dl = std::max(dl, fabs(p010[i] - p000[i]));
			dl = std::max(dl, fabs(p110[i] - p100[i]));
			dl = std::max(dl, fabs(p111[i] - p101[i]));
			dl = std::max(dl, fabs(p011[i] - p001[i]));
		}
		dl *= 3;
		data_in.tol[1] = config.co_domain_tolerance / dl;

		dl = 0;
		for (int i = 0; i < 3; i++)
		{
			dl = std::max(dl, fabs(p001[i] - p000[i]));
			dl = std::max(dl, fabs(p101[i] - p100[i]));
			dl = std::max(dl, fabs(p111[i] - p110[i]));
			dl = std::max(dl, fabs(p011[i] - p010[i]));
		}
		dl *= 3;
		data_in.tol[2] = config.co_domain_tolerance / dl;
	}

	inline void compute_vf_tolerance_and_error_bound_memory_pool(CCDdata &data, const CCDConfig &config)
	{
		compute_face_vertex_tolerance_memory_pool(data, config);

#ifdef CALCULATE_ERROR_BOUND
		get_numerical_error_vf_memory_pool(data);
#endif
	}
	inline void compute_ee_tolerance_and_error_bound_memory_pool(CCDdata &data,
																 const CCDConfig &config)
	{
		compute_edge_edge_tolerance_memory_pool(data, config);

#ifdef CALCULATE_ERROR_BOUND
		get_numerical_error_ee_memory_pool(data);
#endif
	}

	inline std::array<Scalar, 8> function_ee(
		const Scalar &a0s,
		const Scalar &a1s,
		const Scalar &b0s,
		const Scalar &b1s,
		const Scalar &a0e,
		const Scalar &a1e,
		const Scalar &b0e,
		const Scalar &b1e,
		const std::array<Scalar, 8> &t,
		const std::array<Scalar, 8> &u,
		const std::array<Scalar, 8> &v)
	{
		std::array<Scalar, 8> rst;
		for (int i = 0; i < 8; i++)
		{
			Scalar edge0_vertex0 = (a0e - a0s) * t[i] + a0s;
			Scalar edge0_vertex1 = (a1e - a1s) * t[i] + a1s;
			Scalar edge1_vertex0 = (b0e - b0s) * t[i] + b0s;
			Scalar edge1_vertex1 = (b1e - b1s) * t[i] + b1s;

			Scalar edge0_vertex =
				(edge0_vertex1 - edge0_vertex0) * u[i]
				+ edge0_vertex0;
			Scalar edge1_vertex =
				(edge1_vertex1 - edge1_vertex0) * v[i]
				+ edge1_vertex0;
			rst[i] = edge0_vertex - edge1_vertex;
		}
		return rst;
	}

	inline std::array<Scalar, 8> function_vf(
		const Scalar &vs,
		const Scalar &t0s,
		const Scalar &t1s,
		const Scalar &t2s,
		const Scalar &ve,
		const Scalar &t0e,
		const Scalar &t1e,
		const Scalar &t2e,
		const std::array<Scalar, 8> &t,
		const std::array<Scalar, 8> &u,
		const std::array<Scalar, 8> &v)
	{
		std::array<Scalar, 8> rst;
		for (int i = 0; i < 8; i++)
		{
			Scalar ver = (ve - vs) * t[i] + vs;
			Scalar t0 = (t0e - t0s) * t[i] + t0s;
			Scalar t1 = (t1e - t1s) * t[i] + t1s;
			Scalar t2 = (t2e - t2s) * t[i] + t2s;
			Scalar pt = (t1 - t0) * u[i]
						+ (t2 - t0) * v[i] + t0;
			rst[i] = ver - pt;
		}
		return rst;
	}
	inline void convert_tuv_to_array(
		const MP_unit &unit,
		std::array<Scalar, 8> &t,
		std::array<Scalar, 8> &u,
		std::array<Scalar, 8> &v)
	{
		// t order: 0,0,0,0,1,1,1,1
		// u order: 0,0,1,1,0,0,1,1
		// v order: 0,1,0,1,0,1,0,1
		Scalar t0 = unit.itv[0].first,
			   t1 = unit.itv[0].second,
			   u0 = unit.itv[1].first,
			   u1 = unit.itv[1].second,
			   v0 = unit.itv[2].first,
			   v1 = unit.itv[2].second;
		t = {{t0, t0, t0, t0, t1, t1, t1, t1}};

		u = {{u0, u0, u1, u1, u0, u0, u1, u1}};

		v = {{v0, v1, v0, v1, v0, v1, v0, v1}};
	}

	inline bool evaluate_bbox_one_dimension_vector_return_tolerance(
		std::array<Scalar, 8> &t,
		std::array<Scalar, 8> &u,
		std::array<Scalar, 8> &v,
		const CCDdata &data,
		const int dimension,
		const bool is_edge_edge,
		bool &bbox_in_eps,
		Scalar &tol)
	{

		std::array<Scalar, 8> vs;
		bbox_in_eps = false;
		if (is_edge_edge) // ee
		{
			vs = function_ee(
				data.v0s[dimension], data.v1s[dimension], data.v2s[dimension], data.v3s[dimension],
				data.v0e[dimension], data.v1e[dimension], data.v2e[dimension], data.v3e[dimension],
				t, u, v);
		}
		else
		{
			vs = function_vf(
				data.v0s[dimension], data.v1s[dimension], data.v2s[dimension], data.v3s[dimension],
				data.v0e[dimension], data.v1e[dimension], data.v2e[dimension], data.v3e[dimension],
				t, u, v);
		}

		Scalar minv = vs[0], maxv = vs[0];

		for (int i = 1; i < 8; i++)
		{
			if (minv > vs[i])
			{
				minv = vs[i];
			}
			if (maxv < vs[i])
			{
				maxv = vs[i];
			}
		}
		tol = maxv - minv; // this is the real tolerance
		if (minv - data.ms > data.err[dimension] || maxv + data.ms < -data.err[dimension])
			return false;
		if (minv + data.ms >= -data.err[dimension] && maxv - data.ms <= data.err[dimension])
		{
			bbox_in_eps = true;
		}
		return true;
	}
	inline bool Origin_in_inclusion_function_memory_pool_avx2(

		const CCDdata &data,
		const bool is_edge_edge,
		const MP_unit &unit,
		Scalar &true_tol, bool &box_in_eps)
	{

		box_in_eps = false;
		std::array<Scalar, 8> t;
		std::array<Scalar, 8> u;
		std::array<Scalar, 8> v;
		Scalar tmp_tolerance;
		true_tol = 0;

		convert_tuv_to_array(unit, t, u, v);

		bool ck;
		bool box_in[3];
		for (int i = 0; i < 3; i++)
		{
			ck = evaluate_bbox_one_dimension_vector_return_tolerance(
				t, u, v, data, i, is_edge_edge, box_in[i],
				tmp_tolerance);
			true_tol = max(tmp_tolerance, true_tol);
			if (!ck)
				return false;
		}
		if (box_in[0] && box_in[1] && box_in[2])
		{
			box_in_eps = true;
		}
		return true;
	}

	inline bool Origin_in_inclusion_function_memory_pool(const CCDdata &data_in, const bool is_edge, const MP_unit &unit, Scalar &true_tol, bool &box_in)
	{
		box_in = true;
		true_tol = 0;

		BoxPrimitives bp;
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
						if (is_edge)
						{
							value = calculate_ee(data_in, bp);
						}
						else
						{
							value = calculate_vf(data_in, bp);
						}

						vmin = min(vmin, value);
						vmax = max(vmax, value);
					}
				}
			}

			true_tol = max(vmax - vmin, true_tol);

			if (vmin - data_in.ms > data_in.err[bp.dim] || vmax + data_in.ms < -data_in.err[bp.dim])
			{
				return false;
			}

			if (vmin + data_in.ms < -data_in.err[bp.dim] || vmax - data_in.ms > data_in.err[bp.dim])
			{
				box_in = false;
			}
		}
		return true;
	}

	inline void mutex_update_min(tbb::mutex &mutex, Scalar &value, const Scalar &compare)
	{
		mutex.lock();
		value = compare < value ? compare : value; // if compare is smaller, update it
		mutex.unlock();
	}

	// clarified in queue.h
	inline int split_dimension_memory_pool(const CCDdata &data, Scalar width[3])
	{
		int split = 0;
		Scalar res[3];
		res[0] = width[0] / data.tol[0];
		res[1] = width[1] / data.tol[1];
		res[2] = width[2] / data.tol[2];
		if (res[0] >= res[1] && res[0] >= res[2])
		{
			split = 0;
		}
		else if (res[1] >= res[0] && res[1] >= res[2])
		{
			split = 1;
		}
		else //if (res[2] >= res[1] && res[2] >= res[0])
		{
			split = 2;
		}

		return split;
	}

	inline bool bisect_vf_memory_pool(const MP_unit &unit, int split, const CCDConfig &config, std::vector<MP_unit> &out)
	{
		interval_pair halves(unit.itv[split]); // bisected

		if (halves.first.first >= halves.first.second)
		{
			return true;
		}
		if (halves.second.first >= halves.second.second)
		{
			return true;
		}

		// bisected[0] = unit;
		// bisected[1] = unit;

		out.emplace_back(unit);
		out.back().itv[split] = halves.first;

		if (split == 0)
		{
			if (halves.second.first <= config.toi)
			{
				out.emplace_back(unit);
				out.back().itv[split] = halves.second;
			}
		}
		else if (split == 1)
		{
			if (sum_no_larger_1(halves.second.first, unit.itv[2].first)) // check if u+v<=1
			{
				out.emplace_back(unit);
				out.back().itv[1] = halves.second;
			}
		}
		else if (split == 2)
		{
			if (sum_no_larger_1(halves.second.first, unit.itv[1].first)) // check if u+v<=1
			{
				out.emplace_back(unit);
				out.back().itv[2] = halves.second;
			}
		}

		return false;
	}

	inline bool bisect_ee_memory_pool(const MP_unit &unit, int split, const CCDConfig &config, std::vector<MP_unit> &out)
	{
		interval_pair halves(unit.itv[split]); // bisected

		if (halves.first.first >= halves.first.second)
		{
			return true;
		}
		if (halves.second.first >= halves.second.second)
		{
			return true;
		}

		out.emplace_back(unit);
		out.back().itv[split] = halves.first;

		if (split == 0) // split the time interval
		{
			if (halves.second.first <= config.toi)
			{
				out.emplace_back(unit);
				out.back().itv[split] = halves.second;
			}
		}
		else
		{
			out.emplace_back(unit);
			out.back().itv[split] = halves.second;
		}

		return false;
	}

	// input: "refine" is the number of iterations out side this function
	//        "qmutex" is the mutex for the queue
	void ccd_memory_pool_parallel( // parallel with different unit_id
		const bool is_edge, const std::vector<MP_unit> &vec_in, std::vector<MP_unit> &vec_out,
		const std::vector<CCDdata> &data_in, CCDConfig &config,
		int unit_id,
		tbb::mutex &mutex)
	{
		//TODO maybe record collision pairs

		const MP_unit temp_unit = vec_in[unit_id];
		const int box_id = temp_unit.query_id;

		const CCDdata data = data_in[box_id];
		const int query_size = data_in.size();

		Scalar widths[3];
		bool condition;
		const Scalar time_left = temp_unit.itv[0].first; // the time of this unit

		// if the time is larger than toi, return
		if (time_left >= config.toi)
		{
			return;
		}
		//ADD max checks!!!!!

		Scalar true_tol = 0;
		bool box_in;
#ifdef CPUTICCD_USE_AVX2
		const bool zero_in = Origin_in_inclusion_function_memory_pool_avx2(data, is_edge, temp_unit, true_tol, box_in);
#else
		const bool zero_in = Origin_in_inclusion_function_memory_pool(data, is_edge, temp_unit, true_tol, box_in);
#endif
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
				mutex_update_min(mutex, config.toi, time_left);
				return;
			}

			// Condition 2, the box is inside the epsilon box, have a root, return true;
			if (box_in)
			{
				mutex_update_min(mutex, config.toi, time_left);
				return;
			}
			// Condition 3, real tolerance is smaller than the input tolerance, return true

			condition = true_tol <= config.co_domain_tolerance;
			//TODO write into an "atomic" if necessary
			if (condition)
			{
				mutex_update_min(mutex, config.toi, time_left);
				return;
			}

			const int split = split_dimension_memory_pool(data, widths);

			// MP_unit bisected[2];
			// int valid_nbr;

			const bool sure_in = is_edge ? bisect_ee_memory_pool(temp_unit, split, config, vec_out) : bisect_vf_memory_pool(temp_unit, split, config, vec_out);

			if (sure_in) // in this case, the interval is too small that overflow happens. it should be rare to happen
			{
				mutex_update_min(mutex, config.toi, time_left);
				return;
			}
		}
	}

	inline CCDdata array_to_ccd(const std::array<std::array<Scalar, 3>, 8> &a)
	{
		CCDdata data;
#pragma unroll
		for (int i = 0; i < 3; i++)
		{
			data.v0s[i] = a[0][i];
			data.v1s[i] = a[1][i];
			data.v2s[i] = a[2][i];
			data.v3s[i] = a[3][i];
			data.v0e[i] = a[4][i];
			data.v1e[i] = a[5][i];
			data.v2e[i] = a[6][i];
			data.v3e[i] = a[7][i];
		}
		data.ms = 0;
		return data;
	}

	double ccd(const std::vector<std::array<std::array<Scalar, 3>, 8>> &V, bool is_edge, double max_t)
	{
		std::cout << "runnin CPU parallization" << std::endl;

		// Timer timer;
		// timer.start();

		int query_size = V.size();

		std::array<std::vector<MP_unit>, 2> units; // the input units and output units
		std::vector<CCDdata> data(query_size);     // input data list

		// std::cout << "input " << query_size << " " << is_edge << std::endl;

		int vec_in = 0; // the id of the input vec, switch between 0 and 1
		int vec_out = 1;

		units[vec_in].resize(query_size);

		// units[vec_out].resize(query_size * 2); // set the capicity of the output queue as twice large as the input
		CCDConfig config;
		config.err_in[0] = -1;             // the input error bound calculate from the AABB of the whole mesh
		config.co_domain_tolerance = 1e-6; // tolerance of the co-domain
		// config.max_itr = 1e6;              // the maximal nbr of iterations
		config.toi = max_t;
		tbb::mutex qmutex;

		// now initialize the CCDData and the units.

		tbb::parallel_for(tbb::blocked_range<int>(0, query_size), [&](tbb::blocked_range<int> r) {
			for (int i = r.begin(); i < r.end(); ++i)
			{
				units[vec_in][i].init(i);
				data[i] = array_to_ccd(V[i]);
				if (is_edge)
					compute_ee_tolerance_and_error_bound_memory_pool(data[i], config);
				else
					compute_vf_tolerance_and_error_bound_memory_pool(data[i], config);
			}
		});

		// std::cout<<"initialized"<<std::endl;

		while (1)
		{
			// check if the queue is empty
			int remain_unit_size = units[vec_in].size();
			if (remain_unit_size == 0)
				break;

			tbb::enumerable_thread_specific<std::vector<MP_unit>> storage;

			tbb::parallel_for(tbb::blocked_range<int>(0, remain_unit_size), [&](tbb::blocked_range<int> r) {
				std::vector<MP_unit> &local_storage = storage.local();

				// local_storage.reserve(remain_unit_size);
				for (int i = r.begin(); i < r.end(); ++i)
				{
					ccd_memory_pool_parallel(is_edge, units[vec_in], local_storage, data, config, i, qmutex);
				}
			});
			// timer.stop();
			// run_time = timer.getElapsedTimeInMicroSec();
			// std::cout << "timing temp " << run_time / 1000 << std::endl;

			// timer.start();
			std::vector<const std::vector<MP_unit> *> storages(storage.size());
			std::vector<int> offsets(storage.size());
			int index = 0;
			int size = 0;
			for (const auto &local_storage : storage)
			{
				offsets[index] = size;
				storages[index] = &local_storage;

				size += local_storage.size();
				++index;
			}

			units[vec_out].resize(size);
			tbb::parallel_for(0, int(storages.size()), [&](int i) {
				const auto *storage = storages[i];
				const int offset = offsets[i];
				for (int j = 0; j < storage->size(); ++j)
				{
					units[vec_out][offset + j] = storage->at(j);
				}
			});

			// timer.stop();
			// run_time = timer.getElapsedTimeInMicroSec();
			// std::cout << "done : " << units[vec_out].size() << std::endl;

			// timer.start();
			// tbb::parallel_sort(units[vec_out].begin(), units[vec_out].end(), [](const MP_unit &a, const MP_unit &b) {
			// 	if (a.query_id == b.query_id)
			// 		return a.itv[0].first < b.itv[0].first;
			// 	return a.query_id < b.query_id;
			// });
			// timer.stop();
			// run_time = timer.getElapsedTimeInMicroSec();
			// std::cout << "timing sort " << run_time / 1000 << std::endl;

			vec_out = vec_in;
			vec_in = !bool(vec_out); // switch in and out, now the out contains the old queue

			// break;
		}

		// timer.stop();
		// run_time = timer.getElapsedTimeInMicroSec();

		// std::cout << config.toi << std::endl;

		return config.toi;
	}
} // namespace cpu_ti
