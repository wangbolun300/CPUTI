#pragma once

#include <array>
#include <assert.h>
#include <limits>
#include <utility>
#include <atomic>
#include <vector>
#include <float.h>
#include <climits>

// #define GPUTI_USE_DOUBLE_PRECISION
// #define CALCULATE_ERROR_BOUND

namespace cpu_ti
{

#ifdef GPUTI_USE_DOUBLE_PRECISION
	typedef double Scalar;
#define SCALAR_LIMIT DBL_MAX;
#else
	typedef float Scalar;
#define SCALAR_LIMIT INT_MAX;
#endif

	// TODO next when spliting time intervals, check if overlaps the current toi,
	// then decide if we push it into the heap the reason of considerting it is that
	// the limited heap size. token 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

	class Singleinterval
	{
	public:
		inline Singleinterval(){};
		inline Singleinterval(const Scalar &f, const Scalar &s);

		Scalar first;
		Scalar second;

		Singleinterval &operator=(const Singleinterval &x)
		{
			if (this == &x)
				return *this;
			first = x.first;
			second = x.second;
			return *this;
		}
	};

	class MP_unit
	{
	public:
		inline MP_unit(){};

		Singleinterval itv[3];
		int query_id;

		void init(int i)
		{
			itv[0].first = 0;
			itv[0].second = 1;
			itv[1].first = 0;
			itv[1].second = 1;
			itv[2].first = 0;
			itv[2].second = 1;
			query_id = i;
		}

		MP_unit &operator=(const MP_unit &x)
		{
			if (this == &x)
				return *this;
			itv[0] = x.itv[0];
			itv[1] = x.itv[1];
			itv[2] = x.itv[2];
			query_id = x.query_id;
			return *this;
		}
	};

	// the initialized error input, solve tolerance, time interval upper bound, etc.

	class CCDConfig
	{
	public:
		inline CCDConfig(){};
		Scalar co_domain_tolerance; // tolerance of the co-domain
		Scalar ms;                  // minimum separation
		// int max_itr;                // the maximal nbr of iterations, TODO not there yet
		Scalar toi;
	};

	class BoxPrimitives
	{
	public:
		bool b[3];
		int dim;
		Scalar t;
		Scalar u;
		Scalar v;
		inline BoxPrimitives(){};

		inline void calculate_tuv(const MP_unit &unit)
		{
			if (b[0] == 0)
			{ // t0
				t = unit.itv[0].first;
			}
			else
			{ // t1
				t = unit.itv[0].second;
			}

			if (b[1] == 0)
			{ // u0
				u = unit.itv[1].first;
			}
			else
			{ // u1
				u = unit.itv[1].second;
			}

			if (b[2] == 0)
			{ // v0
				v = unit.itv[2].first;
			}
			else
			{ // v1
				v = unit.itv[2].second;
			}
		};
	};

	class CCDdata
	{
	public:
		inline CCDdata(){};
		// CCDdata(const std::array<std::array<Scalar,3>,8>&input);
		Scalar v0s[3];
		Scalar v1s[3];
		Scalar v2s[3];
		Scalar v3s[3];
		Scalar v0e[3];
		Scalar v1e[3];
		Scalar v2e[3];
		Scalar v3e[3];

#ifdef CALCULATE_ERROR_BOUND
		Scalar err[3]; // error bound of each query, calculated from each scene
#endif
		Scalar tol[3]; // domain tolerance that helps to decide which dimension to split

		CCDdata &operator=(const CCDdata &x)
		{
			if (this == &x)
				return *this;
			for (int i = 0; i < 3; i++)
			{
				v0s[i] = x.v0s[i];
				v1s[i] = x.v1s[i];
				v2s[i] = x.v2s[i];
				v3s[i] = x.v3s[i];
				v0e[i] = x.v0e[i];
				v1e[i] = x.v1e[i];
				v2e[i] = x.v2e[i];
				v3e[i] = x.v3e[i];
				err[i] = x.err[i];
				tol[i] = x.tol[i];
			}

			return *this;
		}
	};

	inline Scalar calculate_vf(const CCDdata &data_in, const BoxPrimitives &bp)
	{
		const Scalar v = (data_in.v0e[bp.dim] - data_in.v0s[bp.dim]) * bp.t + data_in.v0s[bp.dim];
		const Scalar t0 = (data_in.v1e[bp.dim] - data_in.v1s[bp.dim]) * bp.t + data_in.v1s[bp.dim];
		const Scalar t1 = (data_in.v2e[bp.dim] - data_in.v2s[bp.dim]) * bp.t + data_in.v2s[bp.dim];
		const Scalar t2 = (data_in.v3e[bp.dim] - data_in.v3s[bp.dim]) * bp.t + data_in.v3s[bp.dim];
		const Scalar pt = (t1 - t0) * bp.u + (t2 - t0) * bp.v + t0;
		return (v - pt);
	}

	inline Scalar calculate_ee(const CCDdata &data_in, const BoxPrimitives &bp)
	{
		const Scalar edge0_vertex0 = (data_in.v0e[bp.dim] - data_in.v0s[bp.dim]) * bp.t + data_in.v0s[bp.dim];
		const Scalar edge0_vertex1 = (data_in.v1e[bp.dim] - data_in.v1s[bp.dim]) * bp.t + data_in.v1s[bp.dim];
		const Scalar edge1_vertex0 = (data_in.v2e[bp.dim] - data_in.v2s[bp.dim]) * bp.t + data_in.v2s[bp.dim];
		const Scalar edge1_vertex1 = (data_in.v3e[bp.dim] - data_in.v3s[bp.dim]) * bp.t + data_in.v3s[bp.dim];
		const Scalar result = ((edge0_vertex1 - edge0_vertex0) * bp.u + edge0_vertex0)
							  - ((edge1_vertex1 - edge1_vertex0) * bp.v + edge1_vertex0);

		return result;
	}

	inline bool sum_no_larger_1(const Scalar num1, const Scalar num2)
	{
#ifdef GPUTI_USE_DOUBLE_PRECISION
		if (num1 + num2 > 1 / (1 - DBL_EPSILON))
		{
			return false;
		}
#else
		if (num1 + num2 > 1 / (1 - FLT_EPSILON))
		{
			return false;
		}
#endif
		return true;
	}

	//TODO maybe remove (or remove comment)
	class interval_pair
	{
	public:
		// this function do the bisection
		inline interval_pair(const Singleinterval &itv)
		{
			const Scalar c = (itv.first + itv.second) / 2;
			first.first = itv.first;
			first.second = c;
			second.first = c;
			second.second = itv.second;
		};
		interval_pair(){};

		Singleinterval first;
		Singleinterval second;
	};

	double ccd(const std::vector<std::array<std::array<Scalar, 3>, 8>> &V, bool is_edge, double max_t, double tolerace, double ms);
} // namespace cpu_ti
