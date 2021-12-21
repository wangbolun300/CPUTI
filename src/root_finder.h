#pragma once
#include<cputi/queue.h>
#include <math.h>       /* fabs */
#include <vector>
#include <array>

void vertexFaceCCD(const CCDdata &data_in,const CCDConfig& config, CCDOut& out);
void edgeEdgeCCD(const CCDdata &data_in,const CCDConfig& config, CCDOut& out);
double return_time();
double return_time_vf();
void memory_pool_ccd_run(
    const std::vector<std::array<std::array<Scalar, 3>, 8>> &V, bool is_edge,
    std::vector<int> &result_list, double &run_time,
    std::vector<Scalar> &time_impact, int parallel_nbr);