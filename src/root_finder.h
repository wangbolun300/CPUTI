#pragma once
#include<cputi/queue.h>
#include <math.h>       /* fabs */

void vertexFaceCCD(const CCDdata &data_in,const CCDConfig& config, CCDOut& out);
void edgeEdgeCCD(const CCDdata &data_in,const CCDConfig& config, CCDOut& out);
double return_time();
