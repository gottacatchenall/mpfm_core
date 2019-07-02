#include "Population.h"

Population::Population(MPFM_Instance* run, int id, int index, double x, double y, double k, std::vector<double> efs){
    this->id = id;
    this->index = index;
    this->x = x;
    this->y = y;
    this->k= k;
    this->efs = efs;
    this->mean_w = 0;
}
