#ifndef POPULATION_H
#define POPULATION_H

#include "include.h"

class Population{
    private:
    public:
        // Constructor
        Population(MPFM_Instance* run, int id, int index, double x, double y, double K, std::vector<double> efs);

        // Data
        int id;
        int index;
        double x;
        double y;
        double k;
        double mean_w;
        std::vector<double> efs;
        MPFM_Instance* mpfm;

        std::vector<Individual*> indivs;

};

#endif
