

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "include.h"

class Individual{
    private:
    public:
        static int id_counter;
        int id;
        int pop_born_in;
        int current_pop;

        double* haplotype0;
        double* haplotype1;
        double w;
        double exp_num_off;
        MPFM_Instance* mpfm;

        bool has_migrated;
        bool parent_was_migrant;

        int parent1_home;
        int parent2_home;

        Individual(MPFM_Instance* run, int pop_born_in, int parent1_home, int parent2_home, bool parent_was_migrant);
        ~Individual();

        void calc_fitness(std::vector<double> efs);
        void gen_haplotype(Individual* parent, int haplo);

};

#endif
