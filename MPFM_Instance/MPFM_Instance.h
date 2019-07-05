#ifndef MPFM_INSTANCE_H
#define MPFM_INSTANCE_H

#include "include.h"

class MPFM_Instance{
    private:
    public:
        MPFM_Instance();

        std::vector<Population*> populations;
        std::vector<Individual*>* indivs;
        std::vector<Individual*>* next_gen;
        std::vector<std::vector<double>> dis_kernel;
        std::vector<std::vector<double>> burn_in_dis_kernel;
        std::vector<std::vector<double>> fragmentation_dis_kernel;
        std::vector<double> selection_strengths;
        std::vector<int> chromosome_map;
        std::vector<double> map_dists;
        std::vector<int> ef_genome_map;
        std::vector<int> init_polymorphism_ct;
        std::vector<double> length_of_each_chromo;

        std::vector<std::vector<Individual*>> indivs_by_pop;
        DataWrangler* data_wrangler;

        // Methods for initialization
        void initialize();
        void read_params_ini();
        void read_pops_ini();
        void read_dis_kernel(std::string path);
        void read_burnin_dis_kernel();
        void read_fragmentation_dis_kernel();
        void read_genome_dict();
        void init_random_generators();
        void init_individuals();
        void init_genomes();

        // Run generations
        void burn_in();
        void fragmentation();
        void run_generation(int gen, int log_freq);
        void logging(int gen);

        // Core methods
        void dispersal();
        void selection();
        void reproduction();

        // Helpers
        void split_indivs_by_pop();
        std::vector<Individual*> pull_n_random_indivs(int pop_number, int num_indivs);
        int get_pop_index_from_id(int id);
        void calc_indivs_fitnesses(int pop);
        void run_beverton_holt_selection(int pop);
        double beverton_holt_prob(double n, double Kprime);
        Individual* get_random_individual(int pop);
        std::vector<Individual*> sample_n_random_indivs(int pop, int n);

        // Parameters
        std::unordered_map<std::string, float> params;
        int N_POPULATIONS;
        int N_INDIVIDUALS;
        int N_LOCI;
        double BASE_MIGRATION_RATE;
        int sample_size;
        double sample_prop;

        // Random generation
        std::mt19937* main_gen;
        double uniform_01();
        int uniform_int(int lo, int hi);
        int binomial(double p, int n);
        int poisson(double lambda);

};

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


#endif
