
#ifndef DATA_WRANGLER_H
#define DATA_WRANGLER_H

#include "include.h"

#define SAMPLE_ALL -1

typedef struct dependent_allele{
    dependent_allele(int l, double val){
        locus = l;
        allele_val = val;
        n_total = 0;
    };
    int locus;
    double allele_val;
    int n_total;
    std::vector<int> ct_map;
} dependent_allele;

typedef struct allele{
    allele(int l, double val, int n_loci){
        locus = l;
        allele_val = val;
        n_total = 0;

        for (int i = 0; i < n_loci; i++){
            loci.push_back(std::vector<dependent_allele*>());
        }
    };

    int locus;
    double allele_val;
    int n_total;
    std::vector<int> ct_map;
    std::vector<std::vector<dependent_allele*>> loci;
} allele;

class DataWrangler{
    private:
    public:
        std::ofstream* individual_populations_file;
        std::ofstream* pairwise_populations_file;
        std::ofstream* genome_file;
        MPFM_Instance* mpfm;

        std::vector<std::vector<allele*>> genotype_database;
        int sample_size;

        DataWrangler(MPFM_Instance* run);
        std::ofstream* init_file(std::string file_path, std::string header);
        void census(int gen, int sample_size);
        void genotype_indiv(Individual* indiv);

        void log_individual_population_data(int gen);
        void log_pairwise_population_data(int gen);
        void log_genome_data(int gen);

        // Functions to compute stats for pairwise populations
        double** get_ld_matrix(int p1, int p2, int sample_size);
        double calc_pairwise_ld(int p1, int p2, int l1, int l2);
        double calc_ld_from_ct(int ct_a, int ct_b, int ct_ab, int eff_pop_size);
        double get_euclidean_distance(int p1, int p2);
        double get_mean_pooled_ld(double** ld_matrix);
        double get_mean_pooled_fst(int p1, int p2);
        double get_eff_migration(int p_from, int p_to);
        double get_ld_clustering(double** ld_matrix, double ld_threshold);


        // Functions to calculate summary stats for individual populations
        double get_effective_migration_single_pop(int p);
        double get_prop_of_loci_fixed(int p);
        double get_mean_polymorphism_ct(int p);

        // Allele table
        void construct_allele_table(int sample_size);
        allele* update_genotype_database(int locus, double allele_val, int pop);
        void add_haplotype_data(allele* primary_allele, int dependent_locus, double dependent_allele_val, int pop);
        void destruct_allele_table();

        // Functions to write to files
        void write_individual_pop_data(int gen, int pop, double x, double y, double w_mean, double prop_of_k, double eff_mig, double prop_loci_fixed, double mean_polymorphism_ct_per_locus, std::vector<double> efs);

        void write_pairwise_pop_data(int gen, int pop1, int pop2, double dist, double mean_ld, double mean_fst, double mig_p1_to_p2, double mig_p2_to_p1, double ld_clustering, double ld_network_threshold);
};

#endif
