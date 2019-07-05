#include "DataWrangler.h"
#include "MPFM_Instance.h"
#include "Individual.h"
#include "Population.h"

DataWrangler::DataWrangler(MPFM_Instance* run){
    this->mpfm = run;
    this->sample_size = SAMPLE_ALL;
    // ===============================================
    // Genome Logging Inititalization
    // ===============================================
    std::string genome_header = "generation,locus,selection_weight,chromosome,global_polymorphism_ct,mean_polymorphism_ct_per_pop,global_f_st,mean_global_ld,sample_size\n";
    this->genome_file = this->init_file("genome.csv", genome_header);

    // ===============================================
    // Individual Populations Logging Inititalization
    // ===============================================
    std::string individual_populations_header = "generation,pop1,x,y,w_mean,prop_of_k,effective_migration,prop_of_loci_fixed,mean_polymorphism_ct_per_locus,";
    int n_ef = this->mpfm->params["EF_NUMBER"];
    for (int i = 0; i < n_ef; i++){
        individual_populations_header +=  "ef" + std::to_string(i) + ",";
    }
    individual_populations_header += "sample_size\n";
    this->individual_populations_file = this->init_file("individual_populations.csv", individual_populations_header);

    // ===============================================
    // Pairwise Populations Logging Inititalization
    // ===============================================
    std::string pairwise_populations_header = "generation,pop1,pop2,euclidean_distance,mean_ld,mean_f_st,eff_migration_pop1_to_pop2,eff_migration_pop2_to_pop1,ld_clustering,ld_network_threshold,sample_size\n";
    this->pairwise_populations_file = this->init_file("pairwise_populations.csv", pairwise_populations_header);

    int n_loci = this->mpfm->N_LOCI;
    for (int l = 0; l < n_loci; l++){
        this->genotype_database.push_back(std::vector<allele*>());
    }
}


std::ofstream* DataWrangler::init_file(std::string file_path, std::string header){
    std::ofstream* file = new std::ofstream;
    file->open(file_path.c_str(), std::ios::app);
    (*file) << header;
    file->flush();
    return file;
}

void DataWrangler::census(int gen, int sample_size){
    // Build allele table by sampling a certain num from each pop, special mode for full.
    this->construct_allele_table(sample_size);
    this->sample_size = sample_size;

    this->log_individual_population_data(gen);
    this->log_pairwise_population_data(gen);
    this->log_genome_data(gen);

    // dump buffer
    this->individual_populations_file->flush();
    this->genome_file->flush();
    this->pairwise_populations_file->flush();

    // clear data
    this->destruct_allele_table();

}

void DataWrangler::log_individual_population_data(int gen){
    int n_pops = this->mpfm->populations.size();

    for (int p = 0; p < n_pops; p++){
        Population* pop = this->mpfm->populations[p];
        double x = pop->x;
        double y = pop->y;
        double w_mean = pop->mean_w;
        double prop_of_k = double(this->mpfm->indivs_by_pop[p].size())/double(pop->k);
        double eff_migration = this->get_effective_migration_single_pop(p);
        double prop_of_loci_fixed = this->get_prop_of_loci_fixed(p);
        double mean_polymorphism_ct = this->get_mean_polymorphism_ct(p);
        std::vector<double> efs = pop->efs;
    //    printf("pop, w_mean: %d %f\n", p, w_mean);
        this->write_individual_pop_data(gen, p, x, y, w_mean, prop_of_k, eff_migration, prop_of_loci_fixed, mean_polymorphism_ct, efs);

    }

}

unsigned long int DataWrangler::nChoosek( unsigned long int n, unsigned long k ){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


void DataWrangler::log_pairwise_population_data(int gen){

    int n_pops = this->mpfm->populations.size();
    int n_loci = this->mpfm->N_LOCI;
    unsigned long int n_choose_3 = this->nChoosek(n_loci, 3);

    for (int p1 = 0; p1 < n_pops; p1++){
        for (int p2 = 0; p2 < n_pops; p2++){
            if  (p1 != p2){
                    //std::string pairwise_populations_header = "generation,pop1,pop2,euclidean_distance,mean_ld,mean_f_st,eff_migration_pop1_to_pop2,eff_migration_pop2_to_pop1,ld_clustering,ld_modularity,ld_network_threshold\n";

                    double** ld_matrix = this->get_ld_matrix(p1, p2, sample_size);

                    double euclidean_distance = get_euclidean_distance(p1, p2);
                    double mean_ld = get_mean_pooled_ld(ld_matrix);
                    double mean_fst = get_mean_pooled_fst(p1, p2);
                    double eff_migration_pop1_to_pop2 = get_eff_migration(p1, p2);
                    double eff_migration_pop2_to_pop1 = get_eff_migration(p1,p2);
                    //for ldnet_thresold in thresholds
                    std::vector<double> thres = {0.01, 0.03, 0.05, 0.1};
                    for (double t : thres){
                        double ld_clustering = get_ld_clustering(ld_matrix, t, n_choose_3);
                        this->write_pairwise_pop_data(gen, p1, p2, euclidean_distance, mean_ld, mean_fst, eff_migration_pop1_to_pop2, eff_migration_pop2_to_pop1, ld_clustering, t);
                    }

                    for (int l1 = 0; l1 < n_loci; l1++){
                        delete ld_matrix[l1];
                    }
                    delete ld_matrix;
            }
        }
    }

}

void DataWrangler::log_genome_data(int gen){
    int n_loci = this->mpfm->N_LOCI;
    for (int l = 0; l < n_loci; l++){
        double weight = this->mpfm->selection_strengths[l];
        int chromo = this->mpfm->chromosome_map[l];
        int poly_ct = this->get_poly_ct(l);
        double mean_poly_per_pop = this->get_mean_poly_per_pop(l);
        double f_st = 0;
        double mean_global_ld = this->get_mean_global_ld(l);

        this->write_genome_data(gen, l, weight, chromo, poly_ct, mean_poly_per_pop, f_st, mean_global_ld);
    }
}


// ================================================================
// Functions to compute stats
// for pairwise populations
// ================================================================
double** DataWrangler::get_ld_matrix(int p1, int p2, int sample_size){
    int n_loci = this->mpfm->N_LOCI;

    // init matrix
    double ** ld_matrix = new double*[n_loci];
    for (int l = 0; l < n_loci; l++){
        ld_matrix[l] = new double[n_loci];
    }
    for (int l1 = 0; l1 < n_loci; l1++){
        for (int l2 = l1+1; l2 < n_loci; l2++){
            ld_matrix[l1][l2] = 0.0;
        }
    }

    for (int l1 = 0; l1 < n_loci; l1++){
        for (int l2 = l1+1; l2 < n_loci; l2++){
            double ld_l1_l2 = this->calc_pairwise_ld(p1, p2, l1, l2);
            ld_matrix[l1][l2] = ld_l1_l2;
        }
    }
    return ld_matrix;
}

double DataWrangler::calc_pairwise_ld(int p1, int p2, int l1, int l2){
    // for each allele combination that exists, calc LD
    //printf("top calc pw (p1, p2, l1, l2): (%d, %d, %d, %d)\n", p1, p2 ,l1,l2);
    std::vector<allele*> l1_alleles = this->genotype_database[l1];
    std::vector<allele*> l2_alleles = this->genotype_database[l2];

    int eff_pop_size = 0;
    if (this->sample_size == SAMPLE_ALL){
        int n1 = this->mpfm->indivs_by_pop[p1].size();
        int n2 = this->mpfm->indivs_by_pop[p2].size();
//        printf("n1, n2: %d, %d\n", n1, n2);
        eff_pop_size = 2*(n1+n2);
    }
    else {
        eff_pop_size = 2*this->sample_size;
    }

    double ld_avg = this->calc_ld_between_two_loci(l1_alleles, l2_alleles, l1, l2, p1, p2, eff_pop_size);
    return ld_avg;
/*
        for (allele* l1_i : l1_alleles){
            int l1_ct = l1_i->ct_map[p1] + l1_i->ct_map[p2];
            std::vector<dependent_allele*> haplotype_data = l1_i->loci[l2];

            int l2_all_ct = 0;
            for (allele* l2_i : l2_alleles){
                bool seen_with_l1_i = false;
                int l2_haplotype_ct = 0;
                for (dependent_allele* l2_haplotype : haplotype_data){
                    double l2_haplo_val = l2_haplotype->allele_val;
                    if (l2_i->allele_val == l2_haplo_val){
                        seen_with_l1_i = true;
                        l2_haplotype_ct = l2_haplotype->ct_map[p1] + l2_haplotype->ct_map[p2];
                        break;
                    }
                }
                if (seen_with_l1_i){
                    l2_all_ct = l2_i->ct_map[p1] + l2_i->ct_map[p2];
                    if (l2_all_ct < l2_haplotype_ct){
                        assert(false);
                    }
                    double ld = calc_ld_from_ct(l1_ct, l2_all_ct, l2_haplotype_ct, eff_pop_size);
                    ld_avg += ld;
                    ct++;
                }
            }
        }

    return ld_avg/double(ct);*/
}

double DataWrangler::get_euclidean_distance(int p1, int p2){
    Population* pop1 = this->mpfm->populations[p1];
    Population* pop2 = this->mpfm->populations[p2];

    double x1 = pop1->x;
    double x2 = pop2->x;
    double y1 = pop1->y;
    double y2 = pop2->y;

    double deltax = x2-x1;
    double deltay = y2-y1;

    double dist =  sqrt(pow(deltax,2) + pow(deltay,2));
    return dist;
}
double DataWrangler::get_mean_pooled_ld(double** ld_matrix){
    int n_loci = this->mpfm->N_LOCI;
    double sum = 0;
    int ct = 0;
    for (int l1 = 0; l1 < n_loci; l1++){
        for (int l2 = l1+1; l2 < n_loci; l2++){
            sum += ld_matrix[l1][l2];
            ct++;
        }
    }

    return double(sum)/double(ct);
}
double DataWrangler::get_mean_pooled_fst(int p1, int p2){
    int n_loci = this->mpfm->N_LOCI;

    int n1 = this->mpfm->indivs_by_pop[p1].size();
    int n2 = this->mpfm->indivs_by_pop[p2].size();

    std::vector<allele*> alleles_this_locus;
    double f_st_sum = 0;

    for (int l = 0; l < n_loci; l++){
        alleles_this_locus = this->genotype_database[l];

        double J_1 = 0.0;
        double J_2 = 0.0;
        double J_12 = 0.0;
        double J_T = 0.0;

        // Let f_xy be the frequency of the yth allele in xth pop
        for (allele* al_i : alleles_this_locus){

            double f_1i = al_i->ct_map[p1] / double(n1);
            double f_2i = al_i->ct_map[p2] / double(n2);

            double w_1 = double(n1)/double(n1+n2);
            double w_2 = 1.0 - w_1;

            J_1 += (f_1i * f_1i);
            J_2 += (f_2i * f_2i);

            J_12 += (f_1i * f_2i);

            double sum_internal =  w_1*f_1i + w_2*f_2i;
            J_T += sum_internal * sum_internal;
        }

        double D_ST = 0.25 * (0.5*(J_1 + J_2) - J_12);
        double F_ST = D_ST / J_T;
        f_st_sum += F_ST;
    }

    double mean_fst = f_st_sum / double(n_loci);
    return mean_fst;
}
double DataWrangler::get_eff_migration(int p_from, int p_to){
    std::vector<Individual*> indivs = this->mpfm->indivs_by_pop[p_to];
    int eff_migrant_ct = 0;
    for (Individual* indiv : indivs){
        if (indiv->parent1_home == p_from || indiv->parent2_home == p_from){
            eff_migrant_ct++;
        }
    }
    return double(eff_migrant_ct)/double(indivs.size());
}

double DataWrangler::get_ld_clustering(double** ld_matrix, double ld_threshold, int nchoose3){
    int n_loci = this->mpfm->N_LOCI;
    int A[n_loci][n_loci];
    int Asquared[n_loci][n_loci];
    int Acubed[n_loci][n_loci];

    for (int i = 0; i < n_loci; i++){
        for (int j = 0; j < n_loci; j++){
            if (fabs(ld_matrix[i][j]) > ld_threshold){
                A[i][j] = 1;
            }
            else{
                A[i][j] = 0;
            }
        }
    }

    // set tmp = A^2
    for (int i = 0; i < n_loci; i++){
        for (int j = 0; j < n_loci; j++){
            Asquared[i][j] = 0;
            for (int k = 0; k < n_loci; k++){
                Asquared[i][j] += A[i][k] * A[k][j];
            }
        }
    }

    for (int i = 0; i < n_loci; i++){
        for (int j = 0; j < n_loci; j++){
            Acubed[i][j] = 0;
            for (int k = 0; k < n_loci; k++){
                Acubed[i][j] += Asquared[i][k] * A[k][j];
            }
        }
    }

    int trA3 = 0;
    int denom = 0;

    for (int i = 0; i < n_loci; i++){
        for (int j = 0; j < n_loci; j++){
            if (i != j){
                denom += Asquared[i][j];
            }
        }
        trA3 += Acubed[i][i];
    }

//      printf("trA3: %d, denom: %d, num trips: %d\n", trA3, denom, 3*nchoose3);

/*    double nchoose3 = nfac / nmkfac;
    nchoose3 = nchoose3 / kfac;
    printf("nfac: %d", nchoose3);
//    printf("denom: %d", denom2);
//    printf("nfac: %d", nfac);
//    printf("finna assert!\n");
    printf("denom: %d\n", denom);
    assert(denom == nchoose3);*/


    double clus;
    if (denom > 0){
        clus = double(trA3)/double(denom);
    }
    else {
        clus = 0;
    }
    return clus;
}

// ================================================================
// Functions to calculate summary stats
// for individual populations
// ================================================================
double DataWrangler::get_effective_migration_single_pop(int p){
    std::vector<Individual*> indivs = this->mpfm->indivs_by_pop[p];
    int ct = 0;
    for (Individual* indiv : indivs){
        if (indiv->parent_was_migrant){
            ct++;
        }
    }
    double eff_mig = double(ct)/double(indivs.size());
    return eff_mig;
}

double DataWrangler::get_prop_of_loci_fixed(int p){
    int n_loci = this->mpfm->N_LOCI;
    int fixed_ct = 0;
    for (int l = 0; l < n_loci; l++){
        int als_seen = 0;
        for (allele* al : this->genotype_database[l]){
            if (al->ct_map[p] > 0){
                als_seen++;
            }
        }
        if (als_seen == 1){
            fixed_ct++;
        }
    }
    return double(fixed_ct)/double(n_loci);
}

double DataWrangler::get_mean_polymorphism_ct(int p){
    int n_loci = this->mpfm->N_LOCI;

    int poly_ct_total = 0;
    for (int l = 0; l < n_loci; l++){
        int als_seen = 0;
        for (allele* al : this->genotype_database[l]){
            if (al->ct_map[p] > 0){
                als_seen++;
            }
        }
        poly_ct_total += als_seen;
    }
    return double(poly_ct_total)/double(n_loci);
}
// ================================================================
// Functions to calculate summary stats
// by loci
// ================================================================
int DataWrangler::get_poly_ct(int l){
    int poly_ct = 0;
    for (allele* al : this->genotype_database[l]){
        if (al->n_total > 0){
            poly_ct++;
        }
    }
    return poly_ct;
}

double DataWrangler::get_mean_poly_per_pop(int l){
    int n_pops = this->mpfm->populations.size();

    double poly_sum = 0;
    for (allele* al : this->genotype_database[l]){
        for (int p = 0; p < n_pops; p++){
            if (al->ct_map[p] > 0){
                poly_sum++;
            }
        }
    }

    return double(poly_sum)/double(n_pops);
}
double DataWrangler::get_f_st_by_locus(int l){
    return 0;
}
double DataWrangler::get_mean_global_ld(int l){
    int n_loci = this->mpfm->N_LOCI;
    int n_pops = this->mpfm->populations.size();
    int eff_pop_size = 0;
    if (this->sample_size == SAMPLE_ALL){
        int global_pop_size = 0;
        for (Individual * indiv : *(this->mpfm->indivs)){
            if (indiv->current_pop >= 0){
                global_pop_size++;
            }
        }
        eff_pop_size = 2*global_pop_size;
    }
    else {
        eff_pop_size = 2*n_pops*this->sample_size;
    }

    double ld_sum = 0;
    double this_locus_combo_ld_avg;
    int ct = 0;

    for (int l2 = 0; l2 < n_loci; l2++){
        if (l < l2){
            std::vector<allele*> l1_alleles = this->genotype_database[l];
            std::vector<allele*> l2_alleles = this->genotype_database[l2];

            this_locus_combo_ld_avg = this->calc_ld_between_two_loci(l1_alleles, l2_alleles, l, l2, -1, -1, eff_pop_size);
            ld_sum += this_locus_combo_ld_avg;
            ct++;
        }
        else if (l > l2){
            std::vector<allele*> l1_alleles = this->genotype_database[l2];
            std::vector<allele*> l2_alleles = this->genotype_database[l];

            this_locus_combo_ld_avg = this->calc_ld_between_two_loci(l1_alleles, l2_alleles, l2, l, -1, -1, eff_pop_size);
            ld_sum += this_locus_combo_ld_avg;
            ct++;
        }
    }
    double avg_ld = ld_sum / double(ct);
    return avg_ld;
}
// ================================================================
// linkage
// disequilibrium
// ================================================================
double DataWrangler::calc_ld_between_two_loci(std::vector<allele*> l1_alleles, std::vector<allele*> l2_alleles, int l1, int l2, int p1, int p2, int eff_pop_size){
    double ld_sum = 0;
    int ct = 0;
    int l1_ct, l2_all_ct, l2_haplotype_ct;

    bool global = false;
    if (p1 == -1 || p2 == -1){
        global = true;
    }


    for (allele* l1_i : l1_alleles){
        if (!global){
            l1_ct = l1_i->ct_map[p1] + l1_i->ct_map[p2];
        }
        else {
            l1_ct = l1_i->n_total;
        }
        std::vector<dependent_allele*> haplotype_data = l1_i->loci[l2];

        l2_all_ct = 0;
        for (allele* l2_i : l2_alleles){
            bool seen_with_l1_i = false;
            l2_haplotype_ct = 0;
            for (dependent_allele* l2_haplotype : haplotype_data){
                double l2_haplo_val = l2_haplotype->allele_val;
                if (l2_i->allele_val == l2_haplo_val){
                    seen_with_l1_i = true;
                    if (!global){
                        l2_haplotype_ct = l2_haplotype->ct_map[p1] + l2_haplotype->ct_map[p2];
                    }
                    else {
                        l2_haplotype_ct = l2_haplotype->n_total;
                    }
                    break;
                }
            }
            if (seen_with_l1_i){
                if (!global){
                    l2_all_ct = l2_i->ct_map[p1] + l2_i->ct_map[p2];
                }
                else {
                    l2_all_ct = l2_i->n_total;
                }

                if (l2_all_ct < l2_haplotype_ct){
                    assert(false);
                }


                double ld = this->calc_ld_from_ct(l1_ct, l2_all_ct, l2_haplotype_ct, eff_pop_size);
                ld_sum += ld;
                ct++;
            }
        }
    }
    return double(ld_sum)/double(ct);
}

double DataWrangler::calc_ld_from_ct(int ct_a, int ct_b, int ct_ab, int eff_pop_size){
    //printf("eff_pop_size: %f\n", eff_pop_size);
    //printf("ct_a, ct_b, ct_ab: %d, %d, %d\n", ct_a, ct_b, ct_ab);

    double p_a = double(ct_a)/double(eff_pop_size);
    double p_b = double(ct_b)/double(eff_pop_size);
    double p_ab = double(ct_ab)/double(eff_pop_size);
    double D;
    D = p_ab - p_a*p_b;


    if (D == 0){
        return 0;
    }

    double denom = p_a*(1.0-p_a)*p_b*(1.0-p_b);
    double d2 = pow(D, 2);

    double LD = d2/denom;

    return LD;
}

// ================================================================
// Functions to manage allele frequencies
// and haplotype data
// ================================================================
void DataWrangler::construct_allele_table(int sample_size){
    int n_pops = this->mpfm->populations.size();
    int n_loci = this->mpfm->N_LOCI;
    for (int l = 0; l < n_loci; l++){
        for (allele* al : this->genotype_database[l]){
            for (int p = 0; p < n_pops; p++){
                al->ct_map[p] = 0;
            }
            al->n_total = 0;

            for (int l2 = l+1; l2 < n_loci; l2++){
                for (dependent_allele* dep_al: al->loci[l2]){
                    dep_al->n_total = 0;
                    for (int p = 0; p < n_pops; p++){
                        dep_al->ct_map[p] = 0;
                    }
                }
            }
        }
    }

    for (int p = 0; p < n_pops; p++){
        // whole pop
        std::vector<Individual*> indivs;

        // choose whether is easier to include sample, or subtrcat, pop - sample
        if (sample_size == -1){
            indivs = this->mpfm->indivs_by_pop[p];
        }
        else {
            indivs = this->mpfm->sample_n_random_indivs(p, sample_size);
        }

//        printf("genotyping indivs from p = %d\n", p);

        for (Individual* indiv : indivs){
            this->genotype_indiv(indiv);
        }
    }
}

void DataWrangler::genotype_indiv(Individual* indiv){
    int n_loci = this->mpfm->N_LOCI;

    int this_pop = indiv->current_pop;
//    printf("genotyping indiv : %x\n", indiv);


    // check this indivs genome
/*    printf("haplo 0: [");
    for (int l = 0; l < n_loci; l++){
        printf("%f ",indiv->haplotype0[l]);
    }
    printf("]\n");
    printf("haplo 1: [");
    for (int l = 0; l < n_loci; l++){
        printf("%f ",indiv->haplotype1[l]);
    }
    printf("]\n");*/



    for (int l1 = 0; l1 < n_loci; l1++){
        double al0_val = indiv->haplotype0[l1];
        double al1_val = indiv->haplotype1[l1];

        //printf("updating allele* at locus %d\n", l1);
        allele* al0 = this->update_genotype_database(l1, al0_val, this_pop);
        allele* al1 = this->update_genotype_database(l1, al1_val, this_pop);;

    /*    if (al0_val == al1_val){
            al1 = al0;
        }
        else{
            allele* al1 = this->update_genotype_database(l1, al1_val, this_pop);
        }*/

        for (int l2 = l1+1; l2 < n_loci; l2++){
            double l2_haplo0_val = indiv->haplotype0[l2];
            double l2_haplo1_val = indiv->haplotype1[l2];

            this->add_haplotype_data(al0, l2, l2_haplo0_val, this_pop);
            this->add_haplotype_data(al1, l2, l2_haplo1_val, this_pop);
        }
    }
}

allele* DataWrangler::update_genotype_database(int locus, double allele_val, int pop){
    bool new_al = true;

    //printf("looking for allele: %f at locus %d in pop %d", allele_val, locus, pop);

    for (allele* al : this->genotype_database[locus]){
        if (al->allele_val == allele_val){
            al->n_total++;
            al->ct_map[pop]++;
            //printf("returning %x\n", al);
            return al;
        }
    }
    allele* new_allele;
    if (new_al){
        int n_loci = this->mpfm->N_LOCI;
        new_allele = new allele(locus, allele_val, n_loci);
        for (int l = 0; l < n_loci; l++){
            new_allele->ct_map.push_back(0);
        }
        new_allele->n_total++;
        new_allele->ct_map[pop]++;

        this->genotype_database[locus].push_back(new_allele);
        //printf("returning %x\n", new_allele);
        return new_allele;
    }
    assert(false && "reached end of update_genotype_database() w/o returning something");
}

void DataWrangler::add_haplotype_data(allele* primary_allele, int dependent_locus, double dependent_allele_val, int pop){
    for (dependent_allele* al2 : primary_allele->loci[dependent_locus]){
        if (al2->allele_val == dependent_allele_val){
            al2->n_total++;
            al2->ct_map[pop]++;
            return;
        }
    }

    dependent_allele* new_dep_allele = new dependent_allele(dependent_locus, dependent_allele_val);
    //int n_pops = this->mpfm->populations.size();
    int n_loci = this->mpfm->N_LOCI;
    for (int l = 0; l < n_loci; l++){
        new_dep_allele->ct_map.push_back(0);
    }
    new_dep_allele->n_total++;
    new_dep_allele->ct_map[pop]++;
    primary_allele->loci[dependent_locus].push_back(new_dep_allele);
}

void DataWrangler::destruct_allele_table(){
    int n_loci = this->mpfm->N_LOCI;
    int n_pops = this->mpfm->populations.size();
    for (int l = 0; l < n_loci; l++){
        for (allele* al : this->genotype_database[l]){
            for (int p = 0; p < n_pops; p++){
                al->ct_map[p] = 0;
            }
            al->n_total = 0;

            for (int l2 = l+1; l2 < n_loci; l2++){
                for (dependent_allele* dep_al: al->loci[l2]){
                    dep_al->n_total = 0;
                    for (int p = 0; p < n_pops; p++){
                        dep_al->ct_map[p] = 0;
                    }
                }
            }
        }
    }
}

// ================================================================
// Write to file
// methods
// ================================================================
void DataWrangler::write_individual_pop_data(int gen, int pop, double x, double y, double w_mean, double prop_of_k, double eff_mig, double prop_loci_fixed, double mean_polymorphism_ct_per_locus, std::vector<double> efs){
    (*this->individual_populations_file) << gen << "," << pop << "," << x << "," << y << "," << w_mean << "," << prop_of_k << "," << eff_mig << "," << prop_loci_fixed << "," << mean_polymorphism_ct_per_locus << ",";
    double sp = this->mpfm->sample_prop;

    for (double val : efs){
        (*this->individual_populations_file) << val << ",";
    }
    (*this->individual_populations_file) << "," << sp << "\n";
}


void DataWrangler::write_pairwise_pop_data(int gen, int pop1, int pop2, double dist, double mean_ld, double mean_fst, double mig_p1_to_p2, double mig_p2_to_p1, double ld_clustering, double ld_network_threshold){
    double sp = this->mpfm->sample_prop;
    (*this->pairwise_populations_file) << gen << "," << pop1 << "," << pop2 << "," << dist << "," << mean_ld << "," << mean_fst << "," << mig_p1_to_p2 << "," << mig_p2_to_p1 << "," << ld_clustering << "," << ld_network_threshold << "," << sp << "\n";
}

void DataWrangler::write_genome_data(int gen, int locus, double weight, int chromo, int poly_ct, double mean_poly_per_pop, double f_st, double mean_global_ld){
    double sp = this->mpfm->sample_prop;
        (*this->genome_file) << gen << "," << locus  << "," << weight  << "," << chromo  << "," << poly_ct  << "," << mean_poly_per_pop  << "," << f_st  << "," << mean_global_ld << "," << sp << "\n";
}
