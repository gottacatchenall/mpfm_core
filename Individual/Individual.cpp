#include "Individual.h"
#include "MPFM_Instance.h"

int Individual::id_counter = 0;

Individual::Individual(MPFM_Instance* run, int pop_born_in, int parent1_home, int parent2_home, bool parent_was_migrant){
    this->id = this->id_counter++;
    this->mpfm = run;
    this->pop_born_in = pop_born_in;
    this->current_pop = pop_born_in;
    this->parent1_home = parent1_home;
    this->parent2_home = parent2_home;
    this->parent_was_migrant = (bool) parent_was_migrant;

    this->has_migrated = false;

    int genome_size = this->mpfm->N_LOCI;

    this->haplotype0 = new double[genome_size];
    this->haplotype1 = new double[genome_size];

    for (int l = 0; l < genome_size; l++){
        this->haplotype0[l] = 0.0;
        this->haplotype1[l] = 0.0;
    }
}

Individual::~Individual(){
    delete this->haplotype0;
    delete this->haplotype1;
}

void Individual::calc_fitness(std::vector<double> efs){
    int n_loci = this->mpfm->selection_strengths.size();

    int n_ef = this->mpfm->params["EF_NUMBER"];
    assert(n_ef == efs.size());

    std::vector<double> sel_strs = this->mpfm->selection_strengths;

    double sel_str, L_l0, L_l1, diff, gaussian, s_i, w_i, theta_i;
    int ef;
    double w = 1.0;
    double sigma_s = 1.0;
    double s_max_i = -1 * this->mpfm->params["SELECTION_STRENGTH"];

    for (int l = 0; l < n_loci; l++){
        sel_str = sel_strs[l];
        ef = this->mpfm->ef_genome_map[l];
        if (ef >= 0){
            theta_i = efs[ef];
            assert(theta_i >= 0 && theta_i <= 1);
            L_l0 = this->haplotype0[l];
            diff = L_l0 - theta_i;
            gaussian = exp( (-1*(pow(diff,2))/sigma_s ) );
            s_i = s_max_i * sel_str * gaussian;
            w_i = 1.0 + s_i;
            w = w * w_i;

            L_l1 = this->haplotype1[l];
            //printf("theta: %f, l: %d h1:%f h2:%f\n", theta_i, l, L_l0, L_l1);
            diff = L_l1 - theta_i;
            gaussian = exp( (-1*(pow(diff,2))/sigma_s ) );
            s_i = s_max_i * sel_str * gaussian;
            w_i = 1.0 + s_i;
            w = w * w_i;
        }
    }
    this->w = w;
}

void Individual::gen_haplotype(Individual* parent, int haplo){
    std::vector<int> chromosome_map = this->mpfm->chromosome_map;
    std::vector<double> map_dists = this->mpfm->map_dists;

    int n_loci = this->mpfm->N_LOCI;

    int current_chr = 0;
    int parent_haplotype = this->mpfm->uniform_int(0,1);

    double mutation_rate = this->mpfm->params["MUTATION_RATE"];

    int n_mutations = this->mpfm->binomial(mutation_rate, n_loci);

    std::vector<int> mutation_sites;
    mutation_sites.reserve(n_mutations);
    for (int i = 0; i < n_mutations; i++){
        int random_site = this->mpfm->uniform_int(0, n_loci-1);
        mutation_sites.push_back(random_site);
    }

    for (int l = 0; l < n_loci; l++){
        bool mutation_this_site = false;
        if (n_mutations > 0){
            for (int site : mutation_sites){
                if (site == l){
                    mutation_this_site = true;
                }
            }
        }
        if (mutation_this_site){
            if (haplo == 0){
                this->haplotype0[l] = this->mpfm->uniform_int(0,1);
            }
            else if (haplo == 1){
                this->haplotype1[l] = this->mpfm->uniform_int(0,1);
            }
        }
        else {
            if (l > 0){
                // Is this new chromo
                if (chromosome_map[l] != chromosome_map[l-1]){
                    parent_haplotype = this->mpfm->uniform_int(0,1);
                }

                // Is there crossing over
                double d = map_dists[l] - map_dists[l-1];
                double p = 0.5*(1 - exp(double(-2*d)/double(100)));

                // draws every locus are very slow
                if (this->mpfm->uniform_01() < p){
                    parent_haplotype = !parent_haplotype;
                }
            }

            if (haplo == 0){
                if (parent_haplotype == 0){
                    this->haplotype0[l] = parent->haplotype0[l];
                }
                else if (parent_haplotype == 1){
                    this->haplotype0[l] = parent->haplotype1[l];
                }
            }

            else if (haplo == 1){
                if (parent_haplotype == 0){
                    this->haplotype1[l] = parent->haplotype0[l];
                }
                else if (parent_haplotype == 1){
                    this->haplotype1[l] = parent->haplotype1[l];
                }
            }
        }
    }

    //check to make sure genome is init

    for (int l = 0; l < n_loci; l++){
        double val;
        if (haplo == 0){
            val = this->haplotype0[l];
        }
        else if (haplo == 1){
            val = this->haplotype1[l];
        }

        if (val < 0 || val > 1){
            assert(false && "error in genome!\n");
        }
    }



}

/*
std::vector<double> Individual::get_crossing_over_points(){
    double cM = this->mpfm->GENOME_LENGTH;

    double exp_num_co = cM / 100.0;
    int n_crossover_events = this->mpfm->poisson(this->mpfm->main_gen, exp_num_co);

    std::vector<double> crossing_over_points;

    int n_loci = this->mpfm->N_LOCI;

    int r;
    bool exists;
    for(int i = 0; i < n_crossover_events; i++){
        r = this->mpfm->uniform_real(this->mpfm->main_gen, 0,n_loci-1);
    }

    std::sort(crossing_over_points.begin(), crossing_over_points.end());
    return crossing_over_points;
}
*/
