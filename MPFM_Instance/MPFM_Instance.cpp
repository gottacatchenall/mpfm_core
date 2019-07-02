#include "MPFM_Instance.h"
#include "DataWrangler.h"
#include "Population.h"
#include "Individual.h"

// ================================================================
// Constructor
// ================================================================

MPFM_Instance::MPFM_Instance(){

}

// ================================================================
// Run generations
// ================================================================
void MPFM_Instance::burn_in(){
    int n_gen_burnin = this->params["BURN_IN_LENGTH"];
    for (int gen = 0; gen < n_gen_burnin; gen++){
        this->run_generation(gen);
        if (gen % 10 == 0){
            printf("%d\n", gen);
        }
    }
}

void MPFM_Instance::fragmentation(){
    int n_gen_burnin = this->params["BURN_IN_LENGTH"];
    int n_gen_fragmentation = this->params["FRAGMENTATION_LENGTH"];
    int stop_gen = n_gen_fragmentation + n_gen_burnin;

    this->dis_kernel = this->fragmentation_dis_kernel;

    for (int gen = n_gen_burnin+1; gen < stop_gen; gen++){
        this->run_generation(gen);
        if (gen % 10 == 0){
            printf("%d\n", gen);
        }
    }
}

void MPFM_Instance::run_generation(int gen){
    int log_freq = this->params["CENSUS_FREQ"];
    this->dispersal();
    this->selection();

    if (gen % log_freq == 0){
        this->logging(gen);
    }

    this->reproduction();
}

void MPFM_Instance::logging(int gen){
    this->split_indivs_by_pop();
    printf("about to census\n");
    this->data_wrangler->census(gen, SAMPLE_ALL);
    printf("census done\n");
}


// ================================================================
// core
// methods
// ================================================================

void MPFM_Instance::dispersal(){
    int n_pops = this->populations.size();

    double base_migration = this->params["BASE_MIGRATION_RATE"];

    this->split_indivs_by_pop();

    for (int p1 = 0; p1 < n_pops; p1++){
        std::vector<Individual*> indivs = indivs_by_pop[p1];
        int n = indivs.size();

        double exp_num_migrants = n * base_migration;

        double rem = exp_num_migrants - int(floor(exp_num_migrants));

        int n_migrants =  int(floor(exp_num_migrants));

        if (this->uniform_01() < rem){
            n_migrants++;
        }

        std::vector<Individual*> migrants = pull_n_random_indivs(p1, n_migrants);
        n_migrants = migrants.size();

        for (int i = 0; i < n_migrants; i++){
            Individual* indiv = migrants[i];

            double draw = this->uniform_01();
            double s = 0;
            for (int p2 = 0; p2 < n_pops; p2++){
                s += this->dis_kernel[p1][p2];
                if (s > draw){
                    indiv->current_pop = p2;
                    indiv->has_migrated = true;
                    break;
                }
            }
        }
    }
}

void MPFM_Instance::selection(){
    // indiv->calc_fitness();
    int n_pops = this->populations.size();
    for (int p = 0; p < n_pops; p++){
        this->calc_indivs_fitnesses(p);
        this->run_beverton_holt_selection(p);
    }
}

void MPFM_Instance::reproduction(){
    int base_fec = this->params["AVG_NUM_OFFSPRING_PER_INDIV"];
    int n_pops = this->populations.size();
    this->split_indivs_by_pop();

    std::vector<Individual*> indivs;

    Individual* random1;
    Individual* random2;
    Individual* offspring;
    Population* pop;
    bool parent_migrated;

    for (int p = 0; p < n_pops; p++){
        indivs = this->indivs_by_pop[p];
        pop = this->populations[p];
        int n = indivs.size();
        double exp_num_off = n * pop->mean_w * base_fec;

        int n_off = int(exp_num_off);

        for (int i = 0; i < n_off; i++ ){
            random1 = this->get_random_individual(p);
            random2 = random1;

            while (random2 == random1){
                random2 = this->get_random_individual(p);
            }

            parent_migrated = (random1->has_migrated || random2->has_migrated);

            int parent1_home = random1->pop_born_in;
            int parent2_home = random2->pop_born_in;

            offspring = new Individual(this, p, parent1_home, parent2_home, parent_migrated);
            offspring->gen_haplotype(random1, 0);
            offspring->gen_haplotype(random2, 1);
            this->next_gen->push_back(offspring);
        }

    }


    // Replace current generation with the next
    for (Individual* indiv : *(this->indivs)){
        delete indiv;
    }
    delete this->indivs;
    //delete[] this->indivs;
    this->indivs = next_gen;
    this->next_gen = new std::vector<Individual*>;
}


// ================================================================
// Helper
// functions
// ================================================================
void MPFM_Instance::split_indivs_by_pop(){
    int n_pops = this->populations.size();

    this->indivs_by_pop.clear();

    if (indivs_by_pop.size() == 0){
        for (int p = 0; p < n_pops; p++){
            this->indivs_by_pop.push_back(std::vector<Individual*>());
        }
    }


    for (Individual* indiv : *(this->indivs)){
        int pop = indiv->current_pop;
        if (pop >= 0){
            this->indivs_by_pop[pop].push_back(indiv);
        }
    }
}

std::vector<Individual*> MPFM_Instance::pull_n_random_indivs(int pop_number, int num_indivs){
    int ct = 0;
    std::vector<Individual*> random_indivs;


    std::vector<Individual*> pool;
    for (Individual* indiv : this->indivs_by_pop[pop_number]){
        if (!indiv->has_migrated){
            pool.push_back(indiv);
        }
    }

    if (num_indivs > pool.size()){
        num_indivs = pool.size();
    }

    while (ct < num_indivs){
        int size = pool.size();
        if (size > 0){
            int index = this->uniform_int(0,size-1);
            Individual* indiv = pool[index];

            indiv->has_migrated = true;
            random_indivs.push_back(indiv);
            ct++;
            pool.erase(pool.begin()+index);
        }
        else {
            return random_indivs;
        }
    }

    return random_indivs;
}

int MPFM_Instance::get_pop_index_from_id(int id){
    int n_pops = this->populations.size();
    for (int p = 0; p < n_pops; p++){
        Population* pop = this->populations[p];
        if (pop->id == id){
            return(p);
        }
    }
    return(-1);
}

void MPFM_Instance::calc_indivs_fitnesses(int pop){
    std::vector<double> this_pops_efs = this->populations[pop]->efs;
    std::vector<Individual*> indivs = this->indivs_by_pop[pop];
    for (Individual* indiv : indivs){
        indiv->calc_fitness(this_pops_efs);
    }
}

void MPFM_Instance::run_beverton_holt_selection(int pop){
    std::vector<Individual*> indivs = this->indivs_by_pop[pop];
    double k = this->populations[pop]->k;
    int n = indivs.size();

    double fitness_sum = 0.0;
    int this_pop_ct = 0;

    for (Individual* indiv : indivs){
        double abs_fitness = indiv->w;
        double k_prime = double(k * abs_fitness);
        double prob = this->beverton_holt_prob(n, k_prime);

        bool surv = false;
        double draw = this->uniform_01();
        if (draw < prob){
            surv =  true;
        }

        if (!surv){
            // current_pop of -1 is the graveyard. rip
            indiv->current_pop = -1;
        }
        else {
            fitness_sum += abs_fitness;
            this_pop_ct++;
        }
    }


    double this_pop_mean_abs_fitness = double(fitness_sum)/double(this_pop_ct);
    this->populations[pop]->mean_w = this_pop_mean_abs_fitness;
    // this->exp_num_off_this_population = this_pop_mean_abs_fitness*K;
}

double MPFM_Instance::beverton_holt_prob(double n, double Kprime){
    int b = this->params["AVG_NUM_OFFSPRING_PER_INDIV"];
    double prop_full = double(n)/double(Kprime);
    double prob = double(1.0) / double(1 + (double(b)/double(2) - 1)*(prop_full));
    return prob;
}

Individual* MPFM_Instance::get_random_individual(int pop){
    std::vector<Individual*> indivs = this->indivs_by_pop[pop];

    int size = indivs.size();
    int index = this->uniform_int(0, size-1);
    return indivs[index];
}

std::vector<Individual*> MPFM_Instance::sample_n_random_indivs(int pop, int n){
    std::vector<Individual*> indivs = this->indivs_by_pop[pop];
    std::vector<Individual*> random_indivs;

    int max_index = indivs.size() - 1;

    int ct = 0;
    while (ct  < n){
        int index = this->uniform_int(0, max_index);
        bool new_indiv = true;
        for (Individual* r_indiv : random_indivs){
            if (r_indiv == indivs[index]){
                new_indiv = false;
            }
        }

        if (new_indiv){
            random_indivs.push_back(indivs[index]);
            ct++;
        }
    }

    return random_indivs;
}


// ================================================================
// Methods for
// initialization
// ================================================================
void MPFM_Instance::initialize(){
    this->read_params_ini();
    this->read_pops_ini();

    this->read_burnin_dis_kernel();
    this->read_fragmentation_dis_kernel();

    this->dis_kernel = this->burn_in_dis_kernel;

    this->read_genome_dict();
    this->init_random_generators();
    this->init_individuals();
    this->init_genomes();

    this->data_wrangler = new DataWrangler(this);

}

void MPFM_Instance::read_params_ini(){
    std::ifstream infile("params.ini");
    std::string line;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }

        std::string name = record[0];
        double val = atof(record[1].c_str());

        this->params.insert(std::pair<std::string, double>(name, val));
    }
}

void MPFM_Instance::read_pops_ini(){
    std::ifstream infile("pops.ini");
    std::string line;

    // order of header
    // x,y,k,ef0,ef1,ef2...

    int n_ef = this->params["EF_NUMBER"];
    int ct = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }
        int record_should_be = 4 + n_ef;

        if (record.size() != record_should_be){
            printf("population file error on line: %d\nwrong number of entries.\n\n", ct);
            exit(-1);
        }


        int num = atoi(record[0].c_str());
        double x = atof(record[1].c_str());
        double y = atof(record[2].c_str());
        double k = atof(record[3].c_str());
        std::vector<double> efs;

        for (int ef = 0; ef < n_ef; ef++){
            int index = ef + 4;
            double ef_val = atof(record[index].c_str());
            efs.push_back(ef_val);
        }

        // construct population
        Population* tmp = new Population(this, num, ct, x, y, k, efs);
        this->populations.push_back(tmp);
        ct++;
    }
}


void MPFM_Instance::read_burnin_dis_kernel(){
    read_dis_kernel("diskern_pre.ini");
    this->burn_in_dis_kernel = this->dis_kernel;

}
void MPFM_Instance::read_fragmentation_dis_kernel(){
    read_dis_kernel("diskern_post.ini");
    this->fragmentation_dis_kernel = this->dis_kernel;
}

void MPFM_Instance::read_dis_kernel(std::string file){
    std::ifstream infile(file.c_str());
    std::string line;

    int n_pops = this->params["N_POPULATIONS"];

    for (int p = 0; p < n_pops; p++){
        this->dis_kernel.push_back(std::vector<double>());
        for (int p2 = 0; p2 < n_pops; p2++){
            this->dis_kernel[p].push_back(0.0);
        }
    }

    // file format:
    // pop1,pop2,val

    // get id_val from pops file and match it to id in ld file
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }

        int id1 = atoi(record[0].c_str());
        int id2 = atoi(record[1].c_str());


        int i = this->get_pop_index_from_id(id1);
        int j = this->get_pop_index_from_id(id2);
        double val = atof(record[2].c_str());

        this->dis_kernel[i][j] = val;
    }
}

void MPFM_Instance::read_genome_dict(){
    std::ifstream infile("genome.ini");
    std::string line;

    int n_ef = this->params["EF_NUMBER"];
    int loci_per_ef = this->params["N_LOCI_PER_EF"];
    int n_neutral_loci = this->params["N_NEUTRAL_LOCI"];

    int n_loci = n_ef*loci_per_ef + n_neutral_loci;

    this->N_LOCI = n_loci;

    for (int l = 0; l < n_loci; l++){
        this->selection_strengths.push_back(0.0);
        this->chromosome_map.push_back(0);
        this->map_dists.push_back(0.0);
        this->ef_genome_map.push_back(0);
        this->init_polymorphism_ct.push_back(0);
    }
    // file format:
    // locus,selection_str,ef,chromosome,map_distance
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }

        int l = atoi(record[0].c_str());
        double sel_str = atof(record[1].c_str());
        int ef = atoi(record[2].c_str());
        int chromo = atoi(record[3].c_str());
        double map_dist = atof(record[4].c_str());
        int init_poly= atoi(record[5].c_str());

        this->selection_strengths[l] = sel_str;
        this->ef_genome_map[l] = ef;
        this->chromosome_map[l] = chromo;
        this->map_dists[l] = map_dist;
        this->init_polymorphism_ct[l] = init_poly;
    }
}

void MPFM_Instance::init_random_generators(){
    int rs = this->params["RS_MAIN"];
    this->main_gen = new std::mt19937(rs);
}
void MPFM_Instance::init_individuals(){
    //this->indivs
    this->indivs = new std::vector<Individual*>;
    this->next_gen = new std::vector<Individual*>;

    int n_pops = this->params["N_POPULATIONS"];
    for (int p = 0; p < n_pops; p++){
        Population* pop = this->populations[p];
        int k = pop->k;

        for (int i = 0; i < k; i++){
            Individual* tmp = new Individual(this, p, -1, -1, 0);
            this->indivs->push_back(tmp);
        }
    }
}

void MPFM_Instance::init_genomes(){
    int n_loci = this->N_LOCI;
    std::vector<std::vector<double>> alleles;
    for (int l = 0; l < n_loci; l++){
        int n_alleles = this->init_polymorphism_ct[l];
        if (n_alleles == 0){
            n_alleles = 1;
        }
        alleles.push_back(std::vector<double>());
        for (int al = 0;  al < n_alleles; al++){
            double al_val = this->uniform_01();
            alleles[l].push_back(al_val);
        }
    }

    int allele_ct, rand_allele_index;
    double val;

    for (Individual* indiv : *(this->indivs)){
        for (int l = 0; l < n_loci; l++){
            allele_ct = alleles[l].size();
            rand_allele_index = this->uniform_int(0, allele_ct-1);
            val = alleles[l][rand_allele_index];
            indiv->haplotype0[l] = val;

            rand_allele_index = this->uniform_int(0, allele_ct-1);
            val = alleles[l][rand_allele_index];
            indiv->haplotype1[l] = val;
        }
    }
}

// ================================================================
// Random generation
// ================================================================
double MPFM_Instance::uniform_01(){
    std::uniform_real_distribution<double> dis(0.0,1.0);
    double r = dis(*this->main_gen);
    return r;
}

int MPFM_Instance::uniform_int(int lo, int hi){
    std::uniform_int_distribution<int> dis(lo, hi);
    return(dis(*this->main_gen));
}
int MPFM_Instance::binomial(double p, int n){
    std::binomial_distribution<int> b_dis(n, p);
    return b_dis(*this->main_gen);
}
