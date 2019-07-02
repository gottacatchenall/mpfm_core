#include "include.h"

#include "MPFM_Instance.h"

int main(int argc, char* argv[]){
    if (argv[1]){
        chdir(argv[1]);
    }

    MPFM_Instance* mpfm = new MPFM_Instance();

    mpfm->initialize();
    printf("Initialization completed...");

    mpfm->burn_in();
    mpfm->fragmentation();

    return(0);
}
