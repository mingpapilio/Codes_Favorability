/**********************************************************************************************************
 * The is the source code for the individual-based simulation in 'Benign environments promote thermal 
 * trait diversity'. In this code file, we set up simulation to model abiotic factor on the evolution of
 * functional trait space (how temperature distribution affect the communities thermal performance curves). 
 * Notes:
 *      1. Please make sure random number generator is located at the right directory (dSFMT)
 *      2. Switches (s0, s1, s2) are implemented to make changing simulation settings easier. In particular
 *         , s0 is for deterministic (ecological model only) or stochastic thermal performance curves; s1
 *         is for ecological model or evolutionary model; s2 is for selecting the hypothesis to test
 *      3. Feel free to email 'ming.liu.ac[at]gmail.com' if encounter problems with execution
 **********************************************************************************************************/
// Including crucial packages and MersenneTwister random number generator
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
// Global functions
    double Mean_array(double p[], int length);
    double SD_array(double p[], int length);
    double sum(double p[], int length);
    int RandFromProb(double p[], int length, double RandNum);
    double normal_dist_BM (double mean, double sd, double u1, double u2);
    double **d_array2d(long size_1, long size_2);
    void free_d_array2d(double **x);
    double ****d_array4d(long size_1, long size_2, long size_3, long size_4);
    void free_d_array4d(double ****x);
    double safe_add(double a, double b);
    double get_double(const char *str);
// Main code
int main( int argc, char *argv[] ){
    // Switches
        int s0= 1;          // 0: Deterministic performance curves, 1: asymmetric thermal performance curves generated from distributions (random), 2: 2-species case
        int s1= 1;          // 0: Ecological model, 1: Evolutionary model with migration of species
        int s2= 0;          // 0: No variation in the volume of species performance curve, 1: weaker trade-off, some species can be stronger than others
        int s3= 1;          // 0: Fecundity is independent with average temperature, 1: dependent with average temperature 
        int s4= 0;          // 0: White noise, 1: Sine wave
        if(s0==2) s1=0;
    // Execution parameters
        int T_limit= 5E3;   // Time limit for number of iterations
        if(s1==1) T_limit= 1E5;
        int N_rep= 50;      // Number of repeats for each environmental condition
        int T_birthlog= 100; // Time points of counting species birth rate (last xx time steops of the simulation)
    // Population parameters
        int N_patch= 200;       // Number of patches in the community
        int N_spcs= 50;         // Maximal number of species in the community
        int N_ini= 3;           // Initial population size per species per patch  
        int N_samp= 10;         // Number of sampling for producing offspring
        int N_spcs_ini= 5;      // Number of initial species, used in evolutionary model
        if(s1==0) N_spcs_ini= N_spcs;
        if(s0==2) N_spcs_ini= 2;
        int N_spcs_samp= 100;   // Number of sampling if a species has extinct, used in evolutionary model
        int N_patch_newsp= 10;  // Number of patch to initialise the new species
        //
        double mortality= 0.04;      // Mortality of each individual
        double prob_migrate= 0.001;  // Probability of getting new species migrating to the community, used in evolutionary model
        // Parameters for generating thermal performance curve
        double avg_topt= 20.0;  // Average of temperature of optimal thermal performance (Topt)
        double sd_topt= 7.5;    // SD of Topt
        double avg_tinc= 5.0;   // Average of Tinc, the distance between critical maximal temperature (CTmax) and Topt
        double sd_tinc= 2.5;    // SD of Tinc
        double avg_tsgm= 5.0;   // Average of sigma, decay coefficient for thermal performance below Topt
        double sd_tsgm= 2.5;    // SD of sigma
    // Temperaure distribution (environmental factors)
        int env_span= 1;        // Time scale of environmental episode (how frequent to fluctuate)
        int env_scale= 50;     // Relative time scale of long-term to short-term variaitons
        double tmp_env;         // Mean temperature of temperature distribution
        double curr_env;        // Current tempreature
        double past_env;        // Temperature of the last time-step, check for updating population growth parameters
        // Normal distribution of temperature
            double mean_env= 18.0;  // Average temperature
            double sd_env=	0.0;// Amplitude of long-term variation (SD of a normal distribution)
            double sd_env_short= 0.5;
            double u1, u2, mean_env_short;
        // Mean temperature settings for fecundity
            double env_max= 38.0;       // Critical temperature that fecundity is zero for all above temperature
            double env_opt= 32.0;       // The most-fecund temperature
            double tmean_sigma= 7.5;    // Decay coefficient for fecundity below optimal temperature
            double fecund_scale= 3.0;   // Scaling coefficient
            double fecundity= 2.0;      // 
            double fecund_intg, fecund_decm, fecund_stoc;
    //********** We do not recommend changing things below this line  **********//
    // Indices of species matrix
        int id_popn= 0;     // Population size of focal species
        int id_offs= 1;     // Number of offspring of focal species produced at current time step
        int id_topt= 2;     // Temperature of optimal thermal performance (Topt)
        int id_tinc= 3;     // Distance of critical maximal temperature (CTmax) from Topt
        int id_tsgm= 4;     // Decay coefficient for thermal performance below Topt
        int id_conv= 5;     // Thermal performance at current temperature
        int id_weig= 6;     // Current weighing coefficient for focal species at focal patch
        int id_thrs= 7;     // Threshold coefficient for density-dependent and density-independent reproduction
        int id_scal= 8;     // Scaling coefficient for thermal performance curve to ensure fair competition
        int id_init= 9;     // Time point of initializing the species
        int id_prnt= 10;    // Whether the species has been printed to the extinction log
        int id_coff= 11;    // Cumulative offspring counting
    // Variables
        int i,j,k,t,rep;
        int idx, off_size, temp_size, weigh_dens, prob_temp, samp_id, pop_temp;
        double pp;
        double Popn[N_patch][N_spcs];
        double PatchComp[N_spcs];
        double Spcs[N_spcs][12]; ///// Depends on how many properties to be logged for each species. /////
    // Initialisation of the random number genertor
        int seed;
        dsfmt_t dsfmt;
        seed= time(NULL);
        if(seed==0)seed= 1;
        dsfmt_init_gen_rand(&dsfmt,seed);
    // Initialising output file
        time_t time_start, time_end;
        time_start= time(NULL);
        FILE *surv;
        surv= fopen("species.txt", "w");
        if(s0==2) fprintf(surv,"Tmean\tTsd\tTsd_short\tOutcome\tRep\n");
        else fprintf(surv,"Tmean\tTsd\tTsd_short\tTmin\tTmax\tTopt\tInc\tSigma\tScale\tTpers\tFec\tDens\tBRate\tRep\n");
        FILE *log;
        log= fopen("extinction_log.txt", "w");
        fprintf(log, "Tmean\tTsd\tTsd_short\tTmin\tTmax\tTopt\tLifespan\tRep\n");
    // Read in variables
    if( argc == 3 ) {
        int i;
        // mean_env= get_double(argv[1]);
        sd_env_short= get_double(argv[1]);
        sd_env= get_double(argv[2]);
        // printf("Avg is %lf, Var is %lf.\n", mean_env, sd_env_short);
    }
    else {
        printf("Wrong number of argumets.\n");
    }
    // Main Code
    for(rep=0; rep<N_rep; ++rep){

    // Initialisation of the population
    for(i=0; i<N_patch; ++i){
        for(j=0; j<N_spcs; ++j){
            if(j<N_spcs_ini) Popn[i][j]= N_ini;
            else Popn[i][j]= 0.0;
        }
    }
    // Initialisation of each species /////// Add standardisation
    if(s0==2){
        Spcs[0][id_thrs]= 10.0;
        Spcs[0][id_topt]= 30.0;
        Spcs[0][id_tinc]= 5.0;
        Spcs[0][id_tsgm]= 5.0;
        Spcs[1][id_thrs]= 10.0;
        Spcs[1][id_topt]= 17.0;
        Spcs[1][id_tinc]= 6.0;
        Spcs[1][id_tsgm]= 2.0;
        for(i=0; i<2; ++i){
            Spcs[i][id_scal]= (2.0/3*5+5*sqrt(M_PI))/(2.0/3*Spcs[i][id_tinc]+Spcs[i][id_tsgm]*sqrt(M_PI));
        }
    }
    else{
        for(i=0; i<N_spcs_ini; ++i){
            if(s0==0){
                Spcs[i][id_thrs]= 10.0;
                Spcs[i][id_tsgm]= 5.0;
                Spcs[i][id_topt]= 10+i;
                Spcs[i][id_tinc]= 5.0;
                Spcs[i][id_scal]= 1.0;
                Spcs[i][id_init]= 0.0;
                Spcs[i][id_prnt]= 0.0;
                Spcs[i][id_coff]= 0.0;
            }
            if(s0==1){
                Spcs[i][id_scal]= -0.1;
                while((!(Spcs[i][id_scal]> 0.0)||!(Spcs[i][id_tinc]> 0.0))||!(Spcs[i][id_tsgm]> 0.0)){
                    Spcs[i][id_thrs]= 10.0;
                    Spcs[i][id_topt]= normal_dist_BM(avg_topt,sd_topt,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                    Spcs[i][id_tinc]= normal_dist_BM(avg_tinc,sd_tinc,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                    Spcs[i][id_tsgm]= normal_dist_BM(avg_tsgm,sd_tsgm,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                    Spcs[i][id_scal]= (2.0/3*5+5*sqrt(M_PI))/(2.0/3*Spcs[i][id_tinc]+Spcs[i][id_tsgm]*sqrt(M_PI));
                    if(s2==1) Spcs[i][id_scal]= Spcs[i][id_scal]*normal_dist_BM(1,0.25,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt)); ///// sd set to 0.25
                    Spcs[i][id_init]= 0.0;
                    Spcs[i][id_prnt]= 0.0;
                    Spcs[i][id_coff]= 0.0;
                }
            }
        }
        for(i=N_spcs_ini; i<N_spcs; ++i) Spcs[i][id_prnt]= 1.0;
    }

    // Simulation
    for(t=0; t<T_limit; ++t){
        // Environmental fluctuation
        if(t%env_span==0){
            if(t%(env_span*env_scale)== 0 && s4==0) mean_env_short= normal_dist_BM(mean_env,sd_env,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
            if(s4==1){
                mean_env_short= mean_env+ sd_env*sin(t*2*M_1_PI/env_span/360);
            }
            curr_env= normal_dist_BM(mean_env_short,sd_env_short,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
            //
            if(s3==1){
                if(curr_env> env_max) fecundity= 0.0;
                else{
                    if(curr_env> env_opt) fecundity= fecund_scale*(1-((curr_env-env_opt)*(curr_env-env_opt)/(env_max-env_opt)/(env_max-env_opt)));
                    else fecundity= fecund_scale*exp(-1*((curr_env-env_opt)/2/tmean_sigma)*((curr_env-env_opt)/2/tmean_sigma));
                }
            }
            //
            for(i=0; i<N_spcs; ++i){
                if(curr_env> Spcs[i][id_topt]) Spcs[i][id_conv]= Spcs[i][id_scal]*(1- (curr_env-Spcs[i][id_topt])*(curr_env-Spcs[i][id_topt])/Spcs[i][id_tinc]/Spcs[i][id_tinc]);
                else Spcs[i][id_conv]= Spcs[i][id_scal]*exp(-1*((curr_env-Spcs[i][id_topt])/2/Spcs[i][id_tsgm])*((curr_env-Spcs[i][id_topt])/2/Spcs[i][id_tsgm]));
            }
        }
        // Reproduction
        for(i=0; i<N_patch;++i){
            // sum the weighted density
            weigh_dens= 0.0;
            for(j=0; j<N_spcs; ++j){
                Spcs[j][id_weig]= Popn[i][j]/Spcs[j][id_thrs];
                if(Spcs[j][id_weig]> 1.0) Spcs[j][id_weig]= 1.0;
                weigh_dens+= Spcs[j][id_weig];
            }
            // Ability for competition
            for(j=0; j<N_spcs; ++j) PatchComp[j]= Popn[i][j]/Spcs[j][id_thrs]*Spcs[j][id_conv];
            // Type 2 functional response for competition/ no competition
                fecund_decm= modf(fecundity, &fecund_intg);
                if(dsfmt_genrand_open_close(&dsfmt)< fecund_decm) fecund_stoc= fecund_intg;
                else fecund_stoc= fecund_intg+ 1;
            // 
            for(j=0; j<fecundity; ++j){
                if(dsfmt_genrand_open_close(&dsfmt)< (weigh_dens/(1+weigh_dens))){
                    idx= RandFromProb(PatchComp,N_spcs,dsfmt_genrand_open_close(&dsfmt));
                    if(Popn[i][idx]> 0.0) Spcs[idx][id_offs]+= 1.0;
                }
                else{
                    //
                    for(k=0; k<N_samp; ++k){
                        idx= floor(dsfmt_genrand_open_close(&dsfmt)*N_spcs);
                        if(dsfmt_genrand_open_close(&dsfmt)< Spcs[idx][id_conv] && Popn[i][idx]> 0.0){
                            Spcs[idx][id_offs]+= 1.0;
                            break;
                        }
                    }
                }
            }
        }
        // Mortality
        for(i=0; i<N_patch;++i){
            for(j=0;j<N_spcs;++j){
                // Popn[i][j]= round(Popn[i][j]*(1-mortality));
                temp_size= (int)Popn[i][j];
                if(temp_size>0){
                    for(k=0; k<temp_size; ++k){
                        pp= dsfmt_genrand_open_close(&dsfmt);
                        if(pp< mortality) Popn[i][j]-= 1.0;
                }}
                if(Popn[i][j]< 0.0) Popn[i][j]= 0.0;
            }
        }
        // Allocating the offspring
        for(i=0; i<N_spcs; ++i){
            off_size= (int)round(Spcs[i][id_offs]);
            // Adding the number of offspring for all species (from T_birthlog to T_limit-1)
            if(t>=(T_limit-T_birthlog)) Spcs[i][id_coff]= safe_add(Spcs[i][id_coff],(double)off_size);
            for(j=0; j<off_size; ++j){
                idx= (int)floor(N_patch*dsfmt_genrand_open_close(&dsfmt));
                Popn[idx][i]+= 1.0;
                Spcs[i][id_offs]-= 1.0;
            }
            Spcs[i][id_offs]= 0.0;
        }
        // Log the extinction events
        // fprintf(log, "Tmin\tTmax\tTopt\tLifespan\n");
        for(i=0; i<N_spcs; ++i){
            if(Spcs[i][id_prnt]== 0.0){
                Spcs[i][id_popn]= 0.0;
                for(j=0;j<N_patch; ++j){
                    Spcs[i][id_popn]+= Popn[j][i];
                    if(Spcs[i][id_popn]> 0.0) break;
                }
                if(Spcs[i][id_popn]== 0.0 && Spcs[i][id_offs]== 0.0){
                    fprintf(log,"%lf\t", mean_env);
                    fprintf(log,"%lf\t", sd_env);
                    fprintf(log,"%lf\t", sd_env_short);
                    fprintf(log,"%lf\t",(Spcs[i][id_topt]-Spcs[i][id_tsgm]*3.46164));
                    fprintf(log,"%lf\t",(Spcs[i][id_topt]+Spcs[i][id_tinc]*0.974679));
                    fprintf(log,"%lf\t", Spcs[i][id_topt]);
                    fprintf(log,"%lf\t",(t-Spcs[i][id_init]));
                    fprintf(log,"%d\n", rep);
                    Spcs[i][id_prnt]= 1.0;
                }
            }
        }
    
        // Immigration of new species through allopatric speciation+ migration
        if(dsfmt_genrand_open_close(&dsfmt)< prob_migrate && s1==1){
            for(i=0; i<N_spcs_samp; ++i){
                samp_id= floor(dsfmt_genrand_open_close(&dsfmt)*N_spcs);
                pop_temp= 0;
                for(j=0; j<N_patch; ++j) pop_temp+= Popn[j][samp_id];
                if(pop_temp==0){
                    // Define the thermal performance curve
                    Spcs[samp_id][id_scal]= -0.1;
                    while((!(Spcs[samp_id][id_scal]> 0.0)||!(Spcs[samp_id][id_tinc]> 0.0))||!(Spcs[samp_id][id_tsgm]> 0.0)){
                        Spcs[samp_id][id_thrs]= 10.0;
                        Spcs[samp_id][id_topt]= normal_dist_BM(avg_topt,sd_topt,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                        Spcs[samp_id][id_tinc]= normal_dist_BM(avg_tinc,sd_tinc,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                        Spcs[samp_id][id_tsgm]= normal_dist_BM(avg_tsgm,sd_tsgm,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                        Spcs[samp_id][id_scal]= (2.0/3*5+5*sqrt(M_PI))/(2.0/3*Spcs[samp_id][id_tinc]+Spcs[samp_id][id_tsgm]*sqrt(M_PI));
                        if(s2==1) Spcs[samp_id][id_scal]= Spcs[samp_id][id_scal]*normal_dist_BM(1,0.1,dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt)); ///// sd set to 0.1
                        Spcs[samp_id][id_init]= (double)t;
                        Spcs[samp_id][id_prnt]= 0.0;
                        Spcs[samp_id][id_coff]= 0.0;
                    }
                    for(j=0; j< N_patch_newsp; ++j){
                        idx= floor(dsfmt_genrand_open_close(&dsfmt)*N_patch);
                        Popn[idx][samp_id]= N_ini;
                    }
                    break;
                }
            }
        }
        // Print
        if(t== (T_limit-1)){
            for(i=0; i<N_spcs; ++i){
                Spcs[i][id_popn]= 0.0;
                for(j=0;j<N_patch; ++j){
                    Spcs[i][id_popn]+= Popn[j][i];
                }
            }
            if(s0==2){
                fprintf(surv,"%lf\t", mean_env);
                fprintf(surv,"%lf\t", sd_env);
                fprintf(surv,"%lf\t", sd_env_short);
                if(Spcs[0][id_popn]> 0){
                    if(Spcs[1][id_popn]> 0) fprintf(surv,"Coexist\t");
                    else fprintf(surv, "Spcs1\t");
                }
                else{
                    if(Spcs[1][id_popn]> 0) fprintf(surv,"Spcs2\t");
                    else fprintf(surv,"Extinct\t");
                }
                fprintf(surv,"%d\n", rep);
            }
            else{
                for(i=0; i<N_spcs; ++i){
                    if(Spcs[i][id_popn]> 0){
                        fprintf(surv,"%lf\t", mean_env);
                        fprintf(surv,"%lf\t", sd_env);
                        fprintf(surv,"%lf\t", sd_env_short);
                        fprintf(surv,"%lf\t",(Spcs[i][id_topt]-Spcs[i][id_tsgm]*3.46164));
                        fprintf(surv,"%lf\t",(Spcs[i][id_topt]+Spcs[i][id_tinc]*0.974679));
                        fprintf(surv,"%lf\t", Spcs[i][id_topt]);
                        fprintf(surv,"%lf\t", Spcs[i][id_tinc]);
                        fprintf(surv,"%lf\t", Spcs[i][id_tsgm]);
                        fprintf(surv,"%lf\t", Spcs[i][id_scal]);
                        fprintf(surv,"%lf\t", t- Spcs[i][id_init]);
                        fprintf(surv,"%lf\t", fecundity);
                        fprintf(surv,"%lf\t", Spcs[i][id_popn]);
                        fprintf(surv,"%lf\t", (Spcs[i][id_coff]/T_birthlog));
                        fprintf(surv,"%d\n", rep);
                    }
                }
            }
        }
    }}
    // Closing file
    time_end= time(NULL);
    printf("This simulation lasted %lf minutes.\n", (time_end-time_start)/ 60.0);
    fclose(surv);
    // End
    return EXIT_SUCCESS;   
}

////////////////////// Functions are defined in below //////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp;
    temp= 0.0;
    if(length> 0){
        for (i=0; i<length; ++i) temp+= p[i]/length;
        return temp;
    }
    else return NAN;
}

double SD_array(double p[], int length){
    if(length> 1){
        int i;
        double avg, temp;
        avg= temp= 0.0;
        for (i=0; i<length; ++i) avg+= p[i]/length;
        for (i=0; i<length; ++i) temp+= pow(p[i]-avg, 2);
        return sqrt(temp/(length-1));
    }
    else return NAN;
}

double sum(double p[], int length){
    int i;
    double sum=0.0;
    for(i=0; i< length; ++i) sum= safe_add(sum, p[i]);
    return sum;
}

int RandFromProb(double p[], int length, double RandNum){
    int i,temp,check;
    double sum, pp;
    sum= pp= 0.0;
    check= temp= 0;
    for(i=0; i< length; ++i){
        if(p[i]> 0.0) sum= safe_add(sum, p[i]/length);
    }
    if(sum> 0.0){
        for(i=0; i< length; ++i){
            if(RandNum> pp/sum) check=1;
            if(p[i]> 0.0) pp+= p[i]/length;
            if(RandNum< pp/sum&& check==1){
                temp= i;
                break;
    }}}
    else{
        temp= floor(length*RandNum);
    }
    return temp;
}

double normal_dist_BM (double mean, double sd, double u1, double u2){
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}

double **d_array2d(long size_1, long size_2){
    double **x;
    long i;
    x= (double **) malloc((size_t)(size_1*sizeof(double)));
    for(i=0; i< size_1; i++) x[i]= (double *) malloc((size_t)(size_2*sizeof(double)));
    return x;
}
void free_d_array2d(double **x){
    free(x[0]);
    free(x);
}

double ****d_array4d(long size_1, long size_2, long size_3, long size_4){
    double ****x;
    long i,j,k;
    x= (double ****) malloc((size_t)(size_1*sizeof(double ***)));
    for(i=0; i< size_1; i++) {
        x[i]= (double ***) malloc((size_t)(size_2*sizeof(double **)));
        for(j=0; j< size_2; j++){
            x[i][j]= (double **) malloc((size_t)(size_3*sizeof(double *)));
            for(k=0; k< size_3; k++) x[i][j][k]= (double *) malloc((size_t)(size_4*sizeof(double)));
        }
    }
    return x;
}

void free_d_array4d(double ****x){
    free(x[0][0][0]);
    free(x[0][0]);
    free(x[0]);
    free(x);
}

double safe_add(double a, double b){
    if (a> 0 && b> DBL_MAX- a) printf("Overflow in safe_add.\n");
    else if (a< 0 && b< - DBL_MAX- a) printf("Underflow in safe_add.\n");
    // DBL_MIN is the smallest (closest value to 0) floating number
    return a+b;
}

double get_double(const char *str)
{
    /* First skip non-digit characters */
    /* Special case to handle negative numbers and the `+` sign */
    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
        str++;

    /* The parse to a double */
    return strtod(str, NULL);
}
