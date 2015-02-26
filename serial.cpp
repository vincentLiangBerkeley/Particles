#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <stack>
#include <set>
#include "common.h"
using namespace std;


bool sanityCheckOfBins(bin_t *bins, int binSize, int n)
{
    // This is for insanity check
    int sum = 0;
    for (int i = 0; i < binSize; ++i){
        //printf("bin %d contains %d particles\n", i, bins[i].numParticles());
        sum += bins[i].numParticles();
    }

    // printf("The total number of particles in the bins is %d\n", sum);
    return sum == n;
}

void assignParticleToBin(bin_t *bins, particle_t *particles, int index, int numBins)
{
    double binLength = bins[0].x_length;
    int row = int(particles[index].x / binLength);
    if (row == numBins) row = numBins - 1;
    int col = int(particles[index].y / binLength);
    if (col == numBins) col = numBins - 1;
    //printf("Assigning particle %d at (%f, %f) to bin (%d, %d)\n", index, particles[index].x, particles[index].y, row, col);
    bins[row + col * numBins].addParticle(index);
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n ); // Calculate the size
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    // Here is a test for the "bin" struct
    bin_t* bins = initBins();
    int numBins = getBinNum();

    // We have to iterate through all the particles and put them in corresponding bins.
    for (int i = 0; i < n; ++i)
        assignParticleToBin(bins, particles, i, numBins);

    for (int step = 0; step < NSTEPS; ++step)
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        /*
        Compute forces, but notice only compute forces from neighboring bins
         */    
        for (int i = 0; i < numBins; ++i)
        {
            for (int j = 0; j < numBins; ++j)
            {
                //printf("Checking bin (%d, %d)\n", i, j);
                if(bins[i + j * numBins].isEmpty()) continue;
        		for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
        		{
        		    particles[*it].ax = 0;
        		    particles[*it].ay = 0;

                    // Forces from within the bin
                    applyForceFromBin(bins[i + j * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from left
                    if (j > 0) applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force From right
                    if (j < numBins - 1) applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from above
                    if (i < numBins - 1) applyForceFromBin(bins[i + 1 + j * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from below
                    if (i > 0) applyForceFromBin(bins[i - 1 + j * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from north-west
                    if ((i < numBins - 1) && (j > 0)) applyForceFromBin(bins[i + 1 + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from north-east
                    if ((i < numBins - 1) && (j < numBins - 1)) applyForceFromBin(bins[i + 1 + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from south-west
                    if ((i > 0) && (j > 0)) applyForceFromBin(bins[i - 1 + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg);

                    // Force from south-east
                    if ((i > 0) && (j < numBins - 1)) applyForceFromBin(bins[i - 1 + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg);

        		}
            }
        }
        for(int i = 0; i < n; i ++) 
            move(particles[i]);

        // Determine the new bins of the particles
        // Try just simply rebinning all the particles
        for (int i = 0; i < numBins; ++i){
            for (int j = 0; j < numBins; ++j){
                bins[i + j * numBins].particleIndices -> clear();
            }
        }
        for (int i = 0; i < n; ++i){
            assignParticleToBin(bins, particles, i, numBins);
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) 
          {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
        
          //    `
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
	
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    free( bins );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
