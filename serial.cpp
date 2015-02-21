#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
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

    printf("The total number of particles in the bins is %d\n", sum);
    return sum == n;
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
    double binLength = bins[0].x_length;

    printf("The number of bins in a row is %d\n", numBins);    
    if (numBins < 0) printf("true\n");
    for (int i = 0; i < numBins; ++i)
    {
        printf("%d\n", i);
        for (int j = 0; j < numBins; ++j)
        {
            printf("The (%d, %d)th bin has corner: (%f, %f), and with xlength: %f, ylength: %f\n", i, j, bins[i + j * numBins].x, bins[i + j * numBins].y
                , bins[i + j * numBins].x_length, bins[i + j * numBins].y_length);
        }
    }

    // We have to iterate through all the particles and put them in corresponding bins.
    for (int i = 0; i < n; ++i)
    {
        int row = int(particles[i].x / binLength);
        int col = int(particles[i].y / binLength);
        bins[row + col * numBins].addParticle(particles[i]);
    }

    //sanityCheckOfBins(bins, numBins * numBins, n);
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
    davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                // This line of code checks every pair of particles and applies the force
                // We need to use a mesh so that only neibors will have forces applied on a particle
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
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
    if( fsave )
        fclose( fsave );
    
    return 0;
}
