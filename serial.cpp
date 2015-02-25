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
    
    // printf("The number of bins in a row is %d\n", numBins);    
    // for (int i = 0; i < numBins; ++i)
    // {
    //     printf("%d\n", i);
    //     for (int j = 0; j < numBins; ++j)
    //     {
    //         printf("The (%d, %d)th bin has corner: (%f, %f), and with xlength: %f, ylength: %f\n", i, j, bins[i + j * numBins].x, bins[i + j * numBins].y
    //             , bins[i + j * numBins].x_length, bins[i + j * numBins].y_length);
    //     }
    // }

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
                if(bins[i + j * numBins].isEmpty()) 
                {
                    //printf("bin (%d, %d) is empty\n", i, j);
                    continue;
                }
                // Apply force from within the bin
                for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                {
                    applyForceFromBin(bins[i + j * numBins], *it, particles, &dmin, &davg, &navg);
                    // printf("within bin absmin = %f\n", absmin);
                }
                // Iterating through all the bins, handle the edge cases first
                if (i == 0)
                {
                    if (j == 0)
                    {
                        // Left bottom corner, apply force from only three neighbors
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[1], *it, particles, &dmin, &davg, &navg);
                            // printf("1 left bottom corner absmin = %f\n", absmin);
                            // printf("1 left bottom corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[numBins], *it, particles, &dmin, &davg, &navg);
                            // printf("2 left bottom corner absmin = %f\n", absmin);
                            // printf("2 left bottom corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[1 + numBins], *it, particles, &dmin, &davg, &navg);
                            // printf("3 left bottom corner absmin = %f\n", absmin);
                            // printf("3 left bottom corner dmin = %f\n", dmin);
                        }
                        printf("1. left bottom corner absmin = %f\n", absmin);
                        printf("1. left bottom corner dmin = %f\n", dmin);
                    }
                    else if (j == numBins - 1)    
                    {
                        // Right bottom corner
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[i + j * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin above
                            // printf("1 right bottom corner absmin = %f\n", absmin);
                            // printf("1 right bottom corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + j * numBins - numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                            // printf("2 right bottom corner absmin = %f\n", absmin);
                            // printf("2 right bottom corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + j * numBins - numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the left top corner
                            // printf("3 right bottom corner absmin = %f\n", absmin);
                            // printf("3 right bottom corner dmin = %f\n", dmin);
                        }
                        printf("2. right bottom corner absmin = %f\n", absmin);
                        printf("2. right bottom corner dmin = %f\n", dmin);
                    }
                    else
                    {
                        // Just normal edge row, have five neighbors
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                            // printf("1 i=0 normal edge row absmin = %f\n", absmin);
                            // printf("1 i=0 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the right
                            // printf("2 i=0 normal edge row absmin = %f\n", absmin);
                            // printf("2 i=0 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + j * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on top of it
                            // printf("3 i=0 normal edge row absmin = %f\n", absmin);
                            // printf("3 i=0 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j - 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the top left corner
                            // printf("4 i=0 normal edge row absmin = %f\n", absmin);
                            // printf("4 i=0 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j + 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the top right corner
                            // printf("5 i=0 normal edge row absmin = %f\n", absmin);
                            // printf("5 i=0 normal edge row dmin = %f\n", dmin);
                        }
                        printf("3. i = 0 case normal edge row absmin = %f\n", absmin);
                        printf("3. i = 0 case normal edge row dmin = %f\n", dmin);
                    }
                }
                else if (i == numBins - 1)
                {
                    if (j == 0)
                    {
                        // Left top corner, three neighbors
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the right
                            // printf("1 left top corner absmin = %f\n", absmin);
                            // printf("1 left top corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i - j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                            // printf("2 left top corner absmin = %f\n", absmin);
                            // printf("2 left top corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j + 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the bottom right corner
                            // printf("3 left top corner absmin = %f\n", absmin);
                            // printf("3 left top corner dmin = %f\n", dmin);
                        }
                        printf("4. left top corner absmin = %f\n", absmin);
                        printf("4. left top corner dmin = %f\n", dmin);
                    }
                    else if (j == numBins - 1)
                    {
                        // Right top corner, three neighbors
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                            // printf("1 right top corner absmin = %f\n", absmin);
                            // printf("1 right top corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                            // printf("2 right top corner absmin = %f\n", absmin);
                            // printf("2 right top corner dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j - 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the left bottom corner
                            // printf("3 right top corner absmin = %f\n", absmin);
                            // printf("3 right top corner dmin = %f\n", dmin);
                        }
                        printf("5. right top corner absmin = %f\n", absmin);
                        printf("5. right top corner dmin = %f\n", dmin);
                    }
                    else
                    {
                        // Just normal edge row, have five neighbors
                        for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                        {
                            applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                            // printf("1 i=numBins-1 normal edge row absmin = %f\n", absmin);
                            // printf("1 i=numBins-1 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the right
                            // printf("2 i=numBins-1 normal edge row absmin = %f\n", absmin);
                            // printf("2 i=numBins-1 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                            // printf("3 i=numBins-1 normal edge row absmin = %f\n", absmin);
                            // printf("3 i=numBins-1 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j - 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the left bottom corner
                            // printf("4 i=numBins-1 normal edge row absmin = %f\n", absmin);
                            // printf("4 i=numBins-1 normal edge row dmin = %f\n", dmin);
                            applyForceFromBin(bins[i + (j + 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the right bottom corner
                            // printf("5 i=numBins-1 normal edge row absmin = %f\n", absmin);
                            // printf("5 i=numBins-1 normal edge row dmin = %f\n", dmin);
                        }
                        printf("6. i=numBins-1 normal edge row absmin = %f\n", absmin);
                        printf("6. i=numBins-1 normal edge row dmin = %f\n", dmin);
                    }
                }
                else if (j == 0 && i > 0 && i < (numBins - 1))
                {
                    for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                    {
                        applyForceFromBin(bins[i + j * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin above
                        // printf("1 j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("1 j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                        // printf("2 j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("2 j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the right
                        // printf("3 j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("3 j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the right top corner
                        // printf("4 j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("4 j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the right bottom corner
                        // printf("5 j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("5 j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                    }
                    printf("7. j == 0 && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                    printf("7. j == 0 && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                }
                else if (j == (numBins - 1) && i > 0 && i < (numBins - 1))
                {
                    for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                    {
                        applyForceFromBin(bins[i + j * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin above
                        // printf("1 j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("1 j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                        // printf("2 j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("2 j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                        // printf("3 j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("3 j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j - 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the left top corner
                        // printf("4 j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("4 j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j - 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the left bottom corner
                        // printf("5 j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                        // printf("5 j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                    }
                    printf("8. j == (numBins - 1) && i > 0 && i < (numBins - 1) absmin = %f\n", absmin);
                    printf("8. j == (numBins - 1) && i > 0 && i < (numBins - 1) dmin = %f\n", dmin);
                }
                else 
                {
                    // Now just inner bins, have nine neighbors 
                    for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); ++it)
                    {
                        applyForceFromBin(bins[i + (j - 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the left
                        // printf("1 inner bins absmin = %f\n", absmin);
                        // printf("1 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins], *it, particles, &dmin, &davg, &navg); // Bin on the right
                        // printf("2 inner bins absmin = %f\n", absmin);
                        // printf("2 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + j * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin above
                        // printf("3 inner bins absmin = %f\n", absmin);
                        // printf("3 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + j * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin below
                        // printf("4 inner bins absmin = %f\n", absmin);
                        // printf("4 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j - 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the left top corner
                        // printf("5 inner bins absmin = %f\n", absmin);
                        // printf("5 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins + 1], *it, particles, &dmin, &davg, &navg); // Bin on the right top corner
                        // printf("6 inner bins absmin = %f\n", absmin);
                        // printf("6 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j - 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the left bottom corner
                        // printf("7 inner bins absmin = %f\n", absmin);
                        // printf("7 inner bins dmin = %f\n", dmin);
                        applyForceFromBin(bins[i + (j + 1) * numBins - 1], *it, particles, &dmin, &davg, &navg); // Bin on the right bottom corner
                        // printf("8 inner bins absmin = %f\n", absmin);
                        // printf("8 inner bins dmin = %f\n", dmin);
                    }
                    printf("9. inner bins absmin = %f\n", absmin);
                    printf("9. inner bins dmin = %f\n", dmin);
                }

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
        }

        printf("step = %d\n", step);

        printf("after force applied absmin = %f\n", absmin);

        for(int i = 0; i < n; i ++) 
        {
           move(particles[i]); // Move particles
           //printf("particle %d is moving at a speed %f\n", i, particles[i].vx);
       }
        // Determine the new bins of the particles
        std::stack<int> movedParticles;
        for (int i = 0; i < numBins; ++i)
        {
            for (int j = 0; j < numBins; ++j)
            {
                if (!bins[i + j * numBins].isEmpty())
                {
                    //printf("Bin (%d, %d) has %d particles.\n", i, j, bins[i + j * numBins].numParticles());
                    for(std::set<int>::iterator it = bins[i + j * numBins].particleIndices -> begin(); it != bins[i + j * numBins].particleIndices -> end(); )
                    {
                        if (bins[i + j * numBins].outOfBound(particles[*it])) 
                        {
                            //printf("particle %d is out of bound\n", *it);
                            movedParticles.push(*it);
                            bins[i + j * numBins].particleIndices -> erase(it++); // If the particle is out of bound, delete it from bin
                        }
                        else ++it;
                    }
                }
            }
        }
        // Reassign the moved particles to their new bins
        int count = 0;
        while(!movedParticles.empty())
        {
            count ++;
            assignParticleToBin(bins, particles, movedParticles.top(), numBins);
            movedParticles.pop();
        }
        sanityCheckOfBins(bins, numBins * numBins, n);
        //printf("The number of moved particles is: %d\n", count);
        //printf("dmin = %f, davg = %f\n", dmin, davg);

        // if( find_option( argc, argv, "-no" ) == -1 )
        // {
        //   //
        //   // Computing statistical data
        //   //
        //   if (navg) {
        //     absavg +=  davg/navg;
        //     nabsavg++;
        //   }
        //   if (dmin < absmin) absmin = dmin;
        
        //   //    `
        //   //  save if necessary
        //   //
        //   if( fsave && (step%SAVEFREQ) == 0 )
        //       save( fsave, n, particles );
        // }
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
