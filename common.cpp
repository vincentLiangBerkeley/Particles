#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;
// Define the extern variables
double bin_size;
int bins_per_side;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//  This will also set constants for bin_t
//
void set_size( int n )
{
    size = sqrt( density * n );

    bin_size = sqrt(5 * density); // This defines the extern double, and we define it this way so that each bin contains about 5 particles
    bins_per_side = floor(size / bin_size) + 1;

    assert(bin_size > 2 * cutoff);
    printf("bin_size = %f, bins_per_side = %d, total num of bins = %d\n", bin_size, bins_per_side, bins_per_side * bins_per_side);
}

bin_t::bin_t()
    :numParticles(0) 
    {}

// This function just clears all the bins, but needs the memory of the bins to be allocated elsewhere
void init_bins(bin_t *bins)
{
    for (int i = 0; i < NUM_BINS; ++i){
        bins[i].numParticles = 0;
        bins[i].particles.clear();
    }
}

void assign_particles_to_bin(particle_t &p, bin_t *bins)
{
    int row = floor(p.x / BIN_SIZE);
    int col = floor(p.y / BIN_SIZE);

    int index = row + col * BINS_PER_SIDE;
    bins[index].particles.push_back(&p);
    bins[index].numParticles ++;
}

// This compute all forces applied to particles in the (row, col)th bin
void compute_forces_for_bin(bin_t *bins, int row, int col, double *dmin, double *davg, int *navg)
{
    int index = row + col * BINS_PER_SIDE;

    for (int p = 0; p < bins[index].numParticles; ++p)
    {
        (*bins[index].particles[p]).ax = 0;
        (*bins[index].particles[p]).ay = 0;

    #define FORCE_FROM_BIN(i) \
        for (int p2 = 0; p2 < bins[i].numParticles; ++p2) \
            apply_force(*bins[index].particles[p], *bins[i].particles[p2], dmin, davg, navg);
        
        // Forces from inside the bin
        FORCE_FROM_BIN(index);

        // Forces from left
        if (col > 0) FORCE_FROM_BIN(index - BINS_PER_SIDE);

        // Forces from right
        if (col < BINS_PER_SIDE - 1) FORCE_FROM_BIN(index + BINS_PER_SIDE);

        // Forces from above
        if (row < BINS_PER_SIDE - 1) FORCE_FROM_BIN(index + 1);

        // Forces from below
        if (row > 0) FORCE_FROM_BIN(index - 1);

        // Forces from north east
        if (row < BINS_PER_SIDE - 1 && col > 0) FORCE_FROM_BIN(index + 1 - BINS_PER_SIDE);

        // Forces from north west
        if (row < BINS_PER_SIDE - 1 && col < BINS_PER_SIDE - 1) FORCE_FROM_BIN(index + 1 + BINS_PER_SIDE);

        // Forces from south east
        if (row > 0 && col > 0) FORCE_FROM_BIN(index - 1 - BINS_PER_SIDE);

        // Forces from south west
        if (row > 0 && col < BINS_PER_SIDE - 1) FORCE_FROM_BIN(index - 1 + BINS_PER_SIDE);
    #undef FORCE_FROM_BIN    
    }
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
        {
       if (r2/(cutoff*cutoff) < *dmin * (*dmin))
          *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
        
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
    
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
