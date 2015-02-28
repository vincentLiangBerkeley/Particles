#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

/*
 * bin data structure
 */
struct bin_t 
{
  std::vector<particle_t*> particles; // Stores pointers to particles

  int numParticles;
  bin_t();
};

// The size of each bin
extern double bin_size;
// The number of bins in each row and column
extern int bins_per_side;
#define BIN_SIZE bin_size
#define BINS_PER_SIDE bins_per_side
#define NUM_BINS bins_per_side * bins_per_side

/*
Functions associated with bins
 */
/// Clear all the bins and remove any points from each
void init_bins(bin_t *bins);
void assign_particles_to_bin(particle_t &p, bin_t *bins);
void compute_forces_for_bin(bin_t *bins, int row, int col, double *dmin, double *davg, int *navg);

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
