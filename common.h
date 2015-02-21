#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <list>

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

//
//  timing routines
//
double read_timer( );

typedef struct 
{
    double x;
    double y;
    double x_length;
    double y_length;
    std::list<particle_t> *particles;
    inline int numParticles() {return particles -> size();}
    // Checks whether the particle is out of this bin
    inline bool outOfBound(particle_t particle) {return (particle.x < x || particle.x > x + x_length || particle.y < y || particle.y > y + y_length);}
    // Add a particle to the bin
    inline void addParticle(particle_t particle) {particles -> push_back(particle);}

} bin_t;

// Initialize all the bins
bin_t* initBins();
int getBinNum();

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
