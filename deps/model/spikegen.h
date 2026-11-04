#ifndef SPIKEGEN_H
#define SPIKEGEN_H

/* Parameters of the refractory function in the Zilany, Bruce, and Carney (2014) spike
   generator: amplitudes (c0, c1) and time constants (s0, s1, in s) of the two exponentials
   of relative refractoriness, and the absolute refractory period (dead, in s) */
#define SPIKEGEN_ZBC2014_C0 0.5
#define SPIKEGEN_ZBC2014_S0 0.001
#define SPIKEGEN_ZBC2014_C1 0.5
#define SPIKEGEN_ZBC2014_S1 0.0125
#define SPIKEGEN_ZBC2014_DEAD 0.00075

/* Random numbers used by spikegen_zbc2014 beyond one per spike */
#define SPIKEGEN_ZBC2014_EXTRA_DRAWS 2

int spikegen_zbc2014_max_spikes(int n, double tdres);
int spikegen_zbc2014(const double *rate, int n, double tdres, const double *rand_unif,
                     int n_rand, double *sptime, int n_sptime);

#endif
