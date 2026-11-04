#include <math.h>

#include "spikegen.h"

/**
 * spikegen_zbc2014_max_spikes
 *
 * Returns the maximum number of spikes that spikegen_zbc2014 can generate from a rate
 * function of n samples at time resolution tdres, for sizing the buffers passed to it
 */
int spikegen_zbc2014_max_spikes(int n, double tdres) {
    int deadtime_index = (int) floor(SPIKEGEN_ZBC2014_DEAD / tdres);
    return (n + deadtime_index) / (deadtime_index + 1);
}

/**
 * spikegen_zbc2014
 *
 * Passes the output of the synapse model through the spike generator
 *
 * The spike generator uses a method coded up by B. Scott Jackson (bsj22@cornell.edu).
 *
 * This is a revision of the SpikeGenerator function in the 2014 model code (model_Synapse.c).
 * The differences from that code are:
 * - The 2014 code computed the end of the dead time of the spike before t = 0 with the
 *   macro __max(0, log(randNums[randBufIndex++]) / synouttmp[0] + dead), which evaluates
 *   its second argument twice. When that argument was positive, the macro consumed a second
 *   random number, returned a value computed from it without clamping it to 0, and shifted
 *   all subsequent random numbers by one. Here, one random number is consumed and the
 *   result is clamped, as intended.
 * - Random numbers are read from a buffer filled by the caller, instead of being requested
 *   through a callback. The buffer must contain at least
 *   spikegen_zbc2014_max_spikes(n, tdres) + SPIKEGEN_ZBC2014_EXTRA_DRAWS values drawn
 *   uniformly from (0, 1).
 * - The caller sizes and owns all buffers, instead of the function allocating them
 *   internally. The spike-time buffer must hold at least spikegen_zbc2014_max_spikes(n,
 *   tdres) values. The function checks both buffer sizes and returns -1 instead of reading
 *   or writing past them.
 *
 * @param rate Synapse output (sp/s), of length n
 * @param n Number of samples in rate
 * @param tdres Time resolution of rate (s)
 * @param rand_unif Uniform random numbers, of length n_rand
 * @param n_rand Number of values in rand_unif
 * @param sptime Output buffer for spike times (s), of length n_sptime
 * @param n_sptime Number of values that sptime can hold
 * @return Number of spikes written to sptime, or -1 if rand_unif or sptime is too small
 */
int spikegen_zbc2014(const double *rate, int n, double tdres, const double *rand_unif,
                     int n_rand, double *sptime, int n_sptime) {
    const double c0 = SPIKEGEN_ZBC2014_C0;
    const double s0 = SPIKEGEN_ZBC2014_S0;
    const double c1 = SPIKEGEN_ZBC2014_C1;
    const double s1 = SPIKEGEN_ZBC2014_S1;
    const double dead = SPIKEGEN_ZBC2014_DEAD;

    if (n <= 0) return 0;
    if (n_rand < SPIKEGEN_ZBC2014_EXTRA_DRAWS) return -1;

    double dur = n * tdres;  /* Total duration of the rate function */
    int n_spikes = 0;
    int idx_rand = 0;

    /* Calculate useful constants */
    int deadtime_index = (int) floor(dead / tdres);  /* Integer number of discrete time bins within deadtime */
    double deadtime_rnd = deadtime_index * tdres;    /* Deadtime rounded down to length of an integer number of discrete time bins */

    double refrac_mult0 = 1.0 - tdres / s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+tdres) = y0(t)*refrac_mult0 */
    double refrac_mult1 = 1.0 - tdres / s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+tdres) = y1(t)*refrac_mult1 */

    /* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
    double end_last_deadtime = fmax(0.0, log(rand_unif[idx_rand++]) / rate[0] + dead);  /* End of last deadtime before t=0 */
    double refrac_value0 = c0 * exp(end_last_deadtime / s0);  /* Value of first exponential in refractory function */
    double refrac_value1 = c1 * exp(end_last_deadtime / s1);  /* Value of second exponential in refractory function */
    double xsum = rate[0] * (-end_last_deadtime + c0 * s0 * (exp(end_last_deadtime / s0) - 1.0)
                                                + c1 * s1 * (exp(end_last_deadtime / s1) - 1.0));
        /* Value of time-warping sum */
        /*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */

    /* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
    double unit_rate_interval = -log(rand_unif[idx_rand++]) / tdres;
        /* NOTE: Both 'unit_rate_interval' and 'xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.
           This way we only need to divide by 'tdres' once per spike (when calculating 'unit_rate_interval'), instead of
           multiplying by 'tdres' once per time bin (when calculating the new value of 'xsum').                         */

    double count_time = tdres;
    for (int k = 0; (k < n) && (count_time < dur);
         ++k, count_time += tdres, refrac_value0 *= refrac_mult0, refrac_value1 *= refrac_mult1) {  /* Loop through rate vector */
        if (rate[k] > 0.0) {  /* Nothing to do for non-positive rates, i.e. xsum += 0 for non-positive rates. */
            xsum += rate[k] * (1.0 - refrac_value0 - refrac_value1);  /* Add synout*(refractory value) to time-warping sum */

            if (xsum >= unit_rate_interval) {  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
                if (n_spikes >= n_sptime || idx_rand >= n_rand) return -1;
                sptime[n_spikes++] = count_time;
                unit_rate_interval = -log(rand_unif[idx_rand++]) / tdres;
                xsum = 0.0;

                /* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
                k += deadtime_index;
                count_time += deadtime_rnd;
                refrac_value0 = c0;
                refrac_value1 = c1;
            }
        }
    } /* End of rate vector loop */

    return n_spikes;  /* Number of spikes that occurred. */
}
