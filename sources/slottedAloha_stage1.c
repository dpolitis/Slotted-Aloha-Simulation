/*******************************************************************************************************
* slottedAloha.c
* Computes simulated throughput, delay, mean busy stations of slotted ALOHA
* protocol, based on certain assumptions.
*
* To compile, run 'gcc slottedAloha.c -O3 -o slottedAloha -lm', on linux or
* other *nix command line.
*
* Politis Dimitrios, Mar 2017
*
* Uses parts from Mersenne Twister Random Number Generator, found at
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
*
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*
* The source code can be found at
* https://github.com/dpolitis/Slotted-Aloha-Simulation
\******************************************************************************************************/

/* header includes */
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

/* defines */
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL          // most significant 33 bits
#define LM 0x7FFFFFFFULL                  // least significant 31 bits

#define M 20                              // number of stations
#define R 100000000                       // number of time-slots per invocation

/* global variable definition */
int i, j;                                 // counters of transmissions per invocation
int B;                                    // counter of backlogged stations
unsigned int C;                           // total time-slot counter
unsigned int S;                           // total succesful time-slot counter
unsigned int busy;                        // total busy stations in each invocation
double bc[M+1][M+1];                      // binomial coefficients array

double A, F;                              // IDF functions
double X1, X2;                            // random numbers [0,1]
double p = 0.0;                           // probability of free stations
double r = 0.3;                           // probability of busy stations

/* The array for the state vector */
static unsigned long long mt[NN];

/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;

/* initializes mt[NN] with a seed */
void init_genrand64(unsigned long long seed) {
  mt[0] = seed;
  for (mti=1; mti<NN; mti++)
    mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void) {
  int i;
  unsigned long long x;
  static unsigned long long mag01[2]={0ULL, MATRIX_A};

  if (mti >= NN) {                        // generate NN words at one time

    /* if init_genrand64() has not been called, a default initial seed is used     */
    if (mti == NN+1)
      init_genrand64(5489ULL);

    for (i=0;i<NN-MM;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    for (;i<NN-1;i++) {
      x = (mt[i]&UM)|(mt[i+1]&LM);
      mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
    }
    x = (mt[NN-1]&UM)|(mt[0]&LM);
    mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

    mti = 0;
  }

  x = mt[mti++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x;
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void) {
  return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/* begin code */
int main () {
  /* fill binomial coefficients table */
  for (i = 0; i < M + 1; i++) {           // fill each line iterativelly
    for (j = 0; j < M + 1; j++) {         // using dynamic programming
      if (j > i) bc[i][j] = 0.0;
      else if ((j == 0) || (j == i)) bc[i][j] = 1.0;
      else bc[i][j] = (j <= 1 + i / 2) ? bc[i - 1][j - 1] + bc[i -1][j]: bc[i][i-j];
    }
  }

  do {

    /* initialize with random seed */
    init_genrand64((unsigned)time(NULL));

    p = p + 0.001;                        // increase probability at each invocation
    B = 0;
    C = 0;
    S = 0;
    busy = 0;

    /* until simulation is over */
    do {
      /* initialize variables */
      A = 0.0, i = 0;
      X1 = genrand64_real1();

      /* compute IDF to find how many free stations transmitted */
      do {
        A = A + bc[M - B][i] * pow(p, (double)i) * pow((1.0 - p), (double)(M - B - i));

        if (X1 < A) { break; }
          i++;
      } while (i <= (M - B));

      /* initialize variables */
      F = 0.0, j = 0;
      X2 = genrand64_real1();

      /* compute IDF to find how many busy stations transmitted */
      do {
        F = F + bc[B][j] * pow(r, (double)j) * pow((1.0 - r), (double)(B - j));

        if (X2 < F) { break; }
          j++;
      } while (j <= B);

      /* define the outcome */
      if ((i == 1) && (j == 0)) {
        /* success, free station transmitted, increase successful timeslot counter */
        S++;
      }
      else if ((i == 0) && (j == 1)) {
        /* success, busy station transmitted */
        S++;
        B--;
      }
      else {
        /* collision occured, increase busy stations by free stations transmitted,
        increase fail counter for busy stations */
        B+= i;
      }

      /* end of timeslot, proceed to next one */
      C++;
      busy+= B;

    } while (C < R);                     // end of simulation cycle

    /* printout simulated throughput, delay, mean busy stations */
    printf("%.16f \t %.16f \t %.16f %.16f\n", p, (double)S / (double)R, 1.0 + (double)(busy) / (double)S, (double)busy / (double)R);

  } while (p < 1.0);

  return 0;
}