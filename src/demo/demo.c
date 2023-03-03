/*
 *  This is a set of testbench code for handling large dynamically allocated arrays in C
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void accept_vector(double *x) {
    /* Print simple confirmation that everything is working */
    printf("Yo");
}   

void accept_matrix(double **x) {
    /* Print simple confirmation that everything is working */
    printf("Yo");
    printf("\n");

    /* Loop through elements and print */
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            printf("%.6f", x[i][j]);
            printf("\n");
        }
    }
}   
