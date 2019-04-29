#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>

#include<gmp.h>
#include<mpfr.h>

void testZipf(double s, int n);
void getDiv(double sp, double sq, int n);
int zipf(double alpha, int n);

int main(void) {
    getDiv(1.0, 2.0, 50);
    getDiv(1.0, 2.5, 50);

    //testZipf(1.0, 100);
    //testZipf(2.0, 100);
    //testZipf(4.0, 100);

    return 0;
}

void getDiv(double sp, double sq, int n) {
    const int expN = 100000;
    const int precision = 200;
    int x, y, i;
    double *pis, *qis;
    double pi, qi;
    // Normal precision.
    double sum_p=0, sum_q=0, div=0;
    // MPFR
    mpfr_t tmp_p, pi_p, qi_p, sum_p_p, sum_q_p, div_p;

    // Init Arrays
    pis = (double*) malloc(sizeof(double) * n);
    qis = (double*) malloc(sizeof(double) * n);
    memset(pis, 0, sizeof(double) * n);
    memset(qis, 0, sizeof(double) * n);

    // Compute incidences (throughing dices).
    for (i=0; i<expN; i++) {
        x = zipf(sp, n) - 1;
        y = zipf(sq, n) - 1;
        pis[x] = pis[x] + 1.0;
        qis[y] = qis[y] + 1.0;
    }

    // Compute divergence (Normal precision)
    for (i=0; i<n; i++) {
        pi = pis[i];
        qi = qis[i];
        sum_p += pi;
        sum_q += qi;
        div += pi * log(pi/qi);
    }
    //
    div /= sum_p;
    div += log(sum_q / sum_p);
    printf("(Normal)\tsp,sq=%f,%f, n=%d\t-- div=%f\n", sp, sq, n, div);

    // Compute divergence (high precision)
    // a. Init
    mpfr_init2(tmp_p, precision);
    mpfr_init2(pi_p, precision);
    mpfr_init2(qi_p, precision);
    mpfr_init2(sum_p_p, precision);
    mpfr_init2(sum_q_p, precision);
    mpfr_init2(div_p, precision);
    //
    mpfr_set_d(sum_p_p, 0.0, MPFR_RNDD);
    mpfr_set_d(sum_q_p, 0.0, MPFR_RNDD);
    mpfr_set_d(div_p, 0.0, MPFR_RNDD);
    // b. Go over incidences.
    for (i=0; i<expN; i++) {
        mpfr_set_d(pi_p, pis[i], MPFR_RNDD);
        mpfr_set_d(qi_p, qis[i], MPFR_RNDD);
    }
    // c. Post loop
    // d. Clearing
    mpfr_clear(tmp_p);
    mpfr_clear(pi_p);
    mpfr_clear(qi_p);
    mpfr_clear(sum_p_p);
    mpfr_clear(sum_q_p);
    mpfr_clear(div_p);
    mpfr_free_cache();

    free(pis);
    free(qis);
}

void testZipf(double s, int n) {
    const int expN = 100000;
    int tot = 0, sqTot = 0;
    double avg, var;
    int i, z;

    for (i=0; i<expN; i++) {
        z = zipf(s, n);
        tot += z;
        sqTot += z * z;
    }
    avg = ((double) tot) / expN;
    var = ((double) sqTot) / expN;
    var -= avg*avg;
    var = sqrt(var);
    printf("%f - %f\n", avg, var);
}

int zipf(double alpha, int n) {
    int first = 1;      // Static first time flag
    double c = 0;          // Normalization constant
    double *sum_probs;     // Pre-calculated sum of probabilities
    double z;                     // Uniform random number (0 < z < 1)
    int zipf_value;               // Computed exponential value to be returned
    int    i;                     // Loop counter
    int low, high, mid;           // Binary-search bounds

    // Compute normalization constant on first call only
    if (first) {
        for (i=1; i<=n; i++)
            c = c + (1.0 / pow((double) i, alpha));
        c = 1.0 / c;

        sum_probs = malloc((n+1)*sizeof(*sum_probs));
        sum_probs[0] = 0;
        for (i=1; i<=n; i++) {
            sum_probs[i] = sum_probs[i-1] + c / pow((double) i, alpha);
        }
        first = 0;
    }

    // Pull a uniform random number (0 < z < 1)
    do {
        z = (double) (rand() / (RAND_MAX + 1.));
    }
    while ((z == 0) || (z == 1));

    // Map z to the value
    low = 1, high = n, mid;
    do {
        mid = floor((low+high)/2);
        if (sum_probs[mid] >= z && sum_probs[mid-1] < z) {
            zipf_value = mid;
            break;
        } else if (sum_probs[mid] >= z) {
            high = mid-1;
        } else {
            low = mid+1;
        }
    } while (low <= high);

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >=1) && (zipf_value <= n));

    return zipf_value;
}
