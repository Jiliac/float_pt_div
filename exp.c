#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

void testZipf(void);
double getDiv(double sp, double sq, int np, int nq);
int zipf(double alpha, int n);

int main(void) {
    double res = getDiv(1.0, 2.0, 20, 20);
    printf("res: %f\n", res);
    return 0;
}

double getDiv(double sp, double sq, int np, int nq) {
    const int expN = 100000;
    int i;
    double pi, qi;
    double sum_p=0, sum_q=0, div=0;

    for (i=0; i<expN; i++) {
        pi = (double) zipf(sp, np);
        qi = (double) zipf(sq, nq);
        sum_p += pi;
        sum_q += qi;
        div += pi * log(pi/qi);
    }

    div /= sum_p;
    div += log(sum_q) / log(sum_p);

    return div;
}

void testZipf(void) {
    const int expN = 100000;
    int tot = 0, sqTot = 0;
    double avg, var;
    int i, z;

    for (i=0; i<expN; i++) {
        z = zipf(1, 10);
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
    static int first = 1;      // Static first time flag
    static double c = 0;          // Normalization constant
    static double *sum_probs;     // Pre-calculated sum of probabilities
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
