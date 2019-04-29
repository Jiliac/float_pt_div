#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>

void testZipf(double s, int n);
double getDiv(double sp, double sq, int n);
int zipf(double alpha, int n);

int main(void) {
    //double res = getDiv(1.0, 4.0, 100);
    //printf("res: %f\n", res);

    testZipf(1.0, 100);
    testZipf(2.0, 100);
    testZipf(4.0, 100);

    return 0;
}

double getDiv(double sp, double sq, int n) {
    const int expN = 1000000;
    int x, y, i;
    double *pis, *qis;
    double pi, qi;
    double sum_p=0, sum_q=0, div=0;

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

    // Compute divergence
    printf("%f\n", div);
    for (i=0; i<n; i++) {
        pi = pis[i];
        qi = qis[i];
        sum_p += pi;
        sum_q += qi;
        div += pi * log(pi/qi);
        printf("i=%d, %f (pi, qi=%f, %f)\n", i, div, pi, qi);
    }

    div /= sum_p;
    printf("%f\n", div);
    div += log(sum_q / sum_p);

    free(pis);
    free(qis);
    return div;
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
