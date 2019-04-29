#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>

#include<gmp.h>
#include<mpfr.h>
#include <cstdlib>
#include <algorithm>
#include <vector>

void testZipf(double s, int n);
void getDiv(double sp, double sq, int n);
int zipf(double alpha, int n);


double kl_normal(int n, double s_p, double s_q, 
          int samples_p, int samples_q,
          std::vector<int> indices_p, std::vector<int> indices_q) {
    int x, y, i;
    double *pis, *qis;
    double pi, qi;
    // Normal precision.
    double sum_p=0, sum_q=0, div=0;

    // Init Arrays
    pis = (double*) malloc(sizeof(double) * n);
    qis = (double*) malloc(sizeof(double) * n);
    memset(pis, 0, sizeof(double) * n);
    memset(qis, 0, sizeof(double) * n);


    // Compute incidences (throughing dices).
    for (i = 0; i < samples_p; i++) {
        x = indices_p[zipf(s_p, n) - 1];
        pis[x] = pis[x] + 1.0;
    }
    for (i = 0; i < samples_q; i++) {
        x = indices_q[zipf(s_q, n) - 1];
        qis[x] = qis[x] + 1.0;
    }

    for (i=0; i<n; i++) {
        pi = pis[i];
        qi = qis[i];
        sum_p += pi;
        sum_q += qi;
        if (pi < .5 || qi < .5) {
            continue;
        }
        div += pi * log(pi/qi);
    }
    //
    div /= sum_p;
    div += log(sum_q / sum_p);

    free(pis);
    free(qis);

    return div;
}

double kl_precise(int n, double s_p, double s_q, 
          int samples_p, int samples_q,
          std::vector<int> indices_p, std::vector<int> indices_q) {
    int x, y, i;
    double *pis, *qis;
    double pi, qi;
    const int precision = 200;

    // MPFR
    mpfr_t tmp_p, pi_p, qi_p, sum_p_p, sum_q_p, div_p;

    // Init Arrays
    pis = (double*) malloc(sizeof(double) * n);
    qis = (double*) malloc(sizeof(double) * n);
    memset(pis, 0, sizeof(double) * n);
    memset(qis, 0, sizeof(double) * n);


    // Compute incidences (throughing dices).
    for (i = 0; i < samples_p; i++) {
        x = indices_p[zipf(s_p, n) - 1];
        pis[x] = pis[x] + 1.0;
    }
    for (i = 0; i < samples_q; i++) {
        x = indices_q[zipf(s_q, n) - 1];
        qis[x] = qis[x] + 1.0;
    }

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
    for (i=0; i<n; i++) {
        mpfr_set_d(pi_p, pis[i], MPFR_RNDD);
        mpfr_set_d(qi_p, qis[i], MPFR_RNDD);
        mpfr_add(sum_p_p, sum_p_p, pi_p, MPFR_RNDD);
        mpfr_add(sum_q_p, sum_q_p, qi_p, MPFR_RNDD);
        //
        if (pis[i] < .5 || qis[i] < .5) {
            continue;
        }
        mpfr_div(tmp_p, pi_p, qi_p, MPFR_RNDD);     // tmp = pi/qi
        mpfr_log(tmp_p, tmp_p, MPFR_RNDD);          // tmp = log(pi/qi)
        mpfr_mul(tmp_p, pi_p, tmp_p, MPFR_RNDD);    // tmp = pi*log(pi/qi)
        mpfr_add(div_p, div_p, tmp_p, MPFR_RNDD);
    }
    // c. Post loop
    mpfr_div(div_p, div_p, sum_p_p, MPFR_RNDD);     // div /= sum_p
    mpfr_div(tmp_p, sum_q_p, sum_p_p, MPFR_RNDD);   // tmp = sum_q / sum_p
    mpfr_log(tmp_p, tmp_p, MPFR_RNDD);              // tmp = log(sum_q / sum_p)
    mpfr_add(div_p, div_p, tmp_p, MPFR_RNDD);       // div += tmp
    //

    double div = mpfr_get_d(div_p, MPFR_RNDD);

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

    return div;
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
        if (pi < .5) {
            continue;
        }
        div += pi * log(pi/qi);
    }
    //
    div /= sum_p;
    div += log(sum_q / sum_p);
    printf("sp,sq=%f,%f, n=%d\n(Normal)\tdiv=%f\n", sp, sq, n, div);

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
    for (i=0; i<n; i++) {
        mpfr_set_d(pi_p, pis[i], MPFR_RNDD);
        mpfr_set_d(qi_p, qis[i], MPFR_RNDD);
        mpfr_add(sum_p_p, sum_p_p, pi_p, MPFR_RNDD);
        mpfr_add(sum_q_p, sum_q_p, qi_p, MPFR_RNDD);
        //
        if (pis[i] < .5) {
            continue;
        }
        mpfr_div(tmp_p, pi_p, qi_p, MPFR_RNDD);     // tmp = pi/qi
        mpfr_log(tmp_p, tmp_p, MPFR_RNDD);          // tmp = log(pi/qi)
        mpfr_mul(tmp_p, pi_p, tmp_p, MPFR_RNDD);    // tmp = pi*log(pi/qi)
        mpfr_add(div_p, div_p, tmp_p, MPFR_RNDD);
    }
    // c. Post loop
    mpfr_div(div_p, div_p, sum_p_p, MPFR_RNDD);     // div /= sum_p
    mpfr_div(tmp_p, sum_q_p, sum_p_p, MPFR_RNDD);   // tmp = sum_q / sum_p
    mpfr_log(tmp_p, tmp_p, MPFR_RNDD);              // tmp = log(sum_q / sum_p)
    mpfr_add(div_p, div_p, tmp_p, MPFR_RNDD);       // div += tmp
    //
    printf("(Precision=%d)\tdiv=", precision);
    mpfr_out_str(stdout, 10, 10, div_p, MPFR_RNDD);
    printf("\n");
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
    double c = 0;          // Normalization constant
    double *sum_probs;     // Pre-calculated sum of probabilities
    double z;                     // Uniform random number (0 < z < 1)
    int zipf_value;               // Computed exponential value to be returned
    int    i;                     // Loop counter
    int low, high, mid;           // Binary-search bounds

    // Compute normalization constant on first call only
    for (i=1; i<=n; i++)
        c = c + (1.0 / pow((double) i, alpha));
    c = 1.0 / c;

    sum_probs = (double*)malloc((n+1)*sizeof(*sum_probs));
    sum_probs[0] = 0;
    for (i=1; i<=n; i++) {
        sum_probs[i] = sum_probs[i-1] + c / pow((double) i, alpha);
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

    free(sum_probs);

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >=1) && (zipf_value <= n));

    return zipf_value;
}


int main(void) {

  int species = 50;
  double div_n;
  double div_p;
  int i;

  // We want to avoid that index i is *always* greater than index i+1
  std::vector<int> indices_p;
  std::vector<int> indices_q;
  for (int i = 0; i < species; ++i) {
     indices_p.push_back(i);
     indices_q.push_back(i);
  }
  std::random_shuffle(indices_p.begin(), indices_p.end());
  std::random_shuffle(indices_q.begin(), indices_q.end());


  printf("Normal, Precise, Difference\n");
  for (i = 0; i < 50; i++) {
    div_n = kl_normal(species, 2.0, 3.0, 100000, 100000, indices_p, indices_q);
    div_p = kl_precise(species, 2.0, 3.0, 100000, 100000, indices_p, indices_q);

    printf("%20f, %20f, %20f\n", div_n, div_p, div_p - div_n);
  } 

  //  getDiv(1.5, 0.5, 50);
    //testZipf(1.0, 100);
    //testZipf(2.0, 100);
    //testZipf(4.0, 100);

    return 0;
}
