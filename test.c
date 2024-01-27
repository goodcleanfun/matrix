#include <stdio.h>
#include <stdlib.h>

#include "greatest/greatest.h"
#include "double_matrix.h"


TEST test_double_matrix(void) {
    size_t n = 10;
    double_matrix *M = double_matrix_new(n, n);
    for (size_t i = 0; i < M->m; i++) {
        for (size_t j = 0; j < M->n; j++) {
            double_matrix_set_at_index(M, i, j, (double)i + 1);
        }
    }
    for (size_t i = 0; i < M->m; i++) {
        for (size_t j = 0; j < M->n; j++) {
            ASSERT_EQ(double_matrix_get(M, i, j), (double)i + 1);
        }
    }

    double_matrix *N = double_matrix_new(n, n);
    for (size_t i = 0; i < N->m; i++) {
        for (size_t j = 0; j < N->n; j++) {
            double_matrix_set_at_index(N, i, j, (double)j + 1);
        }
    }

    double_matrix *P = double_matrix_new(n, n);
    ASSERT(double_matrix_dot_matrix(M, N, P));

    for (size_t i = 0; i < P->m; i++) {
        for (size_t j = 0; j < P->n; j++) {
            ASSERT_EQ(double_matrix_get(P, i, j), (double)(i + 1) * (j + 1) * n);
        }
    }

    double_matrix_zero(P);
    for (size_t i = 0; i < P->m; i++) {
        for (size_t j = 0; j < P->n; j++) {
            ASSERT_EQ(double_matrix_get(P, i, j), 0.0);
        }
    }

    // Not commutative
    ASSERT(double_matrix_dot_matrix(N, M, P));
    double expected_dot = 0.0;
    for (size_t i = 0; i < P->m; i++) {
        expected_dot += (double)((i + 1) * (i + 1));
    }
    for (size_t i = 0; i < P->m; i++) {
        for (size_t j = 0; j < P->n; j++) {
            ASSERT_EQ(double_matrix_get(P, i, j), expected_dot);
        }
    }

    double_matrix_destroy(M);
    double_matrix_destroy(N);
    double_matrix_destroy(P);
    return 0;
}

/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();      /* command-line options, initialization. */

    RUN_TEST(test_double_matrix);

    GREATEST_MAIN_END();        /* display results */
}