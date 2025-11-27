#ifndef FLOAT_MATRIX_H
#define FLOAT_MATRIX_H

#define MATRIX_NAME float_matrix
#define MATRIX_TYPE float
#define MATRIX_BLAS_PRECISION s
#include "vectorized/float_vector.h"
#include "matrix.h"

static inline void float_matrix_log(float_matrix *self) {
    float_vector_log(self->values, self->m * self->n);
}

static inline void float_matrix_exp(float_matrix *self) {
    float_vector_exp(self->values, self->m * self->n);
}

#undef MATRIX_NAME
#undef MATRIX_TYPE
#undef MATRIX_BLAS_PRECISION

#endif
