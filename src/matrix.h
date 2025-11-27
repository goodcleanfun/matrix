#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include "file_utils/file_utils.h"

#endif

#ifdef HAVE_CBLAS
#include <cblas.h>
#else
#warning "No CBLAS"
#endif


#ifndef MATRIX_NAME
#error "MATRIX_NAME must be defined"
#endif

#ifndef MATRIX_TYPE
#error "MATRIX_TYPE must be defined"
#endif

// Prefer aligned memory for matrices
#if !defined(MATRIX_MALLOC) && !defined(MATRIX_REALLOC) && !defined(MATRIX_FREE)
#include "aligned/aligned.h"

#define MATRIX_ALIGNED
#define MATRIX_ALIGNED_DEFINED

#ifndef MATRIX_MALLOC
#if defined(MATRIX_ALIGNMENT)
#define MATRIX_MALLOC(size) aligned_malloc(size, MATRIX_ALIGNMENT)
#else
#define MATRIX_MALLOC default_aligned_malloc
#define MATRIX_MALLOC_DEFINED
#endif
#endif

#ifndef MATRIX_REALLOC
#define MATRIX_REALLOC_NEEDS_PREV_SIZE

#if defined(MATRIX_ALIGNMENT)
#define MATRIX_REALLOC(a, prev_size, new_size) aligned_resize(a, prev_size, new_size, MATRIX_ALIGNMENT)
#else
#define MATRIX_REALLOC default_aligned_resize
#define MATRIX_REALLOC_DEFINED
#endif
#endif

#ifndef MATRIX_FREE
#define MATRIX_FREE default_aligned_free
#define MATRIX_FREE_DEFINED
#endif

#endif

#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define CONCAT3_(a, b, c) a ## b ## c
#define CONCAT3(a, b, c) CONCAT3_(a, b, c)
#define MATRIX_FUNC(name) CONCAT(MATRIX_NAME, _##name)
#define VECTOR_FUNC(name) CONCAT3(MATRIX_TYPE, _vector, _##name)
#define FILE_ARRAY_FUNC(name) CONCAT3(name##_, MATRIX_TYPE, _array)

typedef struct {
    size_t m, n;
    MATRIX_TYPE *values;
} MATRIX_NAME;


static MATRIX_NAME *MATRIX_FUNC(new)(size_t m, size_t n) {
    MATRIX_NAME *matrix = malloc(sizeof(MATRIX_NAME));

    if (matrix == NULL) {
        return NULL;
    }

    matrix->m = m;
    matrix->n = n;

    matrix->values = MATRIX_MALLOC(sizeof(MATRIX_TYPE) * m * n);
    if (matrix->values == NULL) {
        free(matrix);
        return NULL;
    }

    return matrix;

}

static void MATRIX_FUNC(destroy)(MATRIX_NAME *self) {
    if (self == NULL) return;

    if (self->values != NULL) {
        MATRIX_FREE(self->values);
    }

    free(self);
}

static inline void MATRIX_FUNC(zero)(MATRIX_NAME *self) {
    memset(self->values, 0, self->m * self->n * sizeof(MATRIX_TYPE));
}


static inline bool MATRIX_FUNC(resize)(MATRIX_NAME *self, size_t m, size_t n) {
    if (self == NULL) return false;

    if (m * n > (self->m * self->n)) {
        #ifdef MATRIX_REALLOC_NEEDS_PREV_SIZE
        MATRIX_TYPE *ptr = MATRIX_REALLOC(self->values, sizeof(MATRIX_TYPE) * self->m * self->n, sizeof(MATRIX_TYPE) * m * n);
        #else
        MATRIX_TYPE *ptr = MATRIX_REALLOC(self->values, sizeof(MATRIX_TYPE) * m * n);
        #endif
        if (ptr == NULL) {
            return false;
        }
        self->values = ptr;
    }

    self->m = m;
    self->n = n;

    return true;
}

static inline bool MATRIX_FUNC(resize_fill_zeros)(MATRIX_NAME *self, size_t m, size_t n) {
    size_t old_m = self->m;
    bool ret = MATRIX_FUNC(resize)(self, m, n);
    if (ret && m > old_m) {
        memset(self->values + old_m, 0, (m - old_m) * self->n * sizeof(MATRIX_TYPE));
    }
    return ret;
}

static inline MATRIX_NAME *MATRIX_FUNC(new_copy)(MATRIX_NAME *self) {
    MATRIX_NAME *cpy = MATRIX_FUNC(new)(self->m, self->n);
    size_t num_values = self->m * self->n;
    memcpy(cpy->values, self->values, num_values * sizeof(MATRIX_TYPE));

    return cpy;
}

static inline bool MATRIX_FUNC(copy)(MATRIX_NAME *self, MATRIX_NAME *other) {
    if (self->m != other->m || self->n != other->n) {
        return false;
    }
    size_t num_values = self->m * self->n;

    memcpy(other->values, self->values, num_values * sizeof(MATRIX_TYPE));
    return true;
}

static inline void MATRIX_FUNC(init_values)(MATRIX_NAME *self, MATRIX_TYPE *values) {
    size_t num_values = self->m * self->n;
    memcpy(self->values, values, num_values * sizeof(MATRIX_TYPE));
}

static inline void MATRIX_FUNC(set)(MATRIX_NAME *self, MATRIX_TYPE value) {
    VECTOR_FUNC(set)(self->values, value, self->m * self->n);
}

static inline void MATRIX_FUNC(set_row)(MATRIX_NAME *self, size_t index, MATRIX_TYPE *row) {
    size_t offset = index * self->n;
    MATRIX_TYPE *values = self->values;
    size_t n = self->n;
    memcpy(values + offset, row, n * sizeof(MATRIX_TYPE));
}

static inline void MATRIX_FUNC(set_at_index)(MATRIX_NAME *self, size_t row_index, size_t col_index, MATRIX_TYPE value) {
    size_t offset = row_index * self->n + col_index;
    self->values[offset] = value;
}

static inline void MATRIX_FUNC(add_at_index)(MATRIX_NAME *self, size_t row_index, size_t col_index, MATRIX_TYPE value) {
    size_t offset = row_index * self->n + col_index;
    self->values[offset] += value;
}

static inline void MATRIX_FUNC(sub_at_index)(MATRIX_NAME *self, size_t row_index, size_t col_index, MATRIX_TYPE value) {
    size_t offset = row_index * self->n + col_index;
    self->values[offset] -= value;
}

static inline void MATRIX_FUNC(mul_at_index)(MATRIX_NAME *self, size_t row_index, size_t col_index, MATRIX_TYPE value) {
    size_t offset = row_index * self->n + col_index;
    self->values[offset] *= value;
}

static inline void MATRIX_FUNC(div_at_index)(MATRIX_NAME *self, size_t row_index, size_t col_index, MATRIX_TYPE value) {
    size_t offset = row_index * self->n + col_index;
    self->values[offset] /= value;
}

static inline MATRIX_TYPE MATRIX_FUNC(get)(MATRIX_NAME *self, size_t row_index, size_t col_index) {
    size_t index = row_index * self->n + col_index;
    return self->values[index];
}

static inline MATRIX_TYPE *MATRIX_FUNC(get_row)(MATRIX_NAME *self, size_t row_index) {
    size_t index = row_index * self->n;
    return self->values + index;
}

static inline MATRIX_NAME *MATRIX_FUNC(new_value)(size_t m, size_t n, MATRIX_TYPE value) {
    MATRIX_NAME *matrix = MATRIX_FUNC(new)(m, n);
    MATRIX_FUNC(set)(matrix, value);
    return matrix;
}

static inline MATRIX_NAME *MATRIX_FUNC(new_zeros)(size_t m, size_t n) {
    MATRIX_NAME *matrix = MATRIX_FUNC(new)(m, n);
    MATRIX_FUNC(zero)(matrix);
    return matrix;
}

static inline MATRIX_NAME *MATRIX_FUNC(new_ones)(size_t m, size_t n) {
    return MATRIX_FUNC(new_value)(m, n, (MATRIX_TYPE)1);
}

static inline MATRIX_NAME *MATRIX_FUNC(new_values)(size_t m, size_t n, MATRIX_TYPE *values) {
    MATRIX_NAME *matrix = MATRIX_FUNC(new)(m, n);
    memcpy(matrix->values, values, m * n * sizeof(MATRIX_TYPE));
    return matrix;
}

static inline void MATRIX_FUNC(div)(MATRIX_NAME *self, MATRIX_TYPE value) {
    VECTOR_FUNC(div)(self->values, value, self->m * self->n);
}

static inline bool MATRIX_FUNC(div_matrix)(MATRIX_NAME *self, MATRIX_NAME *other) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(div_vector)(self->values, other->values, self->m * self->n);
    return true;
}

static inline bool MATRIX_FUNC(div_matrix_scaled)(MATRIX_NAME *self, MATRIX_NAME *other, MATRIX_TYPE v) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(div_vector_scaled)(self->values, other->values, v, self->m * self->n);
    return true;
}

static inline void MATRIX_FUNC(mul)(MATRIX_NAME *self, MATRIX_TYPE value) {
    VECTOR_FUNC(mul)(self->values, value, self->m * self->n);
}

static inline bool MATRIX_FUNC(mul_matrix)(MATRIX_NAME *self, MATRIX_NAME *other) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(mul_vector)(self->values, other->values, self->m * self->n);
    return true;
}

static inline bool MATRIX_FUNC(mul_matrix_scaled)(MATRIX_NAME *self, MATRIX_NAME *other, MATRIX_TYPE v) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(mul_vector_scaled)(self->values, other->values, v, self->m * self->n);
    return true;
}

static inline void MATRIX_FUNC(add)(MATRIX_NAME *self, MATRIX_TYPE value) {
    VECTOR_FUNC(add)(self->values, self->m * self->n, value);
}


static inline bool MATRIX_FUNC(add_matrix)(MATRIX_NAME *self, MATRIX_NAME *other) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(add_vector)(self->values, other->values, self->m * self->n);
    return true;
}

static inline bool MATRIX_FUNC(add_matrix_scaled)(MATRIX_NAME *self, MATRIX_NAME *other, MATRIX_TYPE v) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(add_vector_scaled)(self->values, other->values, v, self->m * self->n);
    return true;
}

static inline void MATRIX_FUNC(sub)(MATRIX_NAME *self, MATRIX_TYPE value) {
    VECTOR_FUNC(sub)(self->values, value, self->m * self->n);
}

static inline bool MATRIX_FUNC(sub_matrix)(MATRIX_NAME *self, MATRIX_NAME *other) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(sub_vector)(self->values, other->values, self->m * self->n);
    return true;
}

static inline bool MATRIX_FUNC(sub_matrix_scaled)(MATRIX_NAME *self, MATRIX_NAME *other, MATRIX_TYPE v) {
    if (self->m != other->m || self->n != other->n) return false;
    VECTOR_FUNC(sub_vector_scaled)(self->values, other->values, v, self->m * self->n);
    return true;
}

static inline void MATRIX_FUNC(dot_vector)(MATRIX_NAME *self, MATRIX_TYPE *vec, MATRIX_TYPE *result) {
    MATRIX_TYPE *values = self->values;
    size_t m = self->m;
    size_t n = self->n;
    #pragma omp parallel for reduction (+:result[:m])
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i] += values[n * i + j] * vec[j];
        }
    }
}

#if defined(HAVE_CBLAS) && defined(MATRIX_BLAS_PRECISION)

#define CBLAS_FUNC(name) CONCAT3(cblas_, MATRIX_BLAS_PRECISION, name)
#define GEMM CBLAS_FUNC(gemm)
#undef CBLAS_FUNC

static inline bool MATRIX_FUNC(dot_matrix)(MATRIX_NAME *m1, MATRIX_NAME *m2, MATRIX_NAME *result) {
    if (m1->n != m2->m || m1->m != result->m || m2->n != result->n) {
        return false;
    }

    GEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans,
         m1->m, m2->n, m1->n, 1.0,
         m1->values, m1->n,
         m2->values, m2->n, 0.0,
         result->values, result->n
         );

    return true;
}

#undef GEMM

#else
static inline bool MATRIX_FUNC(dot_matrix)(MATRIX_NAME *m1, MATRIX_NAME *m2, MATRIX_NAME *result) {
    if (m1->n != m2->m || m1->m != result->m || m2->n != result->n) {
        return false;
    }

    size_t m1_rows = m1->m;
    size_t m1_cols = m1->n;
    size_t m2_rows = m2->m;
    size_t m2_cols = m2->n;

    MATRIX_TYPE *m1_values = m1->values;
    MATRIX_TYPE *m2_values = m2->values;
    MATRIX_TYPE *result_values = result->values;

    size_t i, j, k;
    #pragma omp parallel for private(i, j, k) shared(m1_values, m2_values, result_values)
    for (i = 0; i < m1_rows; i++) {
        for (k = 0; k < m2_rows; k++) {
            MATRIX_TYPE r = m1_values[m1_cols * i + k];
            for (j = 0; j < m2_cols; j++) {
                result_values[m2_cols * i + j] += r * m2_values[m2_cols * k + j];
            }
        }
    }

    return true;
}
#endif

#undef CONCAT_
#undef CONCAT
#undef CONCAT3_
#undef CONCAT3
#undef MATRIX_FUNC
#undef VECTOR_FUNC
#undef FILE_ARRAY_FUNC
#ifdef MATRIX_ALIGNMENT_DEFINED
#undef MATRIX_ALIGNMENT
#undef MATRIX_ALIGNMENT_DEFINED
#endif

#ifdef MATRIX_ALIGNED_DEFINED
#undef MATRIX_ALIGNED
#undef MATRIX_ALIGNED_DEFINED
#endif

#ifdef MATRIX_MALLOC_DEFINED
#undef MATRIX_MALLOC
#undef MATRIX_MALLOC_DEFINED
#endif
#ifdef MATRIX_REALLOC_DEFINED
#undef MATRIX_REALLOC
#undef MATRIX_REALLOC_DEFINED
#endif
#ifdef MATRIX_FREE_DEFINED
#undef MATRIX_FREE
#undef MATRIX_FREE_DEFINED
#endif