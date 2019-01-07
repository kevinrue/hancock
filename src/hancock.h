#ifndef HANCOCK_H
#define HANCOCK_H

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/logical_matrix.h"
#include "Rcpp.h"

#include <stdexcept>

extern "C" {

SEXP num_detected_markers(SEXP, SEXP, SEXP);

}

#endif 
