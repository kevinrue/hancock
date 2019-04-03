#ifndef HANCOCK_H
#define HANCOCK_H

#include "Rcpp.h"

#include <stdexcept>

extern "C" {

SEXP num_detected_markers(SEXP, SEXP, SEXP);

}

#endif 
