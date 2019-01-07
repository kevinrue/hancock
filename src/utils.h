#ifndef UTILS_H
#define UTILS_H

// Checking for scalar inputs.

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);


#endif
