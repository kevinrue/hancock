#include "hancock.h"
#include "utils.h"

template<class V, typename T, class M>
Rcpp::IntegerVector num_detected_markers_internal(M ptr, Rcpp::IntegerVector markers, T threshold) {
    const size_t Ncells=ptr->get_ncol();
    
    const int Ngenes=ptr->get_nrow();
    for (auto m : markers) {
        if (m<0 || m>=Ngenes) {
            throw std::runtime_error("marker gene index out of range");
        }
    }

    Rcpp::IntegerVector output(Ncells);
    V tmp(Ngenes);

    for (size_t c=0; c<Ncells; ++c) {
        auto cIt=ptr->get_const_col(c, tmp.begin());

        int& counter=output[c];
        for (auto m : markers) {
            if (*(cIt + m) < threshold) {
                break;
            }
            ++counter;
        }
    }

    return output;
}

SEXP num_detected_markers(SEXP mat, SEXP markers, SEXP threshold) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(mat);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(mat);
        int limit=check_integer_scalar(threshold, "threshold");
        return num_detected_markers_internal<Rcpp::IntegerVector>(input.get(), markers, limit);
    } else if (rtype==REALSXP) {
        auto input=beachmat::create_numeric_matrix(mat);
        double limit=check_numeric_scalar(threshold, "threshold");
        return num_detected_markers_internal<Rcpp::NumericVector>(input.get(), markers, limit);
    } else {
        auto input=beachmat::create_logical_matrix(mat);
        return num_detected_markers_internal<Rcpp::LogicalVector, int>(input.get(), markers, 1);
    }
    END_RCPP
}
