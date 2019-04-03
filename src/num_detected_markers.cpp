#include "hancock.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/logical_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

template<class M>
Rcpp::IntegerVector num_detected_markers_internal(Rcpp::RObject mat, Rcpp::IntegerVector markers, typename M::type threshold) {
    auto ptr=beachmat::create_matrix<M>(mat);
    const size_t Ncells=ptr->get_ncol();
    const int Ngenes=ptr->get_nrow();
    for (auto m : markers) {
        if (m<0 || m>=Ngenes) {
            throw std::runtime_error("marker gene index out of range");
        }
    }

    Rcpp::IntegerVector output(Ncells);
    beachmat::const_column<M> col_holder(ptr.get(), false); // need indexing, so no sparsity allowed.

    for (size_t c=0; c<Ncells; ++c) {
        col_holder.fill(c);
        auto cIt=col_holder.get_values();

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
        int limit=check_integer_scalar(threshold, "threshold");
        return num_detected_markers_internal<beachmat::integer_matrix>(mat, markers, limit);
    } else if (rtype==REALSXP) {
        double limit=check_numeric_scalar(threshold, "threshold");
        return num_detected_markers_internal<beachmat::numeric_matrix>(mat, markers, limit);
    } else {
        return num_detected_markers_internal<beachmat::logical_matrix>(mat, markers, 1);
    }
    END_RCPP
}
