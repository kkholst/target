// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/target.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// mlogit_loglik
double mlogit_loglik(arma::vec theta, const arma::uvec& choice, const arma::uvec& alt, unsigned basealt, unsigned nalt, const arma::uvec& id_idx, const arma::mat& z1, const arma::mat& z2, const arma::mat& x, const arma::vec& weights);
RcppExport SEXP _target_mlogit_loglik(SEXP thetaSEXP, SEXP choiceSEXP, SEXP altSEXP, SEXP basealtSEXP, SEXP naltSEXP, SEXP id_idxSEXP, SEXP z1SEXP, SEXP z2SEXP, SEXP xSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type alt(altSEXP);
    Rcpp::traits::input_parameter< unsigned >::type basealt(basealtSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nalt(naltSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_idx(id_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(mlogit_loglik(theta, choice, alt, basealt, nalt, id_idx, z1, z2, x, weights));
    return rcpp_result_gen;
END_RCPP
}
// mlogit_pred
arma::mat mlogit_pred(arma::vec theta, const arma::uvec& alt, unsigned basealt, unsigned nalt, const arma::uvec& id_idx, const arma::mat& z1, const arma::mat& z2, const arma::mat& x);
RcppExport SEXP _target_mlogit_pred(SEXP thetaSEXP, SEXP altSEXP, SEXP basealtSEXP, SEXP naltSEXP, SEXP id_idxSEXP, SEXP z1SEXP, SEXP z2SEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type alt(altSEXP);
    Rcpp::traits::input_parameter< unsigned >::type basealt(basealtSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nalt(naltSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_idx(id_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mlogit_pred(theta, alt, basealt, nalt, id_idx, z1, z2, x));
    return rcpp_result_gen;
END_RCPP
}
// mlogit_expand
Rcpp::List mlogit_expand(const arma::uvec& alt, const arma::mat& x, const arma::vec& weights, arma::uvec alts);
RcppExport SEXP _target_mlogit_expand(SEXP altSEXP, SEXP xSEXP, SEXP weightsSEXP, SEXP altsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type alt(altSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type alts(altsSEXP);
    rcpp_result_gen = Rcpp::wrap(mlogit_expand(alt, x, weights, alts));
    return rcpp_result_gen;
END_RCPP
}
// mlogit_obj
Rcpp::List mlogit_obj(arma::vec theta, const arma::uvec& choice, const arma::uvec& alt, unsigned basealt, unsigned nalt, const arma::uvec& id_idx, const arma::mat& z1, const arma::mat& z2, const arma::mat& x, const arma::vec& weights, bool return_hessian, bool onlyindiv);
RcppExport SEXP _target_mlogit_obj(SEXP thetaSEXP, SEXP choiceSEXP, SEXP altSEXP, SEXP basealtSEXP, SEXP naltSEXP, SEXP id_idxSEXP, SEXP z1SEXP, SEXP z2SEXP, SEXP xSEXP, SEXP weightsSEXP, SEXP return_hessianSEXP, SEXP onlyindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type choice(choiceSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type alt(altSEXP);
    Rcpp::traits::input_parameter< unsigned >::type basealt(basealtSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nalt(naltSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_idx(id_idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_hessian(return_hessianSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyindiv(onlyindivSEXP);
    rcpp_result_gen = Rcpp::wrap(mlogit_obj(theta, choice, alt, basealt, nalt, id_idx, z1, z2, x, weights, return_hessian, onlyindiv));
    return rcpp_result_gen;
END_RCPP
}
// bin_logl
arma::vec bin_logl(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, const arma::vec& weights, std::string type, bool indiv);
static SEXP _target_bin_logl_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_logl(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_logl(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_logl_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_dlogl
arma::mat bin_dlogl(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, const arma::vec& weights, std::string type, bool indiv);
static SEXP _target_bin_dlogl_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_dlogl(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_dlogl(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_dlogl_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_pa
arma::mat bin_pa(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, std::string type);
static SEXP _target_bin_pa_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_pa(y, a, x1, x2, par, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_pa(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_pa_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_dlogl_c
arma::cx_mat bin_dlogl_c(const arma::cx_vec& y, const arma::cx_vec& a, const arma::cx_mat& x1, const arma::cx_mat& x2, const arma::cx_vec par, const arma::cx_vec& weights, std::string type, bool indiv);
static SEXP _target_bin_dlogl_c_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_dlogl_c(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_dlogl_c(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_dlogl_c_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_esteq
arma::mat bin_esteq(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec& pr, arma::vec alpha, arma::vec par, const arma::vec& weights, std::string type);
static SEXP _target_bin_esteq_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP prSEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pr(prSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_esteq(y, a, x1, x2, pr, alpha, par, weights, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_esteq(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP prSEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_esteq_try(ySEXP, aSEXP, x1SEXP, x2SEXP, prSEXP, alphaSEXP, parSEXP, weightsSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_esteq_c
arma::cx_mat bin_esteq_c(const arma::cx_vec& y, const arma::cx_vec& a, const arma::cx_mat& x1, const arma::cx_mat& x2, const arma::cx_mat& x3, arma::cx_vec alpha, arma::cx_vec par, const arma::cx_vec& weights, std::string type);
static SEXP _target_bin_esteq_c_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP x3SEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x3(x3SEXP);
    Rcpp::traits::input_parameter< arma::cx_vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::cx_vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_esteq_c(y, a, x1, x2, x3, alpha, par, weights, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_bin_esteq_c(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP x3SEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_bin_esteq_c_try(ySEXP, aSEXP, x1SEXP, x2SEXP, x3SEXP, alphaSEXP, parSEXP, weightsSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ace_est
Rcpp::List ace_est(const arma::vec& y, const arma::mat& a, const arma::mat& x1, const arma::mat& x2, const arma::vec& theta, const arma::vec& weights, bool binary);
static SEXP _target_ace_est_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP thetaSEXP, SEXP weightsSEXP, SEXP binarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type binary(binarySEXP);
    rcpp_result_gen = Rcpp::wrap(ace_est(y, a, x1, x2, theta, weights, binary));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_ace_est(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP thetaSEXP, SEXP weightsSEXP, SEXP binarySEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_ace_est_try(ySEXP, aSEXP, x1SEXP, x2SEXP, thetaSEXP, weightsSEXP, binarySEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// fast_iid
arma::mat fast_iid(const arma::vec& y, const arma::vec& p, const arma::mat& x1, const arma::vec& weights, bool logistic);
static SEXP _target_fast_iid_try(SEXP ySEXP, SEXP pSEXP, SEXP x1SEXP, SEXP weightsSEXP, SEXP logisticSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type logistic(logisticSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_iid(y, p, x1, weights, logistic));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_fast_iid(SEXP ySEXP, SEXP pSEXP, SEXP x1SEXP, SEXP weightsSEXP, SEXP logisticSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_fast_iid_try(ySEXP, pSEXP, x1SEXP, weightsSEXP, logisticSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// clusterid
Rcpp::List clusterid(const arma::uvec& id);
static SEXP _target_clusterid_try(SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(clusterid(id));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_clusterid(SEXP idSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_clusterid_try(idSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// groupsum
arma::mat groupsum(const arma::mat& x, const arma::uvec& cluster, bool reduce);
static SEXP _target_groupsum_try(SEXP xSEXP, SEXP clusterSEXP, SEXP reduceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type cluster(clusterSEXP);
    Rcpp::traits::input_parameter< bool >::type reduce(reduceSEXP);
    rcpp_result_gen = Rcpp::wrap(groupsum(x, cluster, reduce));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_groupsum(SEXP xSEXP, SEXP clusterSEXP, SEXP reduceSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_groupsum_try(xSEXP, clusterSEXP, reduceSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// softmax
arma::mat softmax(arma::mat& lp, bool ref, bool log);
static SEXP _target_softmax_try(SEXP lpSEXP, SEXP refSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lp(lpSEXP);
    Rcpp::traits::input_parameter< bool >::type ref(refSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(softmax(lp, ref, log));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _target_softmax(SEXP lpSEXP, SEXP refSEXP, SEXP logSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_target_softmax_try(lpSEXP, refSEXP, logSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _target_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::vec(*bin_logl)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,const arma::vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_dlogl)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,const arma::vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_pa)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,std::string)");
        signatures.insert("arma::cx_mat(*bin_dlogl_c)(const arma::cx_vec&,const arma::cx_vec&,const arma::cx_mat&,const arma::cx_mat&,const arma::cx_vec,const arma::cx_vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_esteq)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,arma::vec,arma::vec,const arma::vec&,std::string)");
        signatures.insert("arma::cx_mat(*bin_esteq_c)(const arma::cx_vec&,const arma::cx_vec&,const arma::cx_mat&,const arma::cx_mat&,const arma::cx_mat&,arma::cx_vec,arma::cx_vec,const arma::cx_vec&,std::string)");
        signatures.insert("Rcpp::List(*ace_est)(const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&,bool)");
        signatures.insert("arma::mat(*fast_iid)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,bool)");
        signatures.insert("Rcpp::List(*.clusterid)(const arma::uvec&)");
        signatures.insert("arma::mat(*.groupsum)(const arma::mat&,const arma::uvec&,bool)");
        signatures.insert("arma::mat(*.softmax)(arma::mat&,bool,bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _target_RcppExport_registerCCallable() { 
    R_RegisterCCallable("target", "_target_bin_logl", (DL_FUNC)_target_bin_logl_try);
    R_RegisterCCallable("target", "_target_bin_dlogl", (DL_FUNC)_target_bin_dlogl_try);
    R_RegisterCCallable("target", "_target_bin_pa", (DL_FUNC)_target_bin_pa_try);
    R_RegisterCCallable("target", "_target_bin_dlogl_c", (DL_FUNC)_target_bin_dlogl_c_try);
    R_RegisterCCallable("target", "_target_bin_esteq", (DL_FUNC)_target_bin_esteq_try);
    R_RegisterCCallable("target", "_target_bin_esteq_c", (DL_FUNC)_target_bin_esteq_c_try);
    R_RegisterCCallable("target", "_target_ace_est", (DL_FUNC)_target_ace_est_try);
    R_RegisterCCallable("target", "_target_fast_iid", (DL_FUNC)_target_fast_iid_try);
    R_RegisterCCallable("target", "_target_.clusterid", (DL_FUNC)_target_clusterid_try);
    R_RegisterCCallable("target", "_target_.groupsum", (DL_FUNC)_target_groupsum_try);
    R_RegisterCCallable("target", "_target_.softmax", (DL_FUNC)_target_softmax_try);
    R_RegisterCCallable("target", "_target_RcppExport_validate", (DL_FUNC)_target_RcppExport_validate);
    return R_NilValue;
}

RcppExport SEXP _rcpp_module_boot_MLogit();
RcppExport SEXP _rcpp_module_boot_riskregmodels();

static const R_CallMethodDef CallEntries[] = {
    {"_target_mlogit_loglik", (DL_FUNC) &_target_mlogit_loglik, 10},
    {"_target_mlogit_pred", (DL_FUNC) &_target_mlogit_pred, 8},
    {"_target_mlogit_expand", (DL_FUNC) &_target_mlogit_expand, 4},
    {"_target_mlogit_obj", (DL_FUNC) &_target_mlogit_obj, 12},
    {"_target_bin_logl", (DL_FUNC) &_target_bin_logl, 8},
    {"_target_bin_dlogl", (DL_FUNC) &_target_bin_dlogl, 8},
    {"_target_bin_pa", (DL_FUNC) &_target_bin_pa, 6},
    {"_target_bin_dlogl_c", (DL_FUNC) &_target_bin_dlogl_c, 8},
    {"_target_bin_esteq", (DL_FUNC) &_target_bin_esteq, 9},
    {"_target_bin_esteq_c", (DL_FUNC) &_target_bin_esteq_c, 9},
    {"_target_ace_est", (DL_FUNC) &_target_ace_est, 7},
    {"_target_fast_iid", (DL_FUNC) &_target_fast_iid, 5},
    {"_target_clusterid", (DL_FUNC) &_target_clusterid, 1},
    {"_target_groupsum", (DL_FUNC) &_target_groupsum, 3},
    {"_target_softmax", (DL_FUNC) &_target_softmax, 3},
    {"_rcpp_module_boot_MLogit", (DL_FUNC) &_rcpp_module_boot_MLogit, 0},
    {"_rcpp_module_boot_riskregmodels", (DL_FUNC) &_rcpp_module_boot_riskregmodels, 0},
    {"_target_RcppExport_registerCCallable", (DL_FUNC) &_target_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_target(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
