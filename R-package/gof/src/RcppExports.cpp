// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/gof.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// SupTest
double SupTest(const arma::vec& x);
static SEXP _gof_SupTest_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(SupTest(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _gof_SupTest(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_gof_SupTest_try(xSEXP));
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
// L2Test
double L2Test(const arma::vec& x, const arma::vec& t);
static SEXP _gof_L2Test_try(SEXP xSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(L2Test(x, t));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _gof_L2Test(SEXP xSEXP, SEXP tSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_gof_L2Test_try(xSEXP, tSEXP));
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
static int _gof_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*SupTest)(const arma::vec&)");
        signatures.insert("double(*L2Test)(const arma::vec&,const arma::vec&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _gof_RcppExport_registerCCallable() { 
    R_RegisterCCallable("gof", "_gof_SupTest", (DL_FUNC)_gof_SupTest_try);
    R_RegisterCCallable("gof", "_gof_L2Test", (DL_FUNC)_gof_L2Test_try);
    R_RegisterCCallable("gof", "_gof_RcppExport_validate", (DL_FUNC)_gof_RcppExport_validate);
    return R_NilValue;
}

RcppExport SEXP _rcpp_module_boot_gofmod();

static const R_CallMethodDef CallEntries[] = {
    {"_gof_SupTest", (DL_FUNC) &_gof_SupTest, 1},
    {"_gof_L2Test", (DL_FUNC) &_gof_L2Test, 2},
    {"_rcpp_module_boot_gofmod", (DL_FUNC) &_rcpp_module_boot_gofmod, 0},
    {"_gof_RcppExport_registerCCallable", (DL_FUNC) &_gof_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gof(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}