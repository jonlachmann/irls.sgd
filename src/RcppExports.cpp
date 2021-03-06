// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// logisticG
Eigen::VectorXd logisticG(const Eigen::Map<Eigen::VectorXd>& theta, const Eigen::Map<Eigen::MatrixXd>& data, const Eigen::Map<Eigen::VectorXi>& idx);
RcppExport SEXP _irls_sgd_logisticG(SEXP thetaSEXP, SEXP dataSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(logisticG(theta, data, idx));
    return rcpp_result_gen;
END_RCPP
}
// gaussianG
Eigen::VectorXd gaussianG(const Eigen::Map<Eigen::VectorXd>& theta, const Eigen::Map<Eigen::MatrixXd>& data, const Eigen::Map<Eigen::VectorXi>& idx);
RcppExport SEXP _irls_sgd_gaussianG(SEXP thetaSEXP, SEXP dataSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussianG(theta, data, idx));
    return rcpp_result_gen;
END_RCPP
}
// qrls
Eigen::VectorXd qrls(const Eigen::Map<Eigen::MatrixXd>& x, const Eigen::Map<Eigen::VectorXd>& w, const Eigen::Map<Eigen::VectorXd>& z);
RcppExport SEXP _irls_sgd_qrls(SEXP xSEXP, SEXP wSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(qrls(x, w, z));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_irls_sgd_logisticG", (DL_FUNC) &_irls_sgd_logisticG, 3},
    {"_irls_sgd_gaussianG", (DL_FUNC) &_irls_sgd_gaussianG, 3},
    {"_irls_sgd_qrls", (DL_FUNC) &_irls_sgd_qrls, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_irls_sgd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
