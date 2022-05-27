// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ccd_update_factors_rcpp
arma::mat ccd_update_factors_rcpp(const arma::mat& V, const arma::mat& W, arma::mat& H, double e);
RcppExport SEXP _fastTopics_ccd_update_factors_rcpp(SEXP VSEXP, SEXP WSEXP, SEXP HSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_update_factors_rcpp(V, W, H, e));
    return rcpp_result_gen;
END_RCPP
}
// ccd_update_factors_sparse_rcpp
arma::mat ccd_update_factors_sparse_rcpp(const arma::sp_mat& V, const arma::mat& W, const arma::mat& H, double e);
RcppExport SEXP _fastTopics_ccd_update_factors_sparse_rcpp(SEXP VSEXP, SEXP WSEXP, SEXP HSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_update_factors_sparse_rcpp(V, W, H, e));
    return rcpp_result_gen;
END_RCPP
}
// ccd_update_factors_parallel_rcpp
arma::mat ccd_update_factors_parallel_rcpp(const arma::mat& V, const arma::mat& W, const arma::mat& H, double e);
RcppExport SEXP _fastTopics_ccd_update_factors_parallel_rcpp(SEXP VSEXP, SEXP WSEXP, SEXP HSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_update_factors_parallel_rcpp(V, W, H, e));
    return rcpp_result_gen;
END_RCPP
}
// ccd_update_factors_sparse_parallel_rcpp
arma::mat ccd_update_factors_sparse_parallel_rcpp(const arma::sp_mat& V, const arma::mat& W, const arma::mat& H, double e);
RcppExport SEXP _fastTopics_ccd_update_factors_sparse_parallel_rcpp(SEXP VSEXP, SEXP WSEXP, SEXP HSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_update_factors_sparse_parallel_rcpp(V, W, H, e));
    return rcpp_result_gen;
END_RCPP
}
// cost_rcpp
arma::vec cost_rcpp(const arma::mat& X, const arma::mat& A, const arma::mat& B, double e, bool poisson);
RcppExport SEXP _fastTopics_cost_rcpp(SEXP XSEXP, SEXP ASEXP, SEXP BSEXP, SEXP eSEXP, SEXP poissonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< bool >::type poisson(poissonSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_rcpp(X, A, B, e, poisson));
    return rcpp_result_gen;
END_RCPP
}
// cost_sparse_rcpp
arma::vec cost_sparse_rcpp(const arma::sp_mat& X, const arma::mat& A, const arma::mat& B, double e, bool poisson);
RcppExport SEXP _fastTopics_cost_sparse_rcpp(SEXP XSEXP, SEXP ASEXP, SEXP BSEXP, SEXP eSEXP, SEXP poissonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< bool >::type poisson(poissonSEXP);
    rcpp_result_gen = Rcpp::wrap(cost_sparse_rcpp(X, A, B, e, poisson));
    return rcpp_result_gen;
END_RCPP
}
// le_diff_rcpp
arma::mat le_diff_rcpp(const arma::mat& X);
RcppExport SEXP _fastTopics_le_diff_rcpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(le_diff_rcpp(X));
    return rcpp_result_gen;
END_RCPP
}
// x_over_crossprod_rcpp
arma::vec x_over_crossprod_rcpp(const arma::vec& i, const arma::vec& j, const arma::vec& x, const arma::mat& A, const arma::mat& B, double e);
RcppExport SEXP _fastTopics_x_over_crossprod_rcpp(SEXP iSEXP, SEXP jSEXP, SEXP xSEXP, SEXP ASEXP, SEXP BSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(x_over_crossprod_rcpp(i, j, x, A, B, e));
    return rcpp_result_gen;
END_RCPP
}
// mixem_rcpp
arma::vec mixem_rcpp(const arma::mat& L, const arma::vec& w, const arma::vec& x0, unsigned int numiter);
RcppExport SEXP _fastTopics_mixem_rcpp(SEXP LSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(mixem_rcpp(L, w, x0, numiter));
    return rcpp_result_gen;
END_RCPP
}
// pnmfem_update_factors_rcpp
arma::mat pnmfem_update_factors_rcpp(const arma::mat& X, const arma::mat& F, const arma::mat& L, const arma::vec& j, double numiter);
RcppExport SEXP _fastTopics_pnmfem_update_factors_rcpp(SEXP XSEXP, SEXP FSEXP, SEXP LSEXP, SEXP jSEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(pnmfem_update_factors_rcpp(X, F, L, j, numiter));
    return rcpp_result_gen;
END_RCPP
}
// pnmfem_update_factors_sparse_rcpp
arma::mat pnmfem_update_factors_sparse_rcpp(const arma::sp_mat& X, const arma::mat& F, const arma::mat& L, const arma::vec& j, double numiter);
RcppExport SEXP _fastTopics_pnmfem_update_factors_sparse_rcpp(SEXP XSEXP, SEXP FSEXP, SEXP LSEXP, SEXP jSEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(pnmfem_update_factors_sparse_rcpp(X, F, L, j, numiter));
    return rcpp_result_gen;
END_RCPP
}
// pnmfem_update_factors_parallel_rcpp
arma::mat pnmfem_update_factors_parallel_rcpp(const arma::mat& X, const arma::mat& F, const arma::mat& L, const arma::vec& j, double numiter);
RcppExport SEXP _fastTopics_pnmfem_update_factors_parallel_rcpp(SEXP XSEXP, SEXP FSEXP, SEXP LSEXP, SEXP jSEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(pnmfem_update_factors_parallel_rcpp(X, F, L, j, numiter));
    return rcpp_result_gen;
END_RCPP
}
// pnmfem_update_factors_sparse_parallel_rcpp
arma::mat pnmfem_update_factors_sparse_parallel_rcpp(const arma::sp_mat& X, const arma::mat& F, const arma::mat& L, const arma::vec& j, double numiter);
RcppExport SEXP _fastTopics_pnmfem_update_factors_sparse_parallel_rcpp(SEXP XSEXP, SEXP FSEXP, SEXP LSEXP, SEXP jSEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(pnmfem_update_factors_sparse_parallel_rcpp(X, F, L, j, numiter));
    return rcpp_result_gen;
END_RCPP
}
// poismixem_rcpp
arma::vec poismixem_rcpp(const arma::mat& L, const arma::vec& w, const arma::vec& x0, unsigned int numiter);
RcppExport SEXP _fastTopics_poismixem_rcpp(SEXP LSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(poismixem_rcpp(L, w, x0, numiter));
    return rcpp_result_gen;
END_RCPP
}
// poismixem2_rcpp
arma::vec poismixem2_rcpp(const arma::mat& L1, const arma::vec& w, const arma::vec& u, const arma::vec& x0, unsigned int numiter);
RcppExport SEXP _fastTopics_poismixem2_rcpp(SEXP L1SEXP, SEXP wSEXP, SEXP uSEXP, SEXP x0SEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(poismixem2_rcpp(L1, w, u, x0, numiter));
    return rcpp_result_gen;
END_RCPP
}
// poismixem3_rcpp
arma::vec poismixem3_rcpp(const arma::mat& L1, const arma::vec& w, const arma::vec& u, const arma::uvec& i, const arma::vec& x0, unsigned int numiter);
RcppExport SEXP _fastTopics_poismixem3_rcpp(SEXP L1SEXP, SEXP wSEXP, SEXP uSEXP, SEXP iSEXP, SEXP x0SEXP, SEXP numiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    rcpp_result_gen = Rcpp::wrap(poismixem3_rcpp(L1, w, u, i, x0, numiter));
    return rcpp_result_gen;
END_RCPP
}
// scd_kl_update_rcpp
arma::vec scd_kl_update_rcpp(const arma::mat& L, const arma::vec& w, const arma::vec& x0, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_kl_update_rcpp(SEXP LSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_kl_update_rcpp(L, w, x0, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// scd_kl_update2_rcpp
arma::vec scd_kl_update2_rcpp(const arma::mat& L, const arma::vec& u, const arma::vec& w, const arma::vec& x0, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_kl_update2_rcpp(SEXP LSEXP, SEXP uSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_kl_update2_rcpp(L, u, w, x0, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// ccd_kl_update_rcpp
arma::vec ccd_kl_update_rcpp(const arma::mat& L, const arma::vec& w, const arma::vec& x0, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_ccd_kl_update_rcpp(SEXP LSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_kl_update_rcpp(L, w, x0, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// ccd_kl_update2_rcpp
arma::vec ccd_kl_update2_rcpp(const arma::mat& L, const arma::vec& u, const arma::vec& w, const arma::vec& x0, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_ccd_kl_update2_rcpp(SEXP LSEXP, SEXP uSEXP, SEXP wSEXP, SEXP x0SEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(ccd_kl_update2_rcpp(L, u, w, x0, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// simulate_posterior_poisson_rcpp
List simulate_posterior_poisson_rcpp(const arma::vec& x, const arma::mat& L, const arma::vec& f, const arma::mat& D, const arma::mat& U, const arma::mat& M, double s, double e);
RcppExport SEXP _fastTopics_simulate_posterior_poisson_rcpp(SEXP xSEXP, SEXP LSEXP, SEXP fSEXP, SEXP DSEXP, SEXP USEXP, SEXP MSEXP, SEXP sSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_posterior_poisson_rcpp(x, L, f, D, U, M, s, e));
    return rcpp_result_gen;
END_RCPP
}
// simulate_posterior_poisson_sparse_rcpp
List simulate_posterior_poisson_sparse_rcpp(const arma::vec& x, const arma::mat& L, const arma::vec& w, const arma::vec& f, const arma::mat& D, const arma::mat& U, const arma::mat& M, double s, double e);
RcppExport SEXP _fastTopics_simulate_posterior_poisson_sparse_rcpp(SEXP xSEXP, SEXP LSEXP, SEXP wSEXP, SEXP fSEXP, SEXP DSEXP, SEXP USEXP, SEXP MSEXP, SEXP sSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_posterior_poisson_sparse_rcpp(x, L, w, f, D, U, M, s, e));
    return rcpp_result_gen;
END_RCPP
}
// scd_update_factors_rcpp
arma::mat scd_update_factors_rcpp(const arma::mat& A, const arma::mat& W, const arma::mat& H, const arma::vec& j, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_update_factors_rcpp(SEXP ASEXP, SEXP WSEXP, SEXP HSEXP, SEXP jSEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_update_factors_rcpp(A, W, H, j, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// scd_update_factors_sparse_rcpp
arma::mat scd_update_factors_sparse_rcpp(const arma::sp_mat& A, const arma::mat& W, const arma::mat& H, const arma::vec& j, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_update_factors_sparse_rcpp(SEXP ASEXP, SEXP WSEXP, SEXP HSEXP, SEXP jSEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_update_factors_sparse_rcpp(A, W, H, j, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// scd_update_factors_parallel_rcpp
arma::mat scd_update_factors_parallel_rcpp(const arma::mat& A, const arma::mat& W, const arma::mat& H, const arma::vec& j, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_update_factors_parallel_rcpp(SEXP ASEXP, SEXP WSEXP, SEXP HSEXP, SEXP jSEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_update_factors_parallel_rcpp(A, W, H, j, numiter, e));
    return rcpp_result_gen;
END_RCPP
}
// scd_update_factors_sparse_parallel_rcpp
arma::mat scd_update_factors_sparse_parallel_rcpp(const arma::sp_mat& A, const arma::mat& W, const arma::mat& H, const arma::vec& j, unsigned int numiter, double e);
RcppExport SEXP _fastTopics_scd_update_factors_sparse_parallel_rcpp(SEXP ASEXP, SEXP WSEXP, SEXP HSEXP, SEXP jSEXP, SEXP numiterSEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(scd_update_factors_sparse_parallel_rcpp(A, W, H, j, numiter, e));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastTopics_ccd_update_factors_rcpp", (DL_FUNC) &_fastTopics_ccd_update_factors_rcpp, 4},
    {"_fastTopics_ccd_update_factors_sparse_rcpp", (DL_FUNC) &_fastTopics_ccd_update_factors_sparse_rcpp, 4},
    {"_fastTopics_ccd_update_factors_parallel_rcpp", (DL_FUNC) &_fastTopics_ccd_update_factors_parallel_rcpp, 4},
    {"_fastTopics_ccd_update_factors_sparse_parallel_rcpp", (DL_FUNC) &_fastTopics_ccd_update_factors_sparse_parallel_rcpp, 4},
    {"_fastTopics_cost_rcpp", (DL_FUNC) &_fastTopics_cost_rcpp, 5},
    {"_fastTopics_cost_sparse_rcpp", (DL_FUNC) &_fastTopics_cost_sparse_rcpp, 5},
    {"_fastTopics_le_diff_rcpp", (DL_FUNC) &_fastTopics_le_diff_rcpp, 1},
    {"_fastTopics_x_over_crossprod_rcpp", (DL_FUNC) &_fastTopics_x_over_crossprod_rcpp, 6},
    {"_fastTopics_mixem_rcpp", (DL_FUNC) &_fastTopics_mixem_rcpp, 4},
    {"_fastTopics_pnmfem_update_factors_rcpp", (DL_FUNC) &_fastTopics_pnmfem_update_factors_rcpp, 5},
    {"_fastTopics_pnmfem_update_factors_sparse_rcpp", (DL_FUNC) &_fastTopics_pnmfem_update_factors_sparse_rcpp, 5},
    {"_fastTopics_pnmfem_update_factors_parallel_rcpp", (DL_FUNC) &_fastTopics_pnmfem_update_factors_parallel_rcpp, 5},
    {"_fastTopics_pnmfem_update_factors_sparse_parallel_rcpp", (DL_FUNC) &_fastTopics_pnmfem_update_factors_sparse_parallel_rcpp, 5},
    {"_fastTopics_poismixem_rcpp", (DL_FUNC) &_fastTopics_poismixem_rcpp, 4},
    {"_fastTopics_poismixem2_rcpp", (DL_FUNC) &_fastTopics_poismixem2_rcpp, 5},
    {"_fastTopics_poismixem3_rcpp", (DL_FUNC) &_fastTopics_poismixem3_rcpp, 6},
    {"_fastTopics_scd_kl_update_rcpp", (DL_FUNC) &_fastTopics_scd_kl_update_rcpp, 5},
    {"_fastTopics_scd_kl_update2_rcpp", (DL_FUNC) &_fastTopics_scd_kl_update2_rcpp, 6},
    {"_fastTopics_ccd_kl_update_rcpp", (DL_FUNC) &_fastTopics_ccd_kl_update_rcpp, 5},
    {"_fastTopics_ccd_kl_update2_rcpp", (DL_FUNC) &_fastTopics_ccd_kl_update2_rcpp, 6},
    {"_fastTopics_simulate_posterior_poisson_rcpp", (DL_FUNC) &_fastTopics_simulate_posterior_poisson_rcpp, 8},
    {"_fastTopics_simulate_posterior_poisson_sparse_rcpp", (DL_FUNC) &_fastTopics_simulate_posterior_poisson_sparse_rcpp, 9},
    {"_fastTopics_scd_update_factors_rcpp", (DL_FUNC) &_fastTopics_scd_update_factors_rcpp, 6},
    {"_fastTopics_scd_update_factors_sparse_rcpp", (DL_FUNC) &_fastTopics_scd_update_factors_sparse_rcpp, 6},
    {"_fastTopics_scd_update_factors_parallel_rcpp", (DL_FUNC) &_fastTopics_scd_update_factors_parallel_rcpp, 6},
    {"_fastTopics_scd_update_factors_sparse_parallel_rcpp", (DL_FUNC) &_fastTopics_scd_update_factors_sparse_parallel_rcpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastTopics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}