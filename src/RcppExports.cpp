// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// TMCMC
List TMCMC(List dat, int NN, NumericVector initial, CharacterVector model, CharacterVector mu_fun, NumericMatrix Cov, NumericVector prior, CharacterVector priorDis, double transp);
RcppExport SEXP _SeqBayesDesign_TMCMC(SEXP datSEXP, SEXP NNSEXP, SEXP initialSEXP, SEXP modelSEXP, SEXP mu_funSEXP, SEXP CovSEXP, SEXP priorSEXP, SEXP priorDisSEXP, SEXP transpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type mu_fun(mu_funSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type priorDis(priorDisSEXP);
    Rcpp::traits::input_parameter< double >::type transp(transpSEXP);
    rcpp_result_gen = Rcpp::wrap(TMCMC(dat, NN, initial, model, mu_fun, Cov, prior, priorDis, transp));
    return rcpp_result_gen;
END_RCPP
}
// TMCMCALT
List TMCMCALT(List dat, int NN, NumericVector initial, CharacterVector model, CharacterVector mu_fun, NumericMatrix Cov, NumericVector prior, CharacterVector priorDis, double q);
RcppExport SEXP _SeqBayesDesign_TMCMCALT(SEXP datSEXP, SEXP NNSEXP, SEXP initialSEXP, SEXP modelSEXP, SEXP mu_funSEXP, SEXP CovSEXP, SEXP priorSEXP, SEXP priorDisSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type mu_fun(mu_funSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type priorDis(priorDisSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(TMCMCALT(dat, NN, initial, model, mu_fun, Cov, prior, priorDis, q));
    return rcpp_result_gen;
END_RCPP
}
