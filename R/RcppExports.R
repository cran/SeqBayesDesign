# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

TMCMC <- function(dat, NN, initial, model, mu_fun, Cov, prior, priorDis, transp) {
    .Call(`_SeqBayesDesign_TMCMC`, dat, NN, initial, model, mu_fun, Cov, prior, priorDis, transp)
}

TMCMCALT <- function(dat, NN, initial, model, mu_fun, Cov, prior, priorDis, q) {
    .Call(`_SeqBayesDesign_TMCMCALT`, dat, NN, initial, model, mu_fun, Cov, prior, priorDis, q)
}

