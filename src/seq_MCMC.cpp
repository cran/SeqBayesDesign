#include <Rcpp.h>
#include <stdio.h>

using namespace Rcpp;

NumericVector StandStressFun(double A, double B, NumericVector StandLevels, CharacterVector model, List data)
{
  std::string cmu = Rcpp::as<std::string>(model);
  std::string cmu1 = "linear";
  std::string cmu2 = "Arrhenius";
  std::string cmu3 = "inverse-power";
  std::string cmu4 = "E-C";
  NumericVector mu;
  if((cmu == cmu1) | (cmu == cmu2) | (cmu == cmu3)) {
    mu = A + B*StandLevels;
  }else if(cmu == cmu4) {
    double phi_R, gamma_theta;
    NumericVector ss, cc, DD;
    double Rval=data["Rval"];
    if(Rval>1){phi_R=1/Rval;}else{phi_R=Rval;}
    double angle=data["Angle"];
    NumericVector ff=data["Freq"];
    gamma_theta=1.6-phi_R*sin(angle);
    ss=1/StandLevels;
    cc=(ss-1)*pow(ss,(gamma_theta-1))*pow((1-phi_R),(-gamma_theta));
    DD=B*pow(ff,B)*cc+A;
    mu=(log(DD)-log(A))/B;
  }
  return mu;
}

NumericVector VecDnorm(NumericVector x, NumericVector means, double sd)
{
  int n = x.size();
  NumericVector res(n);
  for(int i=0; i<n; i++) {
    double q = x[i];
    double mu = means[i];
    res[i] = R::dnorm(q, mu, sd, 0);
  }
  return res;
}


NumericVector VecPnorm(NumericVector x, NumericVector means, double sd)
{
  int n = x.size();
  NumericVector res(n);
  for(int i=0; i<n; i++) {
    double q = x[i];
    double mu = means[i];
    res[i] = R::pnorm(q, mu, sd, 1, 0);
  }
  return res;
}

NumericVector dsev(NumericVector z)
{
  int n = z.size();
  NumericVector res(n);
  for(int i=0; i<n; i++) {
    res[i] = exp(z[i]-exp(z[i]));
  }
  return res;
}

NumericVector psev(NumericVector z)
{
  int n = z.size();
  NumericVector res(n);
  for(int i=0; i<n; i++) {
    res[i] = 1-exp(-exp(z[i]));
    if(res[i] > 0.99999999999) {
      res[i] = 0.99999999999;
    }
  }
  return res;
}


NumericVector StandLevel(List use_level, double max_level, NumericVector test_level, CharacterVector mu_fun)
{
  std::string cmu = Rcpp::as<std::string>(mu_fun);
  std::string cmu1 = "linear";
  std::string cmu2 = "Arrhenius";
  std::string cmu3 = "inverse-power";
  std::string cmu4 = "E-C";
  NumericVector res, test;
  test = test_level;
  NumericVector use = use_level["use"];
  if(cmu == cmu1) {
    res=(test-use[0])/(max_level-use[0]);
  }else if(cmu == cmu2) {
    res=(1/(test+273.15)-1/(use[0] +273.15))/(1/(max_level+273.15)-1/(use[0]+273.15));
  }else if(cmu == cmu3) {
    res=(log(test)-log(use[0]))/(log(max_level)-log(use[0]));
  }else if(cmu == cmu4) {
    res=test/max_level;
  }
  return res;
}


double likelihood(List dat, NumericVector pars, CharacterVector model, CharacterVector mu_fun)
{
  std::string cy = Rcpp::as<std::string>(model);
  std::string cy1 = "lnor";
  std::string cy2 = "wei";
  double A, B, nu;
  A = pars[0];
  B = pars[1];
  nu = pars[2];
  NumericVector stlevel, mu;

  NumericVector stress = dat["x"];
  double max = dat["max.level"];
  List use = dat["use.level"];

  stlevel = StandLevel(use, max, stress, mu_fun);
  mu = StandStressFun(A, B, stlevel, mu_fun, dat);

  NumericVector wts = dat["wts"];
  NumericVector Cen = dat["Censored"];
  NumericVector y = dat["Y"];
  NumericVector yy = log(y);
  NumericVector zz, ff, ff1, FF, normcdf, normpdf, ll;
  double res;
  res = NA_REAL;
  if (cy == cy1) {
    normcdf = ifelse(VecPnorm(yy, mu, nu) >= 0.99999999, 0.99999999, VecPnorm(yy, mu, nu));
    normcdf = ifelse(normcdf <= 1-0.99999999, 1-0.99999999, normcdf);
    normpdf = ifelse(VecDnorm(yy, mu, nu)/(y) <= 1-0.99999999, 1-0.99999999, VecDnorm(yy, mu, nu)/(y));
    ll = log(normpdf)*(1-Cen)+log(1-normcdf)*Cen;
    ll = wts*ll;
    res = exp(std::accumulate(ll.begin(), ll.end(), 0.0));
  }else if(cy == cy2) {
    zz = (yy-mu)/nu;
    ff = dsev(zz)/(nu*y);
    ff1 = ifelse(ff <= 1-0.99999999, 1-0.99999999, ff);
    FF = psev(zz);
    ll = wts*(log(ff1)*(1-Cen)+log(1-FF)*Cen);
    res = exp(std::accumulate(ll.begin(), ll.end(), 0.0));
  }
  return(res);
}

// [[Rcpp::export]]
List TMCMC(List dat, int NN, NumericVector initial, CharacterVector model, CharacterVector mu_fun, NumericMatrix Cov, NumericVector prior, CharacterVector priorDis, double transp)
{
  std::string pd = Rcpp::as<std::string>(priorDis);
  std::string pdnorm = "normal";
  std::string pdunif = "uniform";

  List res;
  NumericMatrix par(3, NN);
  NumericMatrix acp(2, NN);
  par(_, 0) = initial;
  double kapa = prior[4];
  double gamma = prior[5];

  //List data = dat["data"];
  NumericVector y = dat["Y"];
  NumericVector Cen = dat["Censored"];
  double Nobs = dat["Nobs"];
  double Ncen = dat["Ncen"];
  NumericVector wts=dat["wts"];
  NumericVector stress = dat["x"];
  NumericVector stlevel, mu;
  List use = dat["use.level"];
  stlevel = StandLevel(use, dat["max.level"], stress, mu_fun);
  for (int i=1; i< NN ; i++) {
    Rcout<<i+1<<std::endl;

    double A = par(0, i-1);
    double B = par(1, i-1);
    double nu = par(2, i-1);

    double sig_sqrtA = sqrt(Cov(0, 0));
    double sig_sqrtB = sqrt(Cov(1, 1));
    double rho = Cov(0, 1)/(sig_sqrtA*sig_sqrtB);

    // A and B //
    double z1 = rnorm(1, 0.0, 1.0)[0];
    double z2 = rnorm(1, 0.0, 1.0)[0];
    double pro2_A = log(B/A) + z1*sig_sqrtA*3;
    double pro2_B = pow(B, transp) + (rho*z1 + sqrt(1-rho*rho)*z2)*sig_sqrtB*3;

    double proA0=pow(pro2_B, 1/transp)*exp(-pro2_A);
    double rr1, rprior1;
    rr1 = NA_REAL;
    rprior1 = NA_REAL;
    if (pd == pdunif) {
      while( proA0 < prior[0] || proA0 > prior[1] ||
             pro2_B < pow(prior[2], transp) || pro2_B > pow(prior[3], transp)) {
        z1 = rnorm(1, 0, 1)[0];
        z2 = rnorm(1, 0, 1)[0];
        pro2_A = log(B/A) + z1*sig_sqrtA*3;
        pro2_B = pow(B, transp) + (rho*z1 + sqrt(1-rho*rho)*z2)*sig_sqrtB*3;
        proA0 = pow(pro2_B, 1/transp)*exp(-pro2_A);
      }
      NumericVector pro1(3);
      pro1[0] = pow(pro2_B, 1/transp)*exp(-pro2_A);
      pro1[1] = pow(pro2_B, 1/transp);
      pro1[2]=nu;
      rr1 = likelihood(dat, pro1, model, mu_fun)/likelihood(dat, par(_, i-1),
                       model, mu_fun);
      rprior1 = 1;
    } else if (pd == pdnorm) {
      NumericVector pro1(3);
      pro1[0] = pow(pro2_B, 1/transp)*exp(-pro2_A);
      pro1[1] = pow(pro2_B, 1/transp);
      pro1[2] = nu;
      rr1 = likelihood(dat, pro1, model, mu_fun)/likelihood(dat, par(_, i-1), model, mu_fun);
      rprior1 = (R::dnorm(pro1[0], prior[0], prior[1], 0)*R::dnorm(pro1[1], prior[2], prior[3], 0))/(R::dnorm(par(0, i-1), prior[0], prior[1], 0)*R::dnorm(par(1, i-1), prior[2], prior[3], 0));
    }
    double p1 = std::min(rr1*rprior1, 1.0);
    double b1 = rbinom(1, 1, p1)[0];
    if(b1 == 1) {
      par(0, i) = pow(pro2_B, 1/transp)*exp(-pro2_A);
      par(1, i) = pow(pro2_B, 1/transp);
      acp(0, i) = 1;
    }else {
      par(0, i) = A;
      par(1, i) = B;
      acp(0, i) = 0;
    }
    mu = StandStressFun(par(0, i), par(1, i), stlevel, mu_fun, dat);

    NumericVector scale1 = (log(y)-mu)*(log(y)-mu)*(1-Cen)*wts;
    double scale2 = std::accumulate(scale1.begin(), scale1.end(),0.0)/2+gamma;
    double shape = (Nobs-Ncen)/2+kapa;
    double nu2 = sqrt(1/rgamma(1, shape, 1/scale2)[0]);

    NumericVector zz0 = (log(y)-mu)/nu;
    NumericVector normcdf0 = ifelse(VecPnorm(zz0, 0, 1)>=0.99999999, 0.99999999, VecPnorm(zz0, 0, 1));
    normcdf0 = ifelse(normcdf0 <= 1-0.99999999, 1-0.99999999, normcdf0);
    NumericVector zz2 = (log(y)-mu)/nu2;
    NumericVector normcdf2 = ifelse(VecPnorm(zz2, 0, 1)>=0.99999999, 0.99999999, VecPnorm(zz2, 0, 1));
    normcdf2 = ifelse(normcdf2 <= 1-0.99999999, 1-0.99999999, normcdf2);
    NumericVector RR= wts*Cen*(log(1-normcdf2)-log(1-normcdf0));
    double r1 = exp(std::accumulate(RR.begin(), RR.end(), 0.0));
    double p2 = std::min(r1, 1.0);
    double b2 = rbinom(1, 1, p2)[0];
    if(b2 == 1) {
      par(2,i) = nu2;
      acp(1 ,i) = 1;
    }else {
      par(2,i) = nu;
      acp(1 ,i) = 0;
    }
  }
  res["par"] = par;
  res["acp"] = acp;
  return(res);
}


// [[Rcpp::export]]
List TMCMCALT(List dat, int NN, NumericVector initial, CharacterVector model, CharacterVector mu_fun, NumericMatrix Cov, NumericVector prior, CharacterVector priorDis, double q)
{
  std::string pd = Rcpp::as<std::string>(priorDis);
  std::string pdnorm = "normal";
  std::string pdunif = "uniform";

  List res;
  NumericMatrix par(3, NN);
  NumericMatrix acp(3, NN);
  par(_,0) = initial;
  double kapa = prior[4];
  double gamma = prior[5];

  NumericVector wts = dat["wts"];
  NumericVector y = dat["Y"];
  NumericVector Cen = dat["Censored"];
  int Nobs = dat["Nobs"];
  int Ncen = dat["Ncen"];
  double zq = R::qnorm(q, 0.0, 1.0, 1, 0);
  NumericVector stress = dat["x"];
  NumericVector stlevel, mu;
  List use = dat["use.level"];
  stlevel = StandLevel(use, dat["max.level"], stress, mu_fun);
  double pro2_A, proA0;
  for (int i=1; i< NN ; i++) {
    Rcout<<i+1<<std::endl;

    double A = par(0, i-1);
    double B = par(1, i-1);
    double nu = par(2, i-1);

    double sig_sqrtA = sqrt(Cov(0,0) + zq*zq*Cov(2,2) + Cov(0,2)*zq);
    double z1 = rnorm(1, 0, 1)[0];
    pro2_A = A + zq*nu + z1*sig_sqrtA;
    proA0 = pro2_A - zq*nu;
    double rr1, rprior1;
    rr1 = NA_REAL;
    rprior1 = NA_REAL;
    if (pd == pdunif) {
      while(proA0 < prior[0] || proA0 > prior[1]) {
        z1 = rnorm(1, 0, 1)[0];
        pro2_A = A + zq*nu + z1*sig_sqrtA;
        proA0 = pro2_A - zq*nu;
      }
      NumericVector pro1(3);
      pro1[0] = proA0;
      pro1[1] = B;
      pro1[2] = nu;
      rr1 = likelihood(dat, pro1, model, mu_fun)/likelihood(dat, par(_, i-1), model, mu_fun);
      rprior1 = 1;
    } else if (pd == pdnorm) {
      NumericVector pro1(3);
      pro1[0] = proA0;
      pro1[1] = B;
      pro1[2] = nu;
      rr1 = likelihood(dat, pro1, model, mu_fun)/likelihood(dat, par(_, i-1), model, mu_fun);
      rprior1 = R::dnorm(pro1[0], prior[0], prior[1], 0)/R::dnorm(par(0, i-1), prior[0], prior[1], 0);
    }
    double p1 = std::min(rr1*rprior1, 1.0);
    double b1 = rbinom(1, 1, p1)[0];
    if(b1 == 1) {
      par(0, i) = proA0;
      acp(0, i)=1;
    }else {
      par(0, i)=A;
      acp(0, i)=0;
    }
    double sig_sqrtB = sqrt(Cov(1, 1));
    double z2 = rnorm(1, 0, 1)[0];
    double pro2_B = B + z2*sig_sqrtB;
    NumericVector pro20(3);
    pro20[0] = par(0, i);
    pro20[1] = B;
    pro20[2] = nu;
    double rr2, rprior2;
    rr2 = NA_REAL;
    rprior2 = NA_REAL;
    if (pd == pdunif) {
      while((pro2_B < prior[2]) || (pro2_B > prior[3])) {
        z2 = rnorm(1, 0, 1)[0];
        pro2_B = B + z2*sig_sqrtB;
      }
      NumericVector pro2(3);
      pro2[0] = par(0, i);
      pro2[1] = pro2_B;
      pro2[2] = nu;
      rr2 = likelihood(dat, pro2, model, mu_fun)/likelihood(dat, pro20, model, mu_fun);
      rprior2 = 1;
    } else if (pd == pdnorm) {
      NumericVector pro2(3);
      pro2[0] = par(0, i);
      pro2[1] = pro2_B;
      pro2[2] = nu;
      rr2 = likelihood(dat, pro2, model, mu_fun)/likelihood(dat, pro20, model, mu_fun);
      rprior2 = R::dnorm(pro2[1], prior[2], prior[3], 0)/R::dnorm(par(1, i-1), prior[2], prior[3], 0);
    }

    double p2 = std::min(rr2*rprior2, 1.0);
    double b2 = rbinom(1, 1, p2)[0];
    if(b2 == 1) {
      par(1, i) = pro2_B;
      acp(1, i) = 1;
    }else {
      par(1, i) = B;
      acp(1, i) = 0;
    }
    mu = StandStressFun(par(0, i), par(1, i), stlevel, mu_fun, dat);
    NumericVector scale1 = (log(y)-mu)*(log(y)-mu)*(1-Cen)*wts;
    double scale2 = std::accumulate(scale1.begin(), scale1.end(), 0.0)/2+gamma;
    double shape = (Nobs-Ncen)/2+kapa;
    double nu2 = sqrt(1/rgamma(1, shape, 1/scale2)[0]);

    NumericVector zz0 = (log(y)-mu)/nu;
    NumericVector normcdf0 = ifelse(VecPnorm(zz0, 0, 1)>=0.99999999, 0.99999999, VecPnorm(zz0, 0, 1));
    normcdf0 = ifelse(normcdf0 <= 1-0.99999999, 1-0.99999999, normcdf0);
    NumericVector zz2 = (log(y)-mu)/nu2;
    NumericVector normcdf2 = ifelse(VecPnorm(zz2, 0, 1)>=0.99999999, 0.99999999, VecPnorm(zz2, 0, 1));
    normcdf2 = ifelse(normcdf2 <= 1-0.99999999, 1-0.99999999, normcdf2);

    NumericVector RR = wts*Cen*(log(1-normcdf2)-log(1-normcdf0));
    double r3 = exp(std::accumulate(RR.begin(), RR.end(),0.0));
    double p3 = std::min(r3,1.0);
    double b3 = rbinom(1, 1, p3)[0];
    if(b3 == 1) {
      par(2,i) = nu2;
      acp(2 ,i) = 1;
    }else {
      par(2,i) = nu;
      acp(2 ,i) = 0;
    }
  }
  res["par"]=par;
  res["acp"]=acp;
  return(res);
}

