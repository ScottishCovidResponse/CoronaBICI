// This file samples from or calculates probailities for distributions (Exp, gamma etc..)

const double shmin = 1, shmax = 20;
const long nsh = 400;
const double gmax = 10;
const long ng = 400;

double tgammaprob(double x, double a, double b);

vector <vector <double> > gammaupst;

double Chain::likenm(long i)         // The likelikelihood of non-Markovian transitions in individual i
{
  long c, e, j, cl;
  long tr, tr2;
  double t, L = 0, mean, kshape, lam;

  c = tra[indev[i][0].tr].cf; 
  for(j = 0; j < nclnm; j++){ cl = clnm[j]; cstnm[cl] = c; tstnm[cl] = indev[i][0].t;}

  for(e = 1; e < nindev[i]; e++){
    tr = indev[i][e].tr; t = indev[i][e].t;
    if(tra[tr].nm == 1){
      cl = tra[tr].cl;
      if(tra[tr].like == 1){
        switch(tra[tr].type){
          case GAMMA_TR:
            mean = nmeq_val[tra[tr].eq];
            kshape = nmeq_val[tra[tr].eqshape];
            L += gammaprob(t-tstnm[cl],kshape,kshape/mean);
            break;

          case WEI_TR:
            lam = nmeq_val[tra[tr].eq];
            kshape = nmeq_val[tra[tr].eqshape];
            L += weibullprob(t-tstnm[cl],lam,kshape);
            break;
        }
      }
      else{  // This is added for the "move" data
        switch(tra[tr].type){
          case GAMMA_TR:
            mean = nmeq_val[tra[tr].eq];
            kshape = nmeq_val[tra[tr].eqshape];
            L += gammaup(t-tstnm[cl],kshape,kshape/mean);
            break;

          case WEI_TR:
            lam = nmeq_val[tra[tr].eq];
            kshape = nmeq_val[tra[tr].eqshape];
            L += weibullup(t-tstnm[cl],lam,kshape);
            break;
        } 
      }
    }

    for(j = 0; j < tra[tr].ntraend; j++){
      tr2 = tra[tr].traend[j]; if(tr2 < 0) emsg("Dist: EC1");
      cl = tra[tr2].cl;
      c = cstnm[cl];

      switch(tra[tr2].type){
        case GAMMA_TR:
          mean = nmeq_val[tra[tr2].eq];
          kshape = nmeq_val[tra[tr2].eqshape];
          L += gammaup(t-tstnm[tra[tr2].cl],kshape,kshape/mean);
          break;

        case WEI_TR:
          lam = nmeq_val[tra[tr2].eq];
          kshape = nmeq_val[tra[tr2].eqshape];
          L += weibullup(t-tstnm[tra[tr2].cl],lam,kshape);
          break;
      }
    }
    c = tra[tr].cf; if(c == NOTALIVE) break;
    cl = tra[tr].cl; cstnm[cl] = c; tstnm[cl] = t;
  }

  return L;
}

void Chain::nmchange(long tr, double dt, double dt2)   // The chnage in Linm from moving an evvent
{
  double mean, kshape, lam;

  switch(tra[tr].type){
    case GAMMA_TR:
      mean = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
      Linm += dgammaprob(dt,dt2,kshape,kshape/mean);
      break;

    case WEI_TR:
      lam = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
      Linm += dweibullprob(dt,dt2,lam,kshape);
      break;
    case FIXED_TR: emsg("Dist: EC2"); break;
  }
}

void Chain::nmcha(vector <EV> &ev, long e, long dir)  // The chnage in Linm from adding an evvent
{
  long ee, cl, n;
  long tr, trr;
  double t, dt, mean, kshape, lam;

  tr = ev[e].tr; t = ev[e].t;
  if(tra[tr].nm == 1){
    cl = tra[tr].cl;
    ee = e-1; while(ee > 0 && tra[ev[ee].tr].cl != cl) ee--;
    dt = t-ev[ee].t;

    switch(tra[tr].type){
      case GAMMA_TR:
        mean = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
        Linm += gammaprob(dt,kshape,kshape/mean)*dir;
        break;

      case WEI_TR:
        lam = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
        Linm += weibullprob(dt,lam,kshape)*dir;
        break;
        case FIXED_TR: emsg("Dist: EC3"); break;
    }
  }

  for(n = 0; n < tra[tr].ntraend; n++){
    trr = tra[tr].traend[n];
    cl = tra[trr].cl;
    ee = e-1; while(ee > 0 && tra[ev[ee].tr].cl != cl) ee--;
    dt = t-ev[ee].t;
    switch(tra[trr].type){
      case GAMMA_TR:
        mean = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
        Linm += gammaup(dt,kshape,kshape/mean)*dir;
        break;

      case WEI_TR:
        lam = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
        Linm += weibullup(dt,lam,kshape)*dir;
        break;
        case FIXED_TR: emsg("Dist: EC4"); break;
    }
  }
}

double gammaup(double x, double a, double b)     // Integral of the gamma function from x up to infinity
{
  int flag[1];
  double val;

  if(x < 0 || a < 0 || b < 0) return -large;
  val = log(1-gammds (x*b, a, flag ));
  if(flag[0] == 2) return -1000 - x*b;  // If in the extremes of distribution
  if(flag[0] != 0) emsg("Dist: EC6");

  return val;
}

double gammasamp(double a, double b)             // Draws a sample from x^(a-1)*exp(-b*x)
{
  if(a < 0 || b < 0) emsg("Dist: EC7");

  if(a < 1){
    double u = ran();
    return gammasamp(1.0 + a, b) * pow (u, 1.0 / a);
  }
  else{
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while(1 == 1){
      do{
        x = sqrt(-2*log(ran()))*cos(2*3.141592654*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      u = ran();

      if (u < 1 - 0.0331*x*x*x*x) break;

      if (log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
}

double gammasamptail(double a, double b, double tmin)     // Samples from the tail of a gamma distribution
{
  long loop = 0;
  double dt;

  if(a < 0 || b < 0 || tmin < 0) emsg("Dist: EC8");

  do{ dt = gammasamp(a,b); loop++;}while(loop < 30 && dt < tmin);
  if(loop == 30) return tmin - log(ran())/b;
  return dt;
}

double gammaprob(double x, double a, double b)              // Log probability of the gamma distribution
{
  if(x < 0 || a < 0 || b < 0) return -large;
  return (a-1)*log(x) - b*x + a*log(b) - lgamma(a);
}

double dgammaprob(double x, double x2, double a, double b)  // Difference in log probability from gamma
{
  return (a-1)*log(x2/x) - b*(x2-x);
}

double tgammaprob(double x, double a, double b)
{
  if(x < 0 || a < 0 || b < 0) return 0;
  return pow(x,(a-1))*exp(- b*x)*pow(a,b)/tgamma(a);
}

double weibullsamp(double lam, double kk)                 // Samples from the Weibull distribution
{
  if(lam < 0 || kk < 0) emsg("Dist: EC9");
  return lam*pow(-log(1-ran()),1.0/kk);
}

double weibullsamptail(double lam, double kk, double tmin)  // Samples from the tail of Weibull
{
  double val, sa;

  if(lam < 0 || kk < 0 || tmin < 0) emsg("Dist: EC10");

  val = exp(-pow(tmin/lam,kk));
  if(val < 0.0000000001) return tmin - log(ran())*(tlea-tent)*0.1;

  sa = lam*pow(-log(val*ran()),1.0/kk);
  return sa;
}

double weibullprob(double x, double lam, double kk)      // Log probability of sampling from Weibull
{
  if(x < 0 || lam < 0 || kk < 0) return -large;
  return log(kk/lam) + (kk-1)*log(x/lam) - pow(x/lam,kk);
}

double dweibullprob(double x, double x2, double lam, double kk) // Differnce in log probability
{
  return  (kk-1)*log(x2/x) - pow(x2/lam,kk) + pow(x/lam,kk);
}

double weibullup(double x, double lam, double kk)       // Log of integral of Weibull from x to infinity
{
  if(x < 0 || lam < 0 || kk < 0) return -large;
  return -pow(x/lam,kk);
}

double dweibullup(double x, double x2, double lam, double kk)  // Difference in the log of tail
{
  return -pow(x2/lam,kk)+pow(x/lam,kk);
}

double normalprob(double x, double mean, double var)    // Log of normal probability distribution
{
  if(var < 0) return -large;

  return -0.5*log(2*3.141592654*var) - (x-mean)*(x-mean)/(2*var);
}

double lognormalprob(double x, double mean, double var) // Log of lognormal probability distribution
{
  if(x < 0 || var < 0) return -large;
  double lx = log(x);
  return -0.5*log(2*3.141592654*var) - lx - (lx-mean)*(lx-mean)/(2*var);
}

double expprob(double x, double rate)                   // Log of exponential probability distribution
{
  if(rate < 0) return -large;

  return log(rate) - rate*x;
}

double betaprob(double x, double a, double b)           // Log of beta probability distribution
{
  if(x < 0 || x > 1 || a < 0 || b < 0) return -large;
 
  return (a-1)*log(x) + (b-1)*log(1-x) + lgamma(a+b) - lgamma(a) - lgamma(b);
}

double betasamp(double a, double b)                     // Sample from the beta probability distribution
{
  double X, Y;

  if(a < 0 || b < 0) emsg("Dist: EC11");
  X = gammasamp(a,1);
  Y = gammasamp(b,1);
  return X/(X+Y);
}

void Chain::setnmeq()                      // Calculates equations for likelihood of non-Markovian events
{
  long k, eq;
  for(k = 0; k < nnmeq; k++){ eq = nmeq[k]; nmeq_val[eq] = calculatenotdep(eq,param);}
}
