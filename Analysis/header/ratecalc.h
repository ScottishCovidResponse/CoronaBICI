 // This file calculates the rates of transitions
 
double ratecalc(long eq, double *popnum, double *param)   // Calculates the rate of an event
{
  double r;

  r = calculate(eq,popnum,param);
  if(r < ratesmall) r = ratesmall;
  return r;
}

double ratecalcdep(long d, double *popnum, double *param) // Calculates the rate of a dep. event
{
  long i, eq;
  double num, num1, num2, r;

  eq = transdepeq[d];
  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
    case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
    case POPNUM: num1 = popnum[eqncalcfastnum1[d][i]]; break;
    case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
    case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
    case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
    case POPNUM: num2 = popnum[eqncalcfastnum2[d][i]]; break;
    case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
    case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
    }

    switch(eqncalcop[eq][i]){
    case ADD: num = num1+num2; break;
    case TAKE: num = num1-num2; break;
    case MULTIPLY: num = num1*num2; break;
    case DIVIDE: if(num2 == 0) emsg("Ratecalc: EC1 A division by zero!"); num = num1/num2; break;
    case POWER: num = pow(num1,num2); break;
    case EXPFUNC: num = exp(num2); break;
    case SINFUNC: num = sin(num2); break;
    case COSFUNC: num = cos(num2); break;
    case LOGFUNC: if(num2 <= 0) emsg("Ratecalc: EC2 Log cannot be calculated"); num = log(num2); break;
    case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
  case PARAM: r = param[eqncalcansnum[eq]]; break;
  case POPNUM: r = popnum[eqncalcansfastnum[d]]; break;
  case REG: r = regcalc[eqncalcansnum[eq]]; break;
  case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

double ratecalcnotdep(long eq, double *paramval)    // Calculates the rate of a non dep. event
{
  long i;
  double num, num1, num2, r;

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
    case PARAM: num1 = paramval[eqncalcnum1[eq][i]]; break;
    case POPNUM: emsg("Ratecalc: EC3"); break;
    case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
    case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
    case PARAM: num2 = paramval[eqncalcnum2[eq][i]]; break;
    case POPNUM: emsg("Ratecalc: EC4"); break;
    case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
    case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
    }

    switch(eqncalcop[eq][i]){
    case ADD: num = num1+num2; break;
    case TAKE: num = num1-num2; break;
    case MULTIPLY: num = num1*num2; break;
    case DIVIDE: if(num2 == 0) emsg("Ratecalc: EC5 A division by zero!"); num = num1/num2; break;
    case POWER: num = pow(num1,num2); break;
    case EXPFUNC: num = exp(num2); break;
    case SINFUNC: num = sin(num2); break;
    case COSFUNC: num = cos(num2); break;
    case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
    case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
  case PARAM: r = paramval[eqncalcansnum[eq]]; break;
  case POPNUM: emsg("Ratecalc: EC6"); break;
  case REG: r = regcalc[eqncalcansnum[eq]]; break;
  case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

// Removes the individual at c when the quantity is calculated
double ratecalcdeptakeoff(long d, double *popnum, double *param, long c)    
{
  long i, eq, p;
  double num, num1, num2, r;

  eq = transdepeq[d];

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
      case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
      case POPNUM:
        p = eqncalcfastnum1[d][i];
        num1 = popnum[p]; if(c != -1) num1 += popnumtake[eq_popnum[eq][p]][c];
        break;
        case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
        case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
      default:
        switch(eqncalcop[eq][i]){
          case ADD: case TAKE: case MULTIPLY: case DIVIDE: case POWER: emsg("Ratecalc: EC7"); break;
        }
        break;
    }

    switch(eqncalc2[eq][i]){
      case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
      case POPNUM: p = eqncalcfastnum2[d][i]; num2 = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
      case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
      case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
      default: emsg("Ratecalc: EC8"); break;
    }

    switch(eqncalcop[eq][i]){
      case ADD: num = num1+num2; break;
      case TAKE: num = num1-num2; break;
      case MULTIPLY: num = num1*num2; break;
      case DIVIDE: if(num2 == 0) emsg("Ratecalc: EC9 A division by zero!"); num = num1/num2; break;
      case POWER: num = pow(num1,num2); break;
      case EXPFUNC: num = exp(num2); break;
      case SINFUNC: num = sin(num2); break;
      case COSFUNC: num = cos(num2); break;
      case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
      case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
      default: emsg("probratecalc1"); break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
    case PARAM: r = param[eqncalcansnum[eq]]; break;
    case POPNUM: p = eqncalcansfastnum[d]; r = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
    case REG: r = regcalc[eqncalcansnum[eq]]; break;
    case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
    default: emsg("Ratecalc: EC10"); break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

double Chain::transra(long eq, double t, long cnow)   // Transition rate but removes individual at cnow
{
  if(transdep[eq] == 1) return transradep(eq,t,cnow);
  else return transnotdepeq_val[transdepref[eq]];
}

double Chain::transradep(long eq, double t, long cnow) // Transition rate but removes individual at cnow
{
  long d, e, ee, kd;

  if(eq == -1 || transdep[eq] != 1) emsg("Ratecalc: EC11");                     

  d = transdepref[eq];
  t *= onefac;

	if(discon == 1){
		kd = d*DX + long(t*dfac);
		if(transdepeqrecalc[d][cnow] == 0) return depeq_disc_evval[kd];
		else{ return ratecalcdeptakeoff(d,depeq_disc_evpopnum[kd],param,cnow);}
	}
	else{	
		e = depeqdiv[d][long(t*depeqdivmax/tmax)];
		do{
			ee = depeq_evnext[e]; if(ee == -1) emsg("Ratecalc: EC12");
			if(ee == -1 || depeq_evt[ee] > t) break;
			e = ee;
		}while(1 == 1);

		if(transdepeqrecalc[d][cnow] == 0) return depeq_evval[e];
		else{ return ratecalcdeptakeoff(d,depeq_evpopnum[e],param,cnow);}
	}
}

double Chain::transrate(long eq, double t)     // Calculates the transition rate 
{
	long d;
  long e, ee;

	if(eq == -1 || transdep[eq] != 1) emsg("Ratecalc: EC13");   
  d = transdepref[eq];
	
	if(discon == 1){
		return depeq_disc_evval[d*DX + long(t*dfac)];
	}
	else{	
		t *= onefac;
		  
		e = depeqdiv[d][long(t*depeqdivmax/tmax)];
		do{
			ee = depeq_evnext[e]; if(ee == -1) emsg("Ratecalc: EC14");
			if(ee == -1 || depeq_evt[ee] > t) break;
			e = ee;
		}while(1 == 1);
		return depeq_evval[e];
	}
}
