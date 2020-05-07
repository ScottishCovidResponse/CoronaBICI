// This file considers proposals which simulate event sequences for individuals

void Chain::indsim_prop(long i)            // Proposal which resimulates the trajectory for an individual
{
  long c, a, e;
  double al, probif, probfi, Ltoti, Ltotf, pr, fac, t;

  partcl = clall;

  if(i < nind) ntr_indsim[i]++; else ntr_indsimuo++;

  getlifespan(i);
  indsim_init(i,indtbirth[i],0);

  c = selectstart(i,tent);
  probif = probstart(i,tent,c); 
	if(probif == -large){ if(i < nind) nfa_indsim[i]++; else nfa_indsimuo++; return;}
  c += ctop;

  pr = indsim(i,c); if(pr == -large){ if(i < nind) nfa_indsim[i]++; else nfa_indsimuo++; return;}
  if(checkon == 1){ double dd = pr-probindsim(evnew); if(dd*dd > tiny) emsg("Indsim: EC1");}
  probif += pr;

  if(checkevdt() == 1) return;

  probfi = probstart(i,indev[i][0].t,tra[indev[i][0].tr].cf%ncomps);
  probfi += probindsim(indev[i]);

  Ltoti = Lwnm();
  indcha(i);
  Ltotf = Lwnm();

  al = exp(Ltotf - Ltoti + probfi - probif);

  if(ran() < al){
    if(i < nind) nac_indsim[i]++; else nac_indsimuo++;
    fac = 1;
  }
  else{
    indrev(i);
    fac = 0.2;
  }

  if(samp < burnin){
    c = tra[indev[i][0].tr].cf; if(c == NOTALIVE) emsg("Indsim: EC2");
    t = indev[i][0].t;
    if(i < nind){
      if(t == 0) sim_initc[i][c%ncomps] += fac; else sim_enterc[i][c%ncomps] += fac;
    }
    else{
      if(t == 0) a = naddremdiv; else a = long(naddremdiv*t/tmax);
      if(t == 0) simuo_initc[a][c%ncomps] += fac; else simuo_enterc[a][c%ncomps] +=fac;
    }
  }
}

double Chain::indsimulate(long i)                        // Simulates event sequence for an individual
{
  long c, e;
  double prob;

  indsim_init(i,indtbirth[i],0);

  c = selectstart(i,tent);
  prob = probstart(i,tent,c); if(prob == -large) return -large;
  c += ctop;

  prob += indsim(i,c);

  for(e = 0; e < long(evnew.size())-1; e++){ if(evnew[e+1].t-evnew[e].t < evdtmin) return -large;}

  return prob;
}

double Chain::probindsimulate(long i)              // Probability of event sequence for an individual
{
  long c;
  double prob;

  indsim_init(i,indtbirth[i],0);

  prob = probstart(i,indev[i][0].t,tra[indev[i][0].tr].cf%ncomps);
  prob += probindsim(indev[i]);

  return prob;
}

double Chain::indsim(long i, long c)                  // Simulates a new sequence for an individual
{
  long e, ee, k, co, j, ch, cf, p, m, n, fl, cl, fev;
  double t, tt, tst, tnext, tnextold, z, R, r, val, kshape, mean, lam, futnext, tmid, tsta, prob = 0;
  long tr;

  t = evp[0].t;
  evnew.clear();
  if(nmfl == 1){ futev.clear(); nfuteref = 0; futeref.clear();}

  if(t == 0) tr = trabeg+c; 
  else{ 
    tr = compiftra[NOTALIVE][c]; if(tr < 0) return -large;
    if(i < nind && indlikeenter[i] == 0) tr += moventra;
  }
  addevnew(tr,t);

  co = NOTALIVE; k = 0; if(k < nexist) tnextold = texist[k]; else tnextold = tmax;
  for(e = 0; e < nevp-1; e++){
    switch(evp[e].ty){
      case CHC:
        tr = compiftra[c][c + evp[e].ch]; if(tr < 0) return -large;
        c = addevnew(tr,t);
        break;

      case FIXEV:
        fev = evp[e].ch; cl = fixev[fev].cl; if(compval[c][cl] != fixev[fev].trai) return -large;
        tr = compiftra[c][c + fixev[fev].dc]; if(tr < 0) return -large;
        if(fixev[fev].like == 0) tr += moventra;
        c = addevnew(tr,t);
        break;
      default:
        emsg("Invalid default on " LINE_STRING " in " __FILE__);
    }

    tnext = evp[e+1].t;
    if(t < tnext){
      do{
        if(compleavesimcldep[partcl][c].size() > 0){
          tt = t+partdt; if(tt > tnext) tt = tnext;
          tmid = 0.5*(t+tt); 
          while(tnextold < tmid){
            co = cexist[k]; k++; if(k > nexist) emsg("Indsim: EC3");
            if(k < nexist) tnextold = texist[k]; else tnextold = tmax;
          }
          getdep(compleavesimcldep[partcl][c],co,tmid);
        }
        else{ Rdep = 0; tt = tnext;}

        R = Rdep + Rnotdep;

        if(nmfl == 1){ if(nfuteref > 0){ futnext = futev[futeref[nfuteref-1]].t; if(tt > futnext) tt = futnext;} else futnext = large;} else futnext = large;

        if(R == 0) t = tt;
        else{
          tst = t;
          t += -log(ran())/R;
          if(t > tt) t = tt;
          else{ // Actual event
            z = ran()*R;
            if(z <= Rdep){
              j = 0; while(j < nRdep && z > Rdepst[j]) j++; if(j == nRdep) emsg("Indsim: EC4");
              r = rdep[j]; tr = compleavesimcldep[partcl][c][j];
            }
            else{
              z -= Rdep;
              j = 0; while(j < nRnotdep && z > Rnotdepst[j]) j++; if(j == nRnotdep) emsg("Indsim: EC5");
              r = rnotdep[j]; tr = compleavesimclnotdep[partcl][c][j];
            }
            prob += log(r);
            c = addevnew(tr,t);
          }
          prob -= (t-tst)*R;
        }

        if(t == futnext){   // Future event
          n = futeref[nfuteref-1]; //futeref.pop_back(); nfuteref--;
          tr = compiftra[c][c + tra[futev[n].tr].dc]; if(tr < 0) return -large;
          c = addevnew(tr,t);
        }
      }while(t < tnext);
    }
  }

  if(tlea == tmax) tr = traend+c; 
  else{
    tr = compiftra[c][NOTALIVE]; if(tr < 0) return -large;
    if(i < nind && indlikeleave[i] == 0) tr += moventra;
  }
  addevnew(tr,tlea);

  return prob;
}

long Chain::addevnew(long tr, double t)           // Adds a new event and considers future events
{
  long n, c, cl, fref, m;
  long trr;
  double tt, val, kshape;

  if(tr < 0) emsg("Indsim: EC6");

  if(nmfl == 1){   // removes relevant future events
    cl = tra[tr].cl;
    for(n = nfuteref-1; n >= 0; n--){
      if(tra[futev[futeref[n]].tr].cl == cl || cl == -1){ 
        if(n == nfuteref-1) futeref.pop_back();
        else futeref.erase(futeref.begin()+n);
        nfuteref--;
      }
    }
  }

  EV ev; ev.tr = tr; ev.t = t; evnew.push_back(ev);
  c = tra[tr].cf;
  getnotdep(compleavesimclnotdep[partcl][c]);

  if(nmfl == 1){   // Adds in new future events
    for(n = 0; n < tra[tr].nnonexptra; n++){
      trr = tra[tr].nonexptra[n];
      if((partcl == clall || tra[trr].cl == partcl) && !(tra[trr].capev >= 0 && capevall[tra[trr].capev] == 1)){
        val = nmeq_val[tra[trr].eq];
        switch(tra[trr].type){
          case GAMMA_TR: kshape = nmeq_val[tra[trr].eqshape]; tt = t+gammasamp(kshape,kshape/val); break;
          case WEI_TR: kshape = nmeq_val[tra[trr].eqshape]; tt = t+weibullsamp(val,kshape); break;
          case FIXED_TR: tt = t+val; break;
          default:
            emsg("Invalid default on " LINE_STRING " in " __FILE__);
        }

        m = 0; while(m < nfuteref && tt < futev[futeref[m]].t) m++;
        if(m == nfuteref) futeref.push_back(futev.size());
        else futeref.insert(futeref.begin()+m,futev.size());
        nfuteref++;
        FUTEV fev; fev.tr = trr; fev.tnow = t; fev.t = tt; futev.push_back(fev);
      }
    }
  }

  return c;
}

double Chain::probindsim(vector<EV> &ev)     // Probability of simulating a new sequence for an individual
{
  long e, k, co, j, ch, cf, p, m, n, fl, c, ee, nev, fev, cl;
  double t, tt, tst, tnext, tnextold, z, R, r, val, tmid, prob = 0, kshape, mean, lam;
  long tr;

  c = tra[ev[0].tr].cf; ee = 1; nev = ev.size();
  getnotdep(compleavesimclnotdep[partcl][c]); 

  co = NOTALIVE; k = 0; if(k < nexist) tnextold = texist[k]; else tnextold = tmax;

  t = evp[0].t;
  for(e = 0; e < nevp-1; e++){
    switch(evp[e].ty){
      case CHC: case OTHCL:
        cf = c + evp[e].ch;
        tr = compiftra[c][cf];
        if(ee == nev || ev[ee].tr != tr) emsg("Indsim: EC7");
        c = cf; ee++;
        getnotdep(compleavesimclnotdep[partcl][c]);
        break;

      case FIXEV:
        fev = evp[e].ch; cl = fixev[fev].cl;
        if(compval[c][cl] != fixev[fev].trai) emsg("Indsim: EC8");
        cf = c + fixev[fev].dc;
        tr = compiftra[c][cf]; if(fixev[fev].like == 0) tr += moventra;
        if(ee == nev || ev[ee].tr != tr) emsg("Indsim: EC9");
        if(fixev[fev].like == 0) tr += moventra;
        c = cf; ee++;
        getnotdep(compleavesimclnotdep[partcl][c]); 
        break;
      default:
        emsg("Invalid default on " LINE_STRING " in " __FILE__);
    }

    tnext = evp[e+1].t;
    if(t < tnext){
      do{
        if(compleavesimcldep[partcl][c].size() > 0){
          tt = t+partdt; if(tt > tnext) tt = tnext;
          tmid = 0.5*(t+tt); while(tnextold < tmid){ co = cexist[k]; k++; if(k < nexist) tnextold = texist[k]; else tnextold = tmax;}
          getdep(compleavesimcldep[partcl][c],co,tmid);
        }
        else{ Rdep = 0; tt = tnext;}

        R = Rdep + Rnotdep;

        tst = t;
        if(ee < nev && ev[ee].t < tt){ // Actual event occurs
          tr = ev[ee].tr; t = ev[ee].t; ee++;
          switch(tra[tr].type){
            case EXP_TR:
              if(transdep[tra[tr].eq] == 1){
                j = 0; while(j < nRdep && compleavesimcldep[partcl][c][j] != tr) j++; 
								if(j == nRdep) emsg("Indsim: EC10");
                r = rdep[j];
              }
              else{
                j = 0; while(j < nRnotdep && compleavesimclnotdep[partcl][c][j] != tr) j++; 
								if(j == nRnotdep) emsg("Indsim: EC11");
                r = rnotdep[j];
              }
              prob += log(r);
              c = tra[tr].cf;
              break;

            default:
              c = tra[tr].cf;
              break;
          }
          getnotdep(compleavesimclnotdep[partcl][c]);
        }
        else t = tt;
        prob -= (t-tst)*R;
      }while(t < tnext);
    }
  }

  if(tra[ev[ee].tr].cf != NOTALIVE) emsg("Indsim: EC12");

  return prob;
}

void Chain::indsim_init(long i, double tbirth, long incobs)   // Initialises division used in simualation
{
  long e, ag, ti, ob, f, tr, fi;
  double dt, t;

  ist = i;

  partdt = tmax/49.99999; dt = (tlea-tent)/9.9999; if(dt < partdt) partdt = dt;

  cexist.clear(); texist.clear();
  for(e = 0; e < nindev[i]; e++){ 
		cexist.push_back(tra[indev[i][e].tr].cf); texist.push_back(indev[i][e].t);
	}
  nexist = cexist.size();

  evp.clear(); nevp = 0;
  evp.push_back(EVPART ()); evp[nevp].ty = BEG; evp[nevp].t = tent; nevp++;

  ctop = 0;
  if(tbirthfl == 1){
    ag = 0; while(ag < nage && tbirth + age[ag] < tent) ag++;
    ctop += ag*classmult[agecl];
    while(ag < nage && tbirth + age[ag] < tlea){ 
			evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[agecl]; 
			evp[nevp].t = tbirth + age[ag]; nevp++; ag++;
		}
  }

  ti = 0; while(ti < nsettime && settime[ti] < tent) ti++;
  ctop += ti*classmult[settimecl];
  while(ti < nsettime && settime[ti] < tlea){ 
		evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[settimecl]; 
		evp[nevp].t = settime[ti]; nevp++; ti++;
	}

  if(i < nind){
    for(e = 0; e < nindfixev[i]; e++){ 
			evp.push_back(EVPART ()); evp[nevp].ty = FIXEV; evp[nevp].ch = indfixev[i][e]; 
			evp[nevp].t = fixev[indfixev[i][e]].t; nevp++;
		}
    if(incobs == 1){
      for(ob = 0; ob < nindobs[i]; ob++){ 
				evp.push_back(EVPART ()); evp[nevp].ty = OBS; evp[nevp].ch = indobs[i][ob]; 
				evp[nevp].t = indobst[i][ob]; nevp++;
			}
    }
  }

  if(partcl < clall){  // Transitions in other classifications
    if(i >= nind) emsg("Indsim: EC13");
    fi = 0;
    for(e = 1; e < nindev[i]-1; e++){
      tr = indev[i][e].tr; t = indev[i][e].t;
      if(tra[tr].cl != partcl && tra[tr].cl < clall){
        while(fi < nindfixev[i] && fixev[indfixev[i][fi]].t < t) fi++;
        if(!(fi < nindfixev[i] && fixev[indfixev[i][fi]].t == t)){
           evp.push_back(EVPART ()); evp[nevp].ty = OTHCL; evp[nevp].ch = tra[tr].dc; 
					 evp[nevp].t = t; nevp++;
        }
      }
    }
  }

  f = 0; while(f < ncapevtrange && capevtrange[f] <= tent) f++;
  ctop += f*classmult[capevcl];
  while(f < ncapevtrange && capevtrange[f] < tlea){ 
    evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[capevcl];
		evp[nevp].t = capevtrange[f]; 
    nevp++; f++;
  }

  evp.push_back(EVPART ()); evp[nevp].ty = END; evp[nevp].t = tlea;
  evp[nevp].ch = 0; if(tlea < tmax && !(i < nind && indlikeleave[i] == 0)) evp[nevp].ch = 1;
  nevp++;

  sort(evp.begin(),evp.end(),evpcomp);
}

long Chain::selectstart(long i, double t)          // Selects the starting compartment for a simulation
{
  long c, a;
  double z;

  if(i < nind){
    if(partcl == clall){
      if(t == 0){
				z = ran()*sim_initcsumtot[i]; 
				c = 0; while(c < ncomps && z > sim_initcsum[i][c]) c++; 
				if(c == ncomps) emsg("Indsim: EC14");
			}
      else{ 
				z = ran()*sim_entercsumtot[i]; 
				c = 0; while(c < ncomps && z > sim_entercsum[i][c]) c++; 
				if(c == ncomps) emsg("Indsim: EC14");
			}
    }
    else{
      long cnow = tra[indev[i][0].tr].cf%ncomps;
      if(t == 0){ 
				z = ran()*simcl_initcsumtot[i][partcl]; 
				c = 0; while(c < nclassval[partcl] && z > simcl_initcsum[i][partcl][c]) c++; 
				if(c == nclassval[partcl]) emsg("Indsim: EC15");
			}
      else{ 
				z = ran()*simcl_entercsumtot[i][partcl]; 
				c = 0; while(c < nclassval[partcl] && z > simcl_entercsum[i][partcl][c]) c++; 
				if(c == nclassval[partcl]) emsg("Indsim: EC16");
			}
      c = cnow + (c-compval[cnow][partcl])*classmult[partcl];
    }
  }
  else{
    if(t == 0) a = naddremdiv; else a = long(naddremdiv*t/tmax);
    if(t == 0){ z = ran()*simuo_initcsumtot[a]; 
		c = 0; while(c < ncomps && z > simuo_initcsum[a][c]) c++; 
		if(c == ncomps) emsg("Indsim: EC17");}
    else{
			z = ran()*simuo_entercsumtot[a];
			c = 0; while(c < ncomps && z > simuo_entercsum[a][c]) c++; 
			if(c == ncomps) emsg("Indsim: EC18");
		}
  }

  return c;
}

double Chain::probstart(long i, double t, long c)   // Probability of starting compartment in a simulation
{
  long a;
  if(i < nind){
    if(partcl == clall){
      if(t == 0) return logsim_initc[i][c]; else return logsim_enterc[i][c];
    }
    else{
      if(t == 0) return logsimcl_initc[i][partcl][compval[c][partcl]]; 
      else return logsimcl_enterc[i][partcl][compval[c][partcl]];
    }
  }
  else{
    if(t == 0) a = naddremdiv; else a = long(naddremdiv*t/tmax);
    if(t == 0) return logsimuo_initc[a][c]; else return logsimuo_enterc[a][c];
  }
}

void Chain::getdep(vector <long> &vec, long co, double t)  // The total rate of dependent transitions
{
  long j, eq;
  Rdep = 0; nRdep = vec.size();
  for(j = 0; j < nRdep; j++){
    eq = tra[vec[j]].eq; if(eq < 0) emsg("Indsim: EC19");
 
    rdep[j] = transradep(eq,t,co);
    Rdep += rdep[j];
    Rdepst[j] = Rdep;
  }
}

void Chain::getnotdep(vector <long> &vec)                  // The total rate of not dependent transitions
{
  long j, eq;

  Rnotdep = 0; nRnotdep = vec.size(); 
  for(j = 0; j < nRnotdep; j++){
    if(vec[j] < 0 || vec[j] >= ntra) emsg("Indsim: EC20");
    eq = tra[vec[j]].eq; if(eq < 0) emsg("Indsim: EC21");
    rnotdep[j] = transnotdepeq_val[transdepref[eq]];
    Rnotdep += rnotdep[j];
    Rnotdepst[j] = Rnotdep;
  }
}

long Chain::getlifespan(long i)                    // Gets the time over which an individual is alive
{
  long e;
  tent = indev[i][0].t; e = nindev[i]-1; 
  while(e > 0 && tra[indev[i][e].tr].cf != NOTALIVE) e--; tlea = indev[i][e].t;
  if(tlea > tmax || e == 0) emsg("Indsim: EC22");
  return e;
}
