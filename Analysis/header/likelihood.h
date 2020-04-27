// This file calculates changes to likelihood under a change to an individual

double Chain::L()               // The posterior probability
{ 
	return logsum[nindtot] + Lpri + (Liinit + Lir - Liexp + Linm)*invTLi 
			+ (Lev_pd + Lob_pd)*invT_pd + Lob_st*invT_st + Lpop*invT_pop + Lder*invT_der;
}

double Chain::Lwnm()           // The posterior probability expliding non-Markovian transitions
{ 
	return logsum[nindtot] + Lpri + (Liinit + Lir - Liexp)*invTLi 
			+ (Lev_pd + Lob_pd)*invT_pd + Lob_st*invT_st + Lpop*invT_pop + Lder*invT_der;
}

void Chain::indcha(long i)     // Makes a change to individual i
{
  indchange(i,0,mlarge,large);
}

void Chain::indrev(long i)     // Reverses the change to individual i
{
  evnew = evrev;
  indchange(i,0,mlarge,large);
}

// Makes a change to the event sequence on individual i
void Chain::indchange(long i, long estart, double t, double tend) 
{
  long e, tr, trold, trnew, ci, cf, nevold, nevnew, eold, enew, cold, cnew, k, ob, o, co, nob, whole;
  double tt, told, tnew, tob;
  vector <long> remlist, addlist;

	//if(i == 3195){  cout << "\n\n"; oe("st",indev[i]); oe("new",evnew);}
  
  timeprop[INDCHANGE] -= clock();

  if(t == mlarge) whole = 1; else whole = 0;

  if(checkon == 1) checkevseqsimp(evnew);

  nevold = nindev[i];
  nevnew = evnew.size();

  if(nmfl == 1){
    if(whole == 1) Linm -= likenm(i);
    else{ 
			e = estart; 
			while(e < nevold && indev[i][e].t < tend){ nmcha(indev[i],e,-1); e++;} nmcha(indev[i],e,-1);
		}
  }

  if(estart == 0) cold = NOTALIVE;
  else cold = tra[indev[i][estart].tr].ci;
  cnew = cold;

  eold = estart; if(eold < nevold) told = indev[i][eold].t; else told = large;
  if(eold >= nevold) emsg("Likelihood: EC1");
  enew = 0; if(enew < nevnew) tnew = evnew[enew].t; else tnew = large;

  if(estart == 0){
    if(told == 0 && indev[i][0].tr != tranull) ci = tra[indev[i][0].tr].cf%ncompswa; else ci = NOTALIVE;  
		if(tnew == 0 && evnew[0].tr != tranull) cf = tra[evnew[0].tr].cf%ncompswa; else cf = NOTALIVE;

    if(ci != cf){
      if(ci != NOTALIVE){ 
				Liinit -= log(probinit[ci]); ncompinit[ci]--; 
				Liinit -= logsum[nindinit]; nindinit--; Liinit += logsum[nindinit];
			}
      if(cf != NOTALIVE){ 
				Liinit += log(probinit[cf]); ncompinit[cf]++; 
				Liinit -= logsum[nindinit]; nindinit++; Liinit += logsum[nindinit];
			}
    }
  }

  if(i < nind){ 
		nob = nindobs[i]; 
		o = 0; while(o < nob && indobst[i][o] < t) o++; 
		if(o  == nob) tob = large; 
		else tob = indobst[i][o];
	}
	else tob = large;
	
  do{
	  if(told < tnew) tt = told; else tt = tnew;
    if(t != tt && cold != cnew) secchange(i,cold,cnew,t,tt);
	
    t = tt;

    while(tob < tt){       // Takes into account observations
      if(cold != cnew){
        ob = indobs[i][o];
        co = obsprob_cope[ob][cold]; 
				if(co < 0) Lob_st -= notobsdL; else{ Lob_st -= cope_val[co]; cope_num[co]--;}
				
        co = obsprob_cope[ob][cnew]; 
				if(co < 0) Lob_st += notobsdL; else{ Lob_st += cope_val[co]; cope_num[co]++;}
      }
      o++; if(o == nindobs[i]) tob = large; else tob = indobst[i][o];
    }

    if(told < tnew){
      remlist.push_back(eold); 
      trold = indev[i][eold].tr;
      cold = tra[trold].cf; eold++; 
      if(eold < nevold) told = indev[i][eold].t; else told = large;
    }
    else{
      if(tnew < told){
        addlist.push_back(enew); trnew = evnew[enew].tr; cnew = tra[trnew].cf; enew++; 
				if(enew < nevnew) tnew = evnew[enew].t; else tnew = large;
      }
      else{             // simultaneous events
        trold = indev[i][eold].tr; trnew = evnew[enew].tr;
        if(trold != trnew){ remlist.push_back(eold); addlist.push_back(enew);}
       
        cold = tra[trold].cf; eold++; if(eold < nevold) told = indev[i][eold].t; else told = large;
        cnew = tra[trnew].cf; enew++; if(enew < nevnew) tnew = evnew[enew].t; else tnew = large;
      }
    }
  }while(enew < nevnew || (eold < nevold && told <= tend));

  if(whole == 1) evrev = indev[i]; 
	else{ evrev.clear(); for(e = estart; e < eold; e++) evrev.push_back(indev[i][e]);}

  for(k = long(remlist.size())-1; k >= 0; k--) indremevent(i,remlist[k]); // Removes exisiting events

  e = estart;                                                             // Adds in new events
  for(k = 0; k < addlist.size(); k++){
    t = evnew[addlist[k]].t;
    while(e < nindev[i] && indev[i][e].t < t) e++;
    indaddevent(i,e,evnew[addlist[k]]);
  }

  if(nmfl == 1){
    if(whole == 1) Linm += likenm(i);
    else{ 
			e = estart; nevold = nindev[i]; 
			while(e < nevold && indev[i][e].t < tend){ nmcha(indev[i],e,1); e++;} nmcha(indev[i],e,1);
		}
  }

  if(checkon == 1) checkevseq(indev[i]);

  timeprop[INDCHANGE] += clock();
}

void Chain::indremevent(long i, long k)                   // Removes an event for an individual
{
  long eq, d, tr, j, ce, fi, kd;
  double t, tt;
  long e, ee;

  tr = indev[i][k].tr;

  if(capevfl == 1){                // Change coming from captured event
    ce = tra[tr].capev;
    if(ce >= 0){
      tt = indev[i][k].t;
      if(i < nind){
        fi = 0; while(fi < nindfixev[i] && fixev[indfixev[i][fi]].t < tt) fi++;
        if((fi < nindfixev[i] && fixev[indfixev[i][fi]].t == tt && fixev[indfixev[i][fi]].capev == ce)
						|| tt == indfixtenter[i] || tt == indfixtleave[i]){ 
          capev_num[ce]--; Lev_pd -= capev_val[ce];
        }
        else{ capev_oneminus[ce]--; Lev_pd -= capev_valoneminus[ce];}
      }
      else{ capev_oneminus[ce]--; Lev_pd -= capev_valoneminus[ce];}
    }
  }

  if(tra[tr].type == EXP_TR && tra[tr].like == 1){
    eq = tra[tr].eq; d = transdepref[eq];
    if(transdep[eq] == 1){
			t = indev[i][k].t;
			if(discon == 1){
				kd = d*DX + long(t*dfac);
				depeq_disc_evn[kd]--;
				Lir -= log(depeq_disc_evval[kd]);
			}
			else{
				e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];
				do{
					ee = depeq_evnext[e]; if(ee == -1) emsg("Likelihood: EC2");
					tt = depeq_evt[ee];
					if(tt >= t) break;
					e = ee;
				}while(1 == 1);

				j = 0; while(j < ndepeq_ev_ev[e] && t != depeq_ev_evt[e][j]) j++;
				if(j == ndepeq_ev_ev[e]) emsg("Likelihood: EC3");

				ndepeq_ev_ev[e]--;
				if(j < ndepeq_ev_ev[e]) depeq_ev_evt[e][j] = depeq_ev_evt[e][ndepeq_ev_ev[e]];
				depeq_ev_evt[e].pop_back();

				Lir -= log(depeq_evval[e]);
			}
    }
    else{
      Lir -= log(transnotdepeq_val[d]); transnotdepeq_num[d]--;
    }
  }

  indev[i].erase(indev[i].begin()+k);
  nindev[i]--;
}

void Chain::indaddevent(long i, long k, EV ev)     // Adds an event for an indivual
{
  long eq, d, tr, ce, fi, kd;
  double t, tt;
  long e, ee;

  tr = ev.tr;

  if(capevfl == 1){                // Change coming from captured event
    ce = tra[tr].capev;
    if(ce >= 0){
      tt = ev.t;
      if(i < nind){
        fi = 0; while(fi < nindfixev[i] && fixev[indfixev[i][fi]].t < tt) fi++;
        if((fi < nindfixev[i] && fixev[indfixev[i][fi]].t == tt && fixev[indfixev[i][fi]].capev == ce) 
						|| tt == indfixtenter[i] || tt == indfixtleave[i]){ 
          capev_num[ce]++; Lev_pd += capev_val[ce];
        }
        else{ capev_oneminus[ce]++; Lev_pd += capev_valoneminus[ce];}
      }
      else{ capev_oneminus[ce]++; Lev_pd += capev_valoneminus[ce];}
    }
  }

  if(tra[tr].type == EXP_TR && tra[tr].like == 1){
    eq = tra[tr].eq; d = transdepref[eq];
    if(transdep[eq] == 1){
			t = ev.t;
			if(discon == 1){
				kd = d*DX + long(t*dfac);
				depeq_disc_evn[kd]++;
				Lir += log(depeq_disc_evval[kd]);
			}
			else{
				e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];
				do{
					ee = depeq_evnext[e]; if(ee == -1) emsg("Likelihood: EC4");
					tt = depeq_evt[ee];
					if(tt >= t) break;
					e = ee;
				}while(1 == 1);
				depeq_ev_evt[e].push_back(t);
				ndepeq_ev_ev[e]++;
				Lir += log(depeq_evval[e]);
			}
    }
    else{
      Lir += log(transnotdepeq_val[d]); transnotdepeq_num[d]++;
    }
  }

  indev[i].insert(indev[i].begin()+k,ev);
  nindev[i]++;
}

// Changes the state of an individual for a particular time span
void Chain::secchange(long i, long ci, long cf, double ti, double tf)         
{
  long k, kmax, j, jmax, p, co, cap, cpe, ob, div, dm, d, ref;
  double tt, dt;
  EQCH eqch;
  NDEQCH ndeqch;

  //if(plfl == 1) cout << ci << " " << cf << " " << ti << " "<< tf <<  "sec cha\n";

  if(capfl == 1 || popfl == 1 || derfl == 1){
    div = long(ncompcapdiv*ti/tmax);         // Calculates how pd chages for a capture campaign

    if(popfl == 1) poplist.clear();
    if(derfl == 1) derlist.clear();

    if(ci != NOTALIVE){
      j = compcapdiv[ci][div]; tt = compcap[ci][j].t; while(tt < ti){ j++; tt = compcap[ci][j].t;}

      while(tt < tf){
        switch(compcap[ci][j].ty){
          case CAP:
            cap = compcap[ci][j].val;
            cpe = capprob_cpe[cap][ci];
            if(i < nind) ob = capindob[cap][i]; else ob = -1;
            if(ob >= 0){   // Observed
              if(i >= nind) emsg("Likelihood: EC4");
              Lob_pd -= cpe_val[cpe]; cpe_num[cpe]--;
            }
            else{
              if(cpe < 0) emsg("Likelihood: EC5");
              Lob_pd -= cpe_valoneminus[cpe]; cpe_oneminus[cpe]--;
            }
            break;

          case POP:
            p = compcap[ci][j].val; if(popmap[p] == 0) poplist.push_back(p); popmap[p]--; 
            break;

          case DER:
            dm = compcap[ci][j].val;
            d = derm[dm]; 
						for(k = 0; k < nderive_dpop[d][ci]; k++){
							dermpopnum[dm][derive_dpop[d][ci][k]] -= derive_dpopweight[d][ci][k];
						}
            if(dermap[dm] == 0){ derlist.push_back(dm); dermap[dm] = 1;}
            break;
        }

        j++; tt = compcap[ci][j].t;
      }
    }

    if(cf != NOTALIVE){
      j = compcapdiv[cf][div]; tt = compcap[cf][j].t; while(tt < ti){ j++; tt = compcap[cf][j].t;}
      while(tt < tf){
        switch(compcap[cf][j].ty){
          case CAP:
            cap = compcap[cf][j].val;
            cpe = capprob_cpe[cap][cf];
            if(i < nind) ob = capindob[cap][i]; else ob = -1;
            if(ob >= 0){   // Observed
              if(i >= nind) emsg("Likelihood: EC6");
              Lob_pd += cpe_val[cpe]; cpe_num[cpe]++;
            }
            else{
              if(cpe < 0) emsg("Likelihood: EC7");
              Lob_pd += cpe_valoneminus[cpe]; cpe_oneminus[cpe]++;
            }
            break;

          case POP:
            p = compcap[cf][j].val; if(popmap[p] == 0) poplist.push_back(p); popmap[p]++;
            break;

          case DER:
            dm = compcap[cf][j].val;
            d = derm[dm]; 
						for(k = 0; k < nderive_dpop[d][cf]; k++){
							dermpopnum[dm][derive_dpop[d][cf][k]] += derive_dpopweight[d][cf][k];
						}
            if(dermap[dm] == 0){ derlist.push_back(dm); dermap[dm] = 1;}
            break;
        }
        j++; tt = compcap[cf][j].t;
      }
    }

    if(popfl == 1){       // Works out change in population likelihood
      for(j = 0; j < poplist.size(); j++){
        p = poplist[j];
        if(popmap[p] != 0){
          pop[p] += popmap[p];
          popmap[p] = 0;
          Lpop -= popL[p]; popL[p] = popcalc(p,pop[p]); Lpop += popL[p];
        }
      }
    }

    if(derfl == 1){       // Works out change in derived likelihood
      for(j = 0; j < derlist.size(); j++){
        dm = derlist[j]; dermap[dm] = 0;
        Lder -= derL[dm]; derL[dm] = dercalc(dm,dermpopnum[dm],param); Lder += derL[dm];
      }
    }
  }

  kmax = transchref[ci][cf].size();         // Goes through basic changes to go from ci to cf
		
  for(k = 0; k < kmax; k++){
    ref = transchref[ci][cf][k];
    jmax = transdepeqch[ref].size();
    for(j = 0; j < jmax; j++){              // Lookup table is used to find the changes
			switch(multich){
			case 0: // If only a single change is made
				eqch = transdepeqch[ref][j];
				secchange2(eqch.d,eqch.n,eqch.popnum,eqch.valch,ti,tf);
				break;
				
			case 1: //If performing multiple changes to individuals
				d = transdepeqch[ref][j].d;
				if(depmap[d] == 0){ depmap[d] = 1; depmaplist.push_back(d); depmaprefl[d].clear();}

				REFL refl; refl.t = ti; refl.ref = ref; refl.j = j; refl.sign = 1;
				depmaprefl[d].push_back(refl);
				refl.t = tf; refl.sign = -1;
				depmaprefl[d].push_back(refl);
				break;
			}
		}

    jmax = transnotdepeqch[ref].size();
    for(j = 0; j < jmax; j++){             // Lookup table is used to find the changes 
      ndeqch = transnotdepeqch[ref][j];
      d = ndeqch.d;
      dt = (tf-ti)*ndeqch.n;
      transnotdepeq_dt[d] += dt;
	    Liexp += transnotdepeq_val[d]*dt;
		}
  }
}

void Chain::secchange2(long d, long dn, double *ch_popnum, long valch, double ti, double tf)
{
	long nold, nnew, efirst, e, ee, enew, jj, eq, kd, kdi, kdf;
	double valold, valnew, dval, t, tt, fac, kdif, kdff, dt;
	
	timeprop[SEC] -= clock();
		
	eq = transdepeq[d];
	n_popnum = neq_popnum[eq];
	
	if(discon == 1){
		kdif = d*DX + ti*dfac; kdi = long(kdif);
		kdff = d*DX + tf*dfac; kdf = long(kdff);
		
		if(kdi == kdf){
			dt = dn*(tf-ti);
			depeq_disc_evdt[kdi] += dt;
			Liexp += depeq_disc_evval[kdi]*dt;
		}
		else{
			//if(valch == 1){
				for(kd = kdi+1; kd <= kdf; kd++){
					for(jj = 0; jj < n_popnum; jj++) depeq_disc_evpopnum[kd][jj] += ch_popnum[jj];
					valold = depeq_disc_evval[kd];
					if(kd == kdi+1) dval = ratecalcdep(d,depeq_disc_evpopnum[kd],param) - valold;
					valnew = valold + dval;
					Lir += depeq_disc_evn[kd]*log(valnew/valold);
					Liexp += depeq_disc_evdt[kd]*(valnew-valold);
					depeq_disc_evval[kd] = valnew;
				}
			//}
			
			dt = dn*(kdi+1-kdif)/dfac;
			depeq_disc_evdt[kdi] += dt;
			Liexp += depeq_disc_evval[kdi]*dt;
					
			for(kd = kdi+1; kd < kdf; kd++){ 
				dt = dn/dfac;
				depeq_disc_evdt[kd] += dt;
				Liexp += depeq_disc_evval[kd]*dt;
			}
			
			dt = dn*(kdff - kdf)/dfac;
			depeq_disc_evdt[kdf] += dt;
			Liexp += depeq_disc_evval[kdf]*dt;
		}
	}
	else{
		t = ti;
		e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];      // Finds where we are on the timeline
		do{
			ee = depeq_evnext[e]; tt = depeq_evt[ee]; if(tt > t) break;
			e = ee;
		}while(1 == 1);

		if(depeq_evt[e] == t){  // Simulataneous event
			valold = depeq_evval[e];
			if(valch == 1){
				for(jj = 0; jj < n_popnum; jj++) depeq_evpopnum[e][jj] += ch_popnum[jj];
				valnew = ratecalcdep(d,depeq_evpopnum[e],param);
			}
			else valnew = valold;
			nold = depeq_evn[e];
			nnew = nold+dn;

			depeq_evval[e] = valnew;
			depeq_evn[e] = nnew;

			efirst = e;
		}
		else{   // Insert new event
			enew = insertevent(d,e,t,eq,dn,ch_popnum,1);
			valold = depeq_evval[e]; nold = depeq_evn[e];
			valnew = depeq_evval[enew]; nnew = depeq_evn[enew];

			e = enew;
			efirst = -1;
		}

		dval = valnew-valold; // THIS WILL NOT ALWAYS WORK
		
		long num = 0;
		if(valch == 0){       // No change in value
			while(tt < tf){     // Does intermediary points
				Liexp += valold*dn*(tt-t);
				e = ee;
				valold = depeq_evval[ee]; 
				depeq_evn[e] += dn;

				t = tt;
				ee = depeq_evnext[e]; tt = depeq_evt[ee];
				num++;
			}
			valnew = valold; nnew = nold+dn;
		}
		else{
			while(tt < tf){     // Does intermediary points
				fac = ndepeq_ev_ev[e]; if(fac != 0) Lir += fac*log(valnew/valold);
				Liexp += (valnew*nnew - valold*nold)*(tt-t);
				e = ee;
				valold = depeq_evval[ee]; 
				
				for(jj = 0; jj < n_popnum; jj++) depeq_evpopnum[e][jj] += ch_popnum[jj];
				valnew = valold + dval;
				//valnew = ratecalcdep(d,depeq_evpopnum[e],param);
				nold = depeq_evn[e];
				nnew = nold+dn;

				depeq_evval[e] = valnew;
				depeq_evn[e] = nnew;

				t = tt;
				ee = depeq_evnext[e]; tt = depeq_evt[ee];
				num++;
			}
		}
			
		if(tf < tt){
			enew = insertevent(d,e,tf,eq,dn,ch_popnum,-1);
		}

		Lir += ndepeq_ev_ev[e]*(log(valnew)-log(valold));
		Liexp += (valnew*nnew - valold*nold)*(tf-t);
	
		if(efirst >= 0) trydeleteevent(d,efirst);
		if(tf == tt) trydeleteevent(d,ee);
	}
	timeprop[SEC] += clock();
}

// inserts event onto equation timeline
long Chain::insertevent(long d, long e, double t, long eq, long dn, double *dpopnum, long sign)  
{
  long jj, div, divf, nn, j;
  long enew, nex;

  if(ndepeq_evstack == 0){ 
    enew = ndepeq_ev;
    ndepeq_ev++;
    depeq_evn.push_back(0);
    depeq_evpopnum.push_back(new double[neq_popnum[eq]]);
    depeq_evt.push_back(0);
    depeq_evnext.push_back(0);
    depeq_evback.push_back(0);
    ndepeq_ev_ev.push_back(0);
    depeq_ev_evt.push_back(vector<double>());
    depeq_evval.push_back(0);
  }
  else{
    ndepeq_evstack--;
    enew = depeq_evstack[ndepeq_evstack];
    depeq_evstack.pop_back();
  }

  nex = depeq_evnext[e];
  depeq_evnext[enew] = nex;
  depeq_evback[nex] = enew;
  depeq_evnext[e] = enew;
  depeq_evback[enew] = e;

  nn = 0;                                              // Splits up the actual events
  j = 0;
  depeq_ev_evt[enew].clear();
  while(j < ndepeq_ev_ev[e]){
    if(depeq_ev_evt[e][j] > t){
      depeq_ev_evt[enew].push_back(depeq_ev_evt[e][j]);
      nn++;
      ndepeq_ev_ev[e]--;
      depeq_ev_evt[e][j] = depeq_ev_evt[e][ndepeq_ev_ev[e]];
      depeq_ev_evt[e].pop_back();
    }
    else j++;
  }
  ndepeq_ev_ev[enew] = nn;

  depeq_evt[enew] = t;
  depeq_evn[enew] = depeq_evn[e] + dn*sign;

  for(jj = 0; jj < neq_popnum[eq]; jj++){
		depeq_evpopnum[enew][jj] = depeq_evpopnum[e][jj] + dpopnum[jj]*sign;
	}
	
  depeq_evval[enew] = ratecalcdep(d,depeq_evpopnum[enew],param);

  div = long(onefac*t*depeqdivmax/tmax)+1;
  divf = long(onefac*depeq_evt[nex]*depeqdivmax/tmax);

  while(div <= divf){
    depeqdiv[d][div] = enew;
    div++;
  }
  return enew;
}

void Chain::trydeleteevent(long d, long e)   // Tries to delete event on the equation time line
{
  long jj, div, divf, j;
  long elast, nex;

  elast = depeq_evback[e];   // Decided if to delete event
  if(elast == -1) return;
  if(depeq_evn[elast] != depeq_evn[e]) return;

  nex = depeq_evnext[e];
  if(nex == -1) return;

  for(jj = 0; jj < n_popnum; jj++) if(depeq_evpopnum[e][jj] != depeq_evpopnum[elast][jj]) return;

  depeq_evstack.push_back(e);
  ndepeq_evstack++;

  depeq_evnext[elast] = nex;
  depeq_evback[nex] = elast;

  for(j = 0; j < ndepeq_ev_ev[e]; j++){                                  // adds on actual events
    depeq_ev_evt[elast].push_back(depeq_ev_evt[e][j]);
    ndepeq_ev_ev[elast]++;
  }
  ndepeq_ev_ev[e] = 0; depeq_ev_evt[e].clear();

  div = long(onefac*depeq_evt[elast]*depeqdivmax/tmax)+1;
  divf = long(onefac*depeq_evt[nex]*depeqdivmax/tmax);

  while(div <= divf){
    depeqdiv[d][div] = elast;
    div++;
  }
}

void Chain::multisec()  // Looks to update sections from multiple individual changes
{
	long k, d, kk, kkmax, si, dn, p, npo;
	double t, tt;
	REFL refl;
	EQCH eqch;
	
	for(k = 0; k < depmaplist.size(); k++){
		d = depmaplist[k];
		 
		if(multich == 1) sort(depmaprefl[d].begin(),depmaprefl[d].end(),comparerefl);
		 
		npo = neq_popnum[transdepeq[d]];
		for(p = 0; p < npo; p++) ch_popn[p] = 0;
		dn = 0;
	
		kkmax = depmaprefl[d].size();
		
		kk = 0;
		t = -large;
		do{
			refl = depmaprefl[d][kk];
			tt = refl.t;
			
			if(t != -large && tt > t)	secchange2(d,dn,ch_popn,1,t,tt);
			
			eqch = transdepeqch[refl.ref][refl.j];
			si = refl.sign; if(multich == -1) si *= -1;
			if(eqch.d != d) emsg("Likelihood: EC8");
			
			dn += si*eqch.n; for(p = 0; p < npo; p++)	ch_popn[p] += si*eqch.popnum[p];
		
			t = tt;
			kk++;
		}while(kk < kkmax && t < tmax);
	}
}

void Chain::multisecend()                           // Tidys variables after multisec change
{
	long k;
	
	multich = 0;
	for(k = 0; k < depmaplist.size(); k++) depmap[depmaplist[k]] = 0;
	depmaplist.clear();
}
