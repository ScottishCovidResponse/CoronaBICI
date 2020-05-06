// This file is used to check every is working correctly

struct EVI { long i; double t; long tr;};
bool compEVI(EVI lhs, EVI rhs) { return lhs.t < rhs.t; }

void Chain::test(short num)
{
}

void Chain::check(long num)                       // Used to check the algorithm is working
{
  long e, nev, j, f, cl, c, p, ce, r,	i, fev, tr, d;
  double t;

	if(simnum == 0) cout << "Check\n";
	
  timeprop[CHECK] -= clock();

  if(isnan(Lpri) || isinf(Lpri)) emsg("Check: EC1");
  if(isnan(Liinit) || isinf(Liinit)) emsg("Check: EC2");
  if(isnan(Lir) || isinf(Lir)) emsg("Check: EC3");
  if(isnan(Liexp) || isinf(Liexp)) emsg("Check: EC4");
  if(isnan(Linm) || isinf(Linm)) emsg("Check: EC5");
  if(isnan(Lev_pd) || isinf(Lev_pd)) emsg("Check: EC6");
  if(isnan(Lob_pd) || isinf(Lob_pd)) emsg("Check: EC7");
  if(isnan(Lob_st) || isinf(Lob_st)) emsg("Check: EC8");
  if(isnan(Lpop) || isinf(Lpop)) emsg("Check: EC9");

  for(i = 0; i < nindtot; i++){
		checkevseqsimp(indev[i]);
	}
	
  for(p = 0; p < npopm; p++) if(popmap[p] != 0) emsg("Check: EC10");
  for(p = 0; p < nderm; p++) if(dermap[p] != 0) emsg("Check: EC11");

  if(nindev.size() != nindtot) emsg("Check: EC12");
  if(indev.size() != nindtot) emsg("Check: EC13");

  for(i = 0; i < nind; i++){                      // Checks fixed events
    for(j = 0; j < nindfixev[i]; j++){
      fev = indfixev[i][j];
      t = fixev[fev].t;
      
      nev = nindev[i]; e = 0; while(e < nev && indev[i][e].t < t) e++;
      if(e == nev) emsg("Check: EC14");
      if(indev[i][e].t != t) emsg("Check: EC15");
      if(tra[indev[i][e].tr].dc != fixev[fev].dc) emsg("Check: EC16");
      if(tra[indev[i][e].tr].like != fixev[fev].like) emsg("Check: EC17");

      ce = fixev[fev].capev;
      if(ce >= 0){
        if(tra[indev[i][e].tr].capev != ce) emsg("Check: EC18");

        c = tra[indev[i][e].tr].ci;
        if(capevfilt[ce][c] != 1) emsg("Check: EC19");
      }
    }

    tr = indev[i][0].tr; if(tra[tr].like != indlikeenter[i]) emsg("Check: EC20");
    tr = indev[i][nindev[i]-1].tr; if(tra[tr].like != indlikeleave[i]) emsg("Check: EC21");

    for(e = 0; e < nindev[i]; e++){
      tr = indev[i][e].tr;
      if(tra[tr].capev >= 0 && tra[tr].type >= 0){
        if(capevall[tra[tr].capev] == 1){
          if(tra[tr].ci == NOTALIVE){
            if(indev[i][e].t != indfixtenter[i]) emsg("Check: EC22");
            if(tra[tr].capev != indfixentercapev[i]) emsg("Check: EC23");
          }
          else{
            if(tra[tr].cf == NOTALIVE){
              if(indev[i][e].t != indfixtleave[i]) emsg("Check: EC24");
              if(tra[tr].capev != indfixleavecapev[i]) emsg("Check: EC25");
            }
            else{
              for(j = 0; j < nindfixev[i]; j++){ 
                if(fixev[indfixev[i][j]].t == indev[i][e].t){
                  if(fixev[indfixev[i][j]].like != tra[tr].like) emsg("Check: EC26");
                  break;
                }
              }
              if(j == nindfixev[i]) emsg("Check: EC27");
            }
          }
        }
      }
    }
  }

  for(i = nind; i < nindtot; i++){
    for(e = 0; e < nindev[i]; e++){
      tr = indev[i][e].tr;
      if(tra[tr].capev >= 0 && tra[tr].type >= 0){
        if(capevall[tra[tr].capev] == 1) emsg("Check: EC27");
      }
    }
  }

  for(i = 0; i < nindtot; i++){
    if(nindev[i] != indev[i].size()) emsg("Check: EC28");

    getlifespan(i);

    if(i < nind){
      if(tent > firstobst[i]) emsg("Check: EC29");
      if(tlea < lastobst[i]) emsg("Check: EC30");
    }

    if(tbirthfl == 1){
      if(indtbirth[i] < tent-age[nage]) emsg("Check: EC31");
      if(indtbirth[i] > tent) emsg("Check: EC32");
    }
    else{ if(indtbirth[i] != tent) emsg("Check: EC33");}

    checkevseq(indev[i]);
  }
 
  checkderive();
  if(corona == 0) checkparamgrad();

  for(i = 0; i < nindtot; i++) if(nindev[i] < 2) emsg("Check: EC34");

  timeprop[CHECK] += clock();
  
	checkinflist(num);

	if(multich != 0) emsg("Check: EC35");
  for(d = 0; d < ntransdepeq; d++) if(depmap[d] != 0) emsg("Check: EC36");
	
	if(discon == 1) checklikedisc(num);
	else checklike(num);
}

void Chain::checkinflist(short num)      // Checks the lists of infected / not infected / obs are correct
{
	long i, r, tot, a, e;
	
	for(i = 0; i < nind; i++){
		r = (indinit[i]/nclassval[0])%nclassval[1];
		a = (indinit[i]/(nclassval[0]*nclassval[1]))%nag;
		
		if(nindev[i] > 2+nsettime){
			for(e = 1; e < nindev[i]; e++){
				if(tra[indev[i][e].tr].ci%nclassval[0] == HH) break;
			}
			
			if(e < nindev[i]){
				if(obsinflist[r][a][indinflist[i]] != i) emsg("Check: EC37");
			}
			else{
				if(inflist[r][a][indinflist[i]] != i) emsg("Check: EC38");
			}
		}
		else{
			if(notinflist[r][a][indinflist[i]] != i) emsg("Check: EC39");
		}
	}
	
	tot = 0;
	for(r = 0; r < nregion; r++){
		for(a = 0; a < nag; a++){
			tot += long(inflist[r][a].size()) + long(notinflist[r][a].size()) + long(obsinflist[r][a].size());
		}
	}
	
	if(tot != nind)	emsg("Check: EC40"); 
}

void checkevseq(vector<EV> &vec)       // Checks that an event sequnece is consistent and correct
{
  long e, ee, cl, inside, c, li, ctop, ti, f, k, n;
  long tr;
  double t;

  if(vec.size() == 2 && vec[0].tr ==  tranull && vec[1].tr == tranull) return;

  tr = vec[0].tr; if(tra[tr].ci != NOTALIVE) emsg("Check: EC41");
  c = tra[tr].cf;
  if(vec[0].t == 0){
    if(tr != trabeg+c) emsg("Check: EC42");
  }
  else{
    if(tr == trabeg+c) emsg("Check: EC43");
    if(tr != compiftra[NOTALIVE][c] && tr != compiftra[NOTALIVE][c]+moventra) emsg("Check: EC44");
  }
  tr = vec[long(vec.size())-1].tr; if(tra[tr].cf != NOTALIVE) emsg("Check: EC45");

  for(e = 0; e < long(vec.size())-1; e++) if(vec[e].tr < 0 || vec[e].tr >= ntra) emsg("Check: EC46");
  for(e = 0; e < long(vec.size())-1; e++){ 
    if(vec[e].t > vec[e+1].t) emsg("Check: EC47");
    if(vec[e].t > vec[e+1].t-evdtmin) emsg("Check: EC48",simnum,samp);
  }
  c = NOTALIVE;
  for(e = 0; e < vec.size(); e++){
    tr = vec[e].tr;
    if(e > 0 && tra[tr].ci != c) emsg("Check: EC49");
    c = tra[tr].cf;

    if(c != NOTALIVE){
      t = vec[e].t;
      ti = 0; while(ti < nsettime && settime[ti] <= t) ti++;
      if(compval[c][settimecl] != ti) emsg("Check: EC50");
			
      f = 0; while(f < ncapevtrange && capevtrange[f] <= t) f++;
      if(compval[c][capevcl] != f) emsg("Check: EC51");
    }
  }

	for(ti = 0; ti < nsettime; ti++){
		t = settime[ti];
		for(e = 0; e < long(vec.size())-1; e++){
			if(vec[e].t < t && vec[e+1].t > t) emsg("Check: EC52");
		}
	}
}

long checkevdt()                           // Checks event times are not too close 
{
  long e;
  for(e = 0; e < long(evnew.size())-1; e++){ if(evnew[e+1].t-evnew[e].t < evdtmin) return 1;}
  return 0;
}

void checkevseqsimp(vector<EV> &vec)       // Checks that an event sequnece is consistent and correct
{
  long e, ee, cl, inside, c, li, ctop, ti, f, k, n;
  long tr;
  double t;

  if(vec.size() == 2 && vec[0].tr ==  tranull && vec[1].tr == tranull) return;

  for(e = 0; e < long(vec.size())-1; e++){
		if(vec[e].tr < 0 || vec[e].tr >= ntra) emsg("Check: EC53");
	}
  for(e = 0; e < long(vec.size())-1; e++){
    if(vec[e].t >= vec[e+1].t) emsg("Check: EC54");
  }

  c = NOTALIVE;
  for(e = 0; e < vec.size(); e++){
    tr = vec[e].tr;
    if(e > 0 && tra[tr].ci != c) emsg("Check: EC55");
    c = tra[tr].cf;
  }
}

void Chain::checklikedisc(long num)             // Checks the likelihood based on a discretises timeline
{
	long i, kd, d, k, eq, c, p, j, jj, e, tr, kdi, kdf, s1, s2, r, a, ii, s, ddt, ddti, ddtf;
	double dd, t, tt, Lr, Lexp, dpop;
	vector<long> depeq_disc_evn_ch;  
	vector<double> depeq_disc_evdt_ch;
	vector<double> depeq_disc_evval_ch;
	vector< double* > depeq_disc_evpopnum_ch;  
	vector< vector<double> > dpopnum;

	dpopnum.resize(ncomp+1);                   // Pre calculates the change in population number
  for(c = 0; c <= ncomp; c++){
    dpopnum[c].resize(npopnum);
    for(p = 0; p < npopnum; p++){ 
      dd = 0;
      for(j = 0; j < popnumterm[p].size(); j++){
        if(popnumterm[p][j] == c) dd += popnumtermweight[p][j];
      }
      dpopnum[c][p] = dd;
    }
  }

	depeq_disc_evn_ch.resize(ntransdepeq*nDD);
	depeq_disc_evdt_ch.resize(ntransdepeq*nDD);
	depeq_disc_evval_ch.resize(ntransdepeq*nDD);
	depeq_disc_evpopnum_ch.resize(ntransdepeq*nDD);
	for(kd = 0; kd < ntransdepeq*nDD; kd++){
		depeq_disc_evn_ch[kd] = 0;  
		depeq_disc_evdt_ch[kd] = 0;
		d = kd/nDD; eq = transdepeq[d];
		depeq_disc_evpopnum_ch[kd] = new double[neq_popnum[eq]];
		for(jj = 0; jj < neq_popnum[eq]; jj++) depeq_disc_evpopnum_ch[kd][jj] = 0;
	}
	
	if(corona == 1){
		for(r = 0; r < nregion; r++){
			for(a = 0; a < nag; a++){
				s1 = obsinflist[r][a].size(); s2 = inflist[r][a].size();
				for(ii = 0; ii < s1+s2; ii++){
					if(ii < s1) i = obsinflist[r][a][ii];
					else i = inflist[r][a][ii-s1];
					
					for(d = 0; d < ntransdepeq; d++){
						eq = transdepeq[d];
						for(jj = 0; jj < neq_popnum[eq]; jj++){
							p = eq_popnum[eq][jj];
							c = tra[indev[i][0].tr].cf; 
							t = indev[i][0].t;
							for(e = 1; e < nindev[i]; e++){
								tt = indev[i][e].t;
								dpop = dpopnum[c][p];
								if(dpop != 0){
									ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
									ddtf = DDconv[long(tt*DDfac)]; if(DDt[ddtf+1] < tt) ddtf++;
									for(kd = long(d*nDD + ddti+1); kd <= long(d*nDD + ddtf); kd++){
										depeq_disc_evpopnum_ch[kd][jj] += dpop;
									}
								}
								t = tt;
								c = tra[indev[i][e].tr].cf; 
							}
						}
					}	
				}
			}
		}
	}
	else{
		for(i = 0; i < nindtot; i++){
			for(d = 0; d < ntransdepeq; d++){
				eq = transdepeq[d];
			
				c = tra[indev[i][0].tr].cf; 
				t = indev[i][0].t;
				for(e = 1; e < nindev[i]; e++){
					tt = indev[i][e].t;
					ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
					ddtf = DDconv[long(tt*DDfac)]; if(DDt[ddtf+1] < tt) ddtf++;
					for(kd = long(d*nDD + ddti+1); kd <= long(d*nDD + ddtf); kd++){
						for(jj = 0; jj < neq_popnum[eq]; jj++){
							p = eq_popnum[eq][jj];
							depeq_disc_evpopnum_ch[kd][jj] += dpopnum[c][p];
						}
					}
					t = tt;
					c = tra[indev[i][e].tr].cf; 
				}
			}	
		}
	}
 
	Lr = 0; Lexp = 0;
	for(i = 0; i < nindtot; i++){
		for(e = 1; e < nindev[i]-1; e++){
			t = indev[i][e].t;
			tr = indev[i][e].tr;
			eq = tra[tr].eq;
			d = transdepref[eq];
			switch(transdep[eq]){
			case 1:
				ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
				kd = d*nDD + ddti;
				depeq_disc_evn_ch[kd]++;
				break;
			case 0:
				Lr += log(transnotdepeq_val[d]);
				break;
			}
		}
		
		c = tra[indev[i][0].tr].cf; 
		t = indev[i][0].t;
		for(e = 1; e < nindev[i]; e++){
			tt = indev[i][e].t;	
			for(j = 0; j < ncompleave[c]; j++){ 
				tr = compleave[c][j]; 
				if(tra[tr].type == EXP_TR){
					eq = tra[tr].eq;
					d = transdepref[eq];
					switch(transdep[eq]){
					case 1:
						ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
						ddtf = DDconv[long(tt*DDfac)]; if(DDt[ddtf+1] < tt) ddtf++;
						kdi = d*nDD + ddti;
						kdf = d*nDD + ddtf;	
						if(kdi == kdf){
							depeq_disc_evdt_ch[kdi] += tt-t;
						}
						else{
							depeq_disc_evdt_ch[kdi] += DDt[ddti+1]-t;
							for(ddt = ddti+1; ddt < ddtf; ddt++){
								kd = d*nDD + ddt;
								depeq_disc_evdt_ch[kd] += DDdt[ddt];
							}
							depeq_disc_evdt_ch[kdf] += tt-DDt[ddtf];
						}
						break;
					case 0:
						Lexp += (tt-t)*transnotdepeq_val[d];
						break;
					default: emsg("Check: EC56"); break;
					}
				}
			}
			c = tra[indev[i][e].tr].cf; 
			t = tt;
		}
	}
	
	for(kd = 0; kd < ntransdepeq*nDD; kd++){
		d = kd/nDD;
		if(DDcalc[kd] == 1) depeq_disc_evval_ch[kd] = ratecalcdep(d,depeq_disc_evpopnum[kd],param);
		else depeq_disc_evval_ch[kd] = 0;
	}
	
	for(kd = 0; kd < ntransdepeq*nDD; kd++){
		d = kd/nDD;
		if(DDcalc[kd] == 1){
			Lr += depeq_disc_evn_ch[kd]*log(depeq_disc_evval_ch[kd]);
			Lexp += depeq_disc_evdt_ch[kd]*depeq_disc_evval_ch[kd];
		}
	}
	
	for(d = 0; d < ntransdepeq; d++){
		eq = transdepeq[d];
		for(i= 0; i < nDD; i++){
			kd = d*nDD+i;
		
			jj = 0;
			if(DDcalc[kd] == 1) dd = depeq_disc_evpopnum[kd][jj] - depeq_disc_evpopnum_ch[kd][jj];
			else dd = depeq_disc_evpopnum[kd][jj];
				
			if(dd > tiny || dd < -tiny) emsg("Check: EC57");
		}
	}
	
	for(kd = 0; kd < ntransdepeq*nDD; kd++){
		d = kd/nDD; eq = transdepeq[d];
		
		if(depeq_disc_evn_ch[kd] != depeq_disc_evn[kd]) emsg("Check: EC57");
		
		dd = depeq_disc_evdt_ch[kd] - depeq_disc_evdt[kd];		
		if(dd*dd > tiny) emsg("Check: EC58");
		
		depeq_disc_evdt[kd] = depeq_disc_evdt_ch[kd]; 
	
		for(jj = 0; jj < neq_popnum[eq]; jj++){
			if(DDcalc[kd] == 1){
				dd = depeq_disc_evpopnum[kd][jj] - depeq_disc_evpopnum_ch[kd][jj];
				if(dd > tiny || dd < -tiny) emsg("Check: EC59");
				depeq_disc_evpopnum[kd][jj] = depeq_disc_evpopnum_ch[kd][jj];
			}
			else{	
				if(depeq_disc_evpopnum[kd][jj] != 0) emsg("Check: EC59b");
			}
		}
		
		dd = depeq_disc_evval_ch[kd] - depeq_disc_evval[kd];
		if(dd*dd > 0.001) emsg("Check: EC60",dd);
		depeq_disc_evval[kd] = depeq_disc_evval_ch[kd];
	}

	dd = Lir-Lr; if(dd*dd > 0.001) emsg("Check: EC61",dd);
	Lir = Lr;
	
	dd = Liexp-Lexp; if(dd*dd > 0.001) emsg("Check: EC62",dd);
	Liexp = Lexp;
}

void Chain::checklike(long num)                  // Checks that the likelihood is correctly calculated
{
  long cap, tr, c, p, flag, j, ob, eq, k, d, cf, ci, stim, ce, fl;
  long nev, e, i, ee, ed, nindinitch, ncompinitch[ncompswa];
	long	numd, numnd, cev_num[ncapev], cev_oneminus[ncapev];                
  double popnum[npopnum], t, tt, pd, R, dd;
  double Lr=0, Lexp=0, Linit=0, Lpr=0, Lpopch=0, Lnm=0, Le_pd=0, Lo_pd=0, Lo_st=0;
  vector <EVI> ev;
  vector <long> ncompind;
  vector <long> stat;
  vector< vector<long> > depeq_ev_evflag;

  for(ce = 0; ce < ncapev; ce++){ cev_num[ce] = 0; cev_oneminus[ce] = 0;}

  for(i = 0; i < nindtot; i++){
    for(e = 1; e < nindev[i]-1; e++){
      EV eve = indev[i][e];
      EVI evi; evi.i = i; evi.t = eve.t; evi.tr = eve.tr;
      if(tra[eve.tr].type >= 0) ev.push_back(evi);
    }
  }
  sort(ev.begin(),ev.end(),compEVI);

  ncompind.resize(ncomp+1); for(c = 0; c < ncomp; c++) ncompind[c] = 0; ncompind[NOTALIVE] = nindtot;
  stat.resize(nindtot); for(i = 0; i < nindtot; i++) stat[i] = NOTALIVE;
 
	for(e = 0; e < ndepeq_ev; e++) depeq_ev_evflag.push_back(vector<long>());

	for(d = 0; d < ntransdepeq; d++){
		e = depeqdiv[d][0];
		eq = transdepeq[d];

		do{
			for(k = 0; k < ndepeq_ev_ev[e]; k++) depeq_ev_evflag[e].push_back(0);
			e = depeq_evnext[e];
		}while(e != -1);
	}

  cap = 0; t = 0; stim = 0;
  nev = ev.size();

  numd = 0; numnd = 0;
  for(e = 0; e <= nev; e++){
		if(e == nev) tt = tmax; else tt = ev[e].t;
    if(stim < nsettime && t >= settime[stim]) stim++;

    for(p = 0; p < npopnum; p++){
			popnum[p] = 0;
			for(i = 0; i < popnumterm[p].size(); i++){
				popnum[p] += ncompind[popnumterm[p][i]]*popnumtermweight[p][i];
			}
		}
  
    while(cap < ncap && capt[cap] < tt){          // Considers the observations
      for(i = 0; i < nindtot; i++){
        flag = 0;
        if(i < nind){
          for(j = 0; j < nindobs[i]; j++){ ob = indobs[i][j]; if(obscap[ob] == cap){ flag = 1; break;}}
        }
        c = stat[i];

        if(stat[i] != NOTALIVE && flag == 0){
          eq = capprobeqn[cap][c];
          if(eq >= 0){
            pd = calcobsprob(eq);
            if(pd == 1) emsg("Check: EC61");
          }
        }

        if(stat[i] == NOTALIVE){
        }
        else{
          eq = capprobeqn[cap][c];
          if(eq >= 0){
            pd = calcobsprob(eq);
            if(flag == 1){ numd++; Lo_pd += log(pd);}
            else{ numnd++; Lo_pd += log(1-pd);}
          }
          else{
          }
        }
      }
      cap++;
    }
 
    R = 0;
    for(i = 0; i < nindtot; i++){                    // Transitions
      c = stat[i];
      if(c != NOTALIVE){
        for(j = 0; j < ncompleave[c]; j++){ tr = compleave[c][j]; if(tra[tr].type == EXP_TR) R += ratecalc(tra[tr].eq,popnum,param);}
      }
    }

    c = NOTALIVE;                                    // Sources
    for(j = 0; j < ncompleave[c]; j++){
      tr = compleave[c][j];
      if(compval[tra[tr].cf][settimecl] == stim) R += ratecalc(tra[tr].eq,popnum,param);
    }

    Lexp += R*(tt-t);

    t = tt;
    if(e == nev) break;

    i = ev[e].i;
    tr = ev[e].tr;

    ce = tra[tr].capev;                        // Capture of events

    if(ce >= 0){
      if(i < nind){
        fl = 0;
        for(j = 0; j < nindfixev[i]; j++){ if(fixev[indfixev[i][j]].t == t && fixev[indfixev[i][j]].capev == ce){ fl = 1; break;}}
        if(t == indfixtenter[i] && indfixentercapev[i] == ce) fl = 1;
        if(t == indfixtleave[i] && indfixleavecapev[i] == ce) fl = 1;

        if(fl == 1) cev_num[ce]++; else cev_oneminus[ce]++;
      }
      else cev_oneminus[ce]++;
    }

    if(tra[tr].type == EXP_TR){
      eq = tra[tr].eq;
      if(transdep[eq] == 1){
        if(tra[tr].like == 1) Lr += log(ratecalc(eq,popnum,param));

        d = transdepref[eq];
        ed = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];
        do{
          ee = depeq_evnext[ed]
              ; if(ee == -1) emsg("Check: EC62");
          if(ee == -1 || depeq_evt[ee] >= t) break;
          ed = ee;
        }while(1 == 1);

        for(k = 0; k < ndepeq_ev_ev[ed]; k++) if(depeq_ev_evt[ed][k] == t) break;
        if(k == ndepeq_ev_ev[ed]){
          if(tra[tr].like == 1) emsg("Check: EC63");
        }
        else{
          if(tra[tr].like == 0) emsg("Check: EC64");

          dd = ratecalc(eq,popnum,param)-depeq_evval[ed];
          if(dd*dd > tiny) emsg("Check: EC65");
          if(depeq_ev_evflag[ed][k] == 1) emsg("Check: EC66");
          depeq_ev_evflag[ed][k] = 1;
        }
      }
      else{
        if(tra[tr].like == 1) Lr += log(ratecalc(eq,popnum,param));
      }
    }

    ci = stat[i];
    if(tra[tr].ci != ci) emsg("Check: EC67",num,samp);
	
    cf = tra[tr].cf; stat[i] = cf;
    ncompind[ci]--; if(ncompind[ci] < 0) emsg("Check: EC68");
    ncompind[cf]++;
  }

  if(capevfl == 0){
    if(Lev_pd != 0) emsg("Check: EC69"); 
  }
  else{
    for(ce = 0; ce < ncapev; ce++){
      if(cev_num[ce] != capev_num[ce]) emsg("Check: EC70");
      if(cev_oneminus[ce] != capev_oneminus[ce]) emsg("Check: EC71");
      pd = calcobsprob(capevprobeqn[ce]);
      if(capev_val[ce] != log(pd)) emsg("Check: EC72");
      if(capev_valoneminus[ce] != log(1-pd)) emsg("Check: EC73");

      if(capev_num[ce] > 0){
        if(pd == 0) emsg("Check: EC74");
        Le_pd += capev_num[ce]*log(pd);
      }
      if(capev_oneminus[ce] > 0){
        if(pd == 1) emsg("Check: EC75");
        Le_pd += capev_oneminus[ce]*log(1-pd);
      }
    }
  }

  for(d = 0; d < ntransdepeq; d++){
    e = depeqdiv[d][0];
    eq = transdepeq[d];

    do{
      for(k = 0; k < ndepeq_ev_ev[e]; k++){
        if(depeq_ev_evflag[e][k] == 0) emsg("Check: EC76");
      }
      e = depeq_evnext[e];
    }while(e != -1);
  }

  Lpr = priorcalc();

  nwrong = 0;
  Lo_st = 0;
  for(i = 0; i < nind; i++){
    if(i < nindtot){
      c = NOTALIVE; e = 0;
      for(j = 0; j < nindobs[i]; j++){
        ob = indobs[i][j];

        t = obst[ob]; 
        while(e < nindev[i] && indev[i][e].t < t){ 
					if(tra[indev[i][e].tr].type >= 0) c = tra[indev[i][e].tr].cf; e++;
				}
        if(c == NOTALIVE){ Lo_st += notobsdL; nwrong++;}
        else{
          eq = obsprobeqn[ob][c];
          if(eq == -1){ Lo_st += notobsdL; nwrong++;}
          else Lo_st += log(calcobsprob(eq));
        }
      }
    }
  }
 
  Lpopch = 0;                                               // Checks popm
  for(p = 0; p < npopm; p++){ 
    short numm = 0;
    t = popmt[p];
    for(i = 0; i < nindtot; i++){
      c = NOTALIVE; e = 0; while(e < nindev[i] && indev[i][e].t < t){ c = tra[indev[i][e].tr].cf; e++;}
      if(popmcomp[p][c] == 1) numm++;
    }
    if(numm != pop[p]) emsg("Check: EC77");

    dd = popL[p]-popcalc(p,numm); if(dd > tiny) emsg("Check: EC78");
    Lpopch += popcalc(p,numm);
  }

  for(k = 0; k < nnmeq; k++){
    eq = nmeq[k];
    if(nmeq_val[eq] != calculatenotdep(eq,param)) emsg("Check: EC79");
  }

  Lnm = 0;                               // Caluculates the likelihood for the non Markovian transitions
	for(i = 0; i < nindtot; i++) Lnm += likenm(i);   

  nindinitch = 0;
  for(c = 0; c < ncompswa; c++) ncompinitch[c] = 0;
  for(i = 0; i < nindtot; i++){ 
		if(indev[i][0].t == 0){ nindinitch++; ncompinitch[tra[indev[i][0].tr].cf%ncompswa]++;}
	}
  for(c = 0; c < ncompswa; c++){ if(ncompinitch[c] != ncompinit[c]) emsg("Check: EC80");}
  if(nindinitch != nindinit) emsg("Check: EC81");

  Linit = logsum[nindinit]; for(c = 0; c < ncompswa; c++) Linit += ncompinit[c]*log(probinit[c]);

  dd = Lir-Lr; if(dd*dd > 0.00001) emsg("Check: EC82");
  dd = Linm-Lnm; if(dd*dd > 0.0001) emsg("Check: EC83");
  Linm = Lnm;

  dd = Liinit-Linit; if(dd*dd > tiny) emsg("Check: EC84");
  dd = Lev_pd-Le_pd; if(dd*dd > tiny) emsg("Check: EC85");
  dd = Lob_pd-Lo_pd; if(dd*dd > tiny) emsg("Check: EC86");
  dd = Lob_st-Lo_st; if(dd*dd > 0.0001) emsg("Check: EC87");
  dd = Lpop-Lpopch; if(dd*dd > 0.0001) emsg("Check: EC88");
  dd = Lpri-Lpr; if(dd*dd > tiny) emsg("Check: EC89");
	dd = Liexp-Lexp; if(dd*dd > tiny) emsg("Check: EC90");
}

void Chain::checkderive()               // Checks the likelihood for derived observations is correct
{
  long d, dm, i, j, k, e, tr, p, c;
  double t, popnum[npopnum], sum, L, LL, dd;
  vector <EV> ev;
  vector <long> nc;

  for(i = 0; i < nindtot; i++){
    for(e = 0; e < nindev[i]; e++) ev.push_back(indev[i][e]);
  }
  sort(ev.begin(),ev.end(),compareev);

  nc.resize(ncomp+1); for(c = 0; c < ncomp; c++) nc[c] = 0; nc[NOTALIVE] = nindtot_sim;

  L = 0; j = 0;
  for(dm = 0; dm < nderm; dm++){
    t = dermt[dm];
    while(j < ev.size() && ev[j].t <= t){
      tr = ev[j].tr; nc[tra[tr].ci]--; nc[tra[tr].cf]++;
      j++;
    }

    for(p = 0; p < npopnum; p++){
      sum = 0; 
			for(k = 0; k < popnumterm[p].size(); k++){
				sum += popnumtermweight[p][k]*nc[popnumterm[p][k]];
			}
      popnum[p] = sum;
    }

    d = derm[dm];
    LL = dercalc(dm,popnum,param);
    dd = LL-derL[dm]; if(dd*dd > tiny) emsg("Check: EC91");
    L += LL;
  }
  dd = L-Lder; if(dd*dd > tiny) emsg("Check: EC92");

  for(dm = 0; dm < nderm; dm++){ if(dermap[dm] != 0) emsg("Check: EC92");}
}

void Chain::checkparamgrad()               // Checks that paramgrad is working correctly
{
  long p;
  double val, dd, Lst, Lup, Ldo, d, pmean2, pvar2, pgrad, pcurve;

  Lst = L();
  for(p = 0; p < nparam; p++){
    val = param[p]; dd = val/1000; if(dd == 0) dd = 0.0000001;
    paramgrad(p);
    changeparam(1,p,val+dd); Lup = L();
    changeparam(1,p,val-dd); Ldo = L();

    pgrad = (Lup-Ldo)/(2*dd);
    pcurve = (Lup+Ldo-2*Lst)/(dd*dd);
   
    if(pcurve > 0){ pmean2 = 0; pvar2 = 0;}
    else{
      pmean2 = val - pgrad/pcurve;
      pvar2 = -1.0/pcurve;
      if(pmean2*pmean2 > 100*val*val){ pmean2 = 0; pvar2 = 0;}
    }
    d = (pmean2 - pmean)/val; 
    if(d*d > 0.0001) emsg("Check: EC93");
    
    d = (pvar2 - pvar)/(val*val); 
    if(d*d > 0.0001) emsg("Check: EC94");

    changeparam(1,p,val);
  }
}
