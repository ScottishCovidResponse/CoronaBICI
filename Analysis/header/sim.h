// This file simlates from the model

double *sim_popnum;

vector< vector< vector<long> > > sim_depref, sim_notdepref;    // [i][d][#] Reference to sim_transdeplist
vector< vector<long> > sim_depref_d, sim_depref_k;
vector< vector<long> > sim_notdepref_d, sim_notdepref_k;

vector< vector<long> > sim_transdeplist, sim_transnotdeplist;  // [d] List of ind. attached to an eq

vector<long> ncomppopnumch;                                    // How a compartment affects the pop nums
vector< vector<long> > comppopnumch;
vector< vector<double> > comppopnumchweight;

vector<double> sim_Rtransdep, sim_Rtransnotdep;                // Stores transition rates

const long ntdiv = 1000;                                        // Discretises the future events

double tdivmin, tdivmax;
vector <vector <long> > tdiv;

struct FUTEV { long tr; long cl; long i; double tnow; double t;}; // Structure for a future event

vector <FUTEV> futev;                                           // Stores future events

vector< vector<long> > sourcetra;                               // [settime][#] A list of all the sources
vector<long> nsourcetra;

long nindtot_sim;                                               // The total number of simulated inds

vector <double> indtbirth_sim;                                  // Individual data
vector<long> nindev_sim;
vector <vector <EV> > indev_sim;

vector<long> sim_statind;                                       // The compartments in which inds reside

void add_transnotdeplist(long d, long i);
void rem_transnotdeplist(long d, long k);
void add_transdeplist(long d, long i);
void rem_transdeplist(long d, long k);

void Chain::sim(double tmin)              // Simulates from the model starting at time tmin
{
  long ob, c, cc, ag, p, j, d, cl, ti, index, ts, fl, flag, e, ee, dep, k, ci, cf, flindex=0;
  long i, tr, trr, eq, kmax, di, per, pernew;
  double  t, R, z, tmi, tma, tadd, tt, kshape, mean, lam, probsum[ncompswa], sum, indext;
  vector <double> Rst;
  vector <long> postr;

	timeprop[SIMTIME] -= clock();
	
  tdivmin = tmin; tdivmax = tmax2; tdiv.resize(ntdiv); for(i = 0; i < ntdiv; i++) tdiv[i].clear();

	//for(p = 0; p < nparam; p++) cout << paramname[p] << " " << param[p] << " \n";
	
  futev.clear();

  switch(simon){
    case 2:
			emsg("Sim: EC1");
      break;

    case 1:             // Sets up a simulation
      nindtot_sim = nind;
      indev_sim.clear();
      nindev_sim.resize(nind); indev_sim.resize(nind); 
			indtbirth_sim.resize(nind); sim_statind.resize(nind);
			
			if(corona == 1){ tmin = 40; indext = 39.9;}
			
      for(i = 0; i < nind; i++){
				if(corona == 1) cc = indinit[i];
				else{
					if(nindobs[i] != 1) emsg("Sim: EC2");
					ob = indobs[i][0];
					cc = -1;
					for(c = 0; c < ncomp; c++){
						if(obsprobeqn[ob][c] != -1){
							if(cc == -1) cc = c; else emsg("Sim: EC3");
						}
					}
				}
				
        nindev_sim[i] = 1;
        EV evbeg; evbeg.t = 0; evbeg.tr = trabeg+cc; indev_sim[i].push_back(evbeg);

				if(corona == 1 && flindex < 3 && compname[cc] == "S_T<80"){
					EV evindex; evindex.t = indext; evindex.tr = compiftra[cc][cc+EE]; indev_sim[i].push_back(evindex);
					flindex++;
					cc = cc+EE;
				}
        sim_statind[i] = cc;

        ag = compval[cc][agecl];
        if(ag == 0) indtbirth_sim[i] = -ran()*age[ag]; 
				else indtbirth_sim[i] = -(age[ag-1] + ran()*(age[ag]-age[ag-1]));

        if(cc != NOTALIVE){ for(cl = 0; cl < nclass-1; cl++) addfuture(i,cc,cl,tmin);}
      }
			break;

    case 0:
      if(tmin > 0){     // When doing inference
        nindtot_sim = nindtot;
        nindev_sim.resize(nindtot_sim); indev_sim.resize(nindtot_sim); 
				indtbirth_sim.resize(nindtot_sim); sim_statind.resize(nindtot_sim);
        for(i = 0; i < nindtot_sim; i++){
          indev_sim[i] = indev[i]; nindev_sim[i] = nindev[i]; indtbirth_sim[i] = indtbirth[i];

          e = nindev_sim[i]-1;
          tr = indev_sim[i][e].tr;

          if(tra[tr].type == END_TR){
            for(j = 0; j < tra[tr].ntraend; j++){
              trr = tra[tr].traend[j];
              cl = tra[trr].cl;
              ee = e-1; while(ee > 0 && tra[indev_sim[i][ee].tr].cl != cl) ee--;
              tt = indev_sim[i][ee].t;

              switch(tra[trr].type){
                case GAMMA_TR:
                  mean = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
                  tadd = tt + gammasamptail(kshape,kshape/mean,tmin-tt);
                  break;

                case WEI_TR:
                  lam = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
                  tadd = tt + weibullsamptail(lam,kshape,tmin-tt);
                  break;

                case FIXED_TR: emsg("Sim: EC4"); break;
              }
              if(tadd < tmin) emsg("Sim: EC5");
              if(tadd < tmax2) addfuturetra(i,tr,tt,tadd);
            }
            indev_sim[i].pop_back(); e--; nindev_sim[i]--;
          }

          c = tra[indev_sim[i][nindev_sim[i]-1].tr].cf;
          sim_statind[i] = c;

          if(c != NOTALIVE){  // adds in any age of set time future transitions
					  for(cl = nclass-3; cl < nclass-1; cl++) addfuture(i,c,cl,tmin);
			    }
        }
      }
      else{     // When doing PPC
        nindtot_sim = 0;
        nindev_sim.clear(); indev_sim.clear(); indtbirth_sim.clear(); sim_statind.clear();
        for(i = 0; i < nindtot; i++){
          if(indev[i][0].t == 0){
            indev_sim.push_back(vector<EV> ());
            indev_sim[nindtot_sim].push_back(indev[i][0]);
            nindev_sim.push_back(1);
            indtbirth_sim.push_back(indtbirth[i]);

            cc = tra[indev[i][0].tr].cf;
            sim_statind.push_back(cc);

            if(cc != NOTALIVE){ for(cl = 0; cl < nclass-1; cl++) addfuture(nindtot_sim,cc,cl,tmin);}
            nindtot_sim++;
          }
        }
      }
      break;
  }

  for(p = 0; p < npopnum; p++) sim_popnum[p] = 0;

	if(corona == 0){
		sim_depref.clear(); sim_depref.resize(nindtot_sim);
		sim_notdepref.clear(); sim_notdepref.resize(nindtot_sim);
		for(i = 0; i < nindtot_sim; i++){  
			sim_depref[i].resize(ntransdepeq);
			sim_notdepref[i].resize(ntransnotdepeq);
		}
	}
	else{
		sim_depref_d.clear(); sim_depref_d.resize(nindtot_sim);
		sim_depref_k.clear(); sim_depref_k.resize(nindtot_sim);
		sim_notdepref_d.clear(); sim_notdepref_d.resize(nindtot_sim);
		sim_notdepref_k.clear(); sim_notdepref_k.resize(nindtot_sim);
	}
	
  for(i = 0; i < nindtot_sim; i++){        // Finds the composition in population at tmax
    c = sim_statind[i];
    if(c != NOTALIVE){
      for(j = 0; j < ncomppopnumch[c]; j++) sim_popnum[comppopnumch[c][j]] += comppopnumchweight[c][j];
    }
  }

  // The total rate of exponential transitions is worked out for each of the equations 

  sim_transdeplist.clear(); sim_transnotdeplist.clear();
  sim_transdeplist.resize(ntransdepeq); sim_transnotdeplist.resize(ntransnotdepeq);

	sim_Rtransdep.resize(ntransdepeq);
  for(d = 0; d < ntransdepeq; d++) sim_Rtransdep[d] = ratecalc(transdepeq[d],sim_popnum,param);

  sim_Rtransnotdep.resize(ntransnotdepeq);
  for(d = 0; d < ntransnotdepeq; d++) sim_Rtransnotdep[d] = ratecalc(transnotdepeq[d],sim_popnum,param);

  for(i = 0; i < nindtot_sim; i++){         // Finds the composition in popluation at tmax
    c = sim_statind[i]; if(c != NOTALIVE) sim_indchange(i,NOTALIVE,c);
  }
		 
	if(1 == 0){
		for(d = 0; d < ntransnotdepeq; d++){
			cout << sim_transnotdeplist[d].size() << "   "  << eqnstr[transnotdepeq[d]] << " " 
						<< sim_Rtransnotdep[d] << "  notdep\n";
		}

		for(d = 0; d < ntransdepeq; d++){
			cout << sim_transdeplist[d].size() << " "  << eqnstr[transdepeq[d]] << "  dep\n";
		}
	}
	
  ti = 0; index = 0;

  Rst.resize(ntransdepeq+ntransnotdepeq);

  tmi = tmin;
  ts = 0; while(ts < settime.size() && tmin >= settime[ts]) ts++;

  per = 0;
  while(ts <= settime.size()){
    if(ts < settime.size()) tma = settime[ts]; else tma = tmax2;

    if(ts > 0){
      for(d = 0; d < ntransdepeq; d++){
        k = 0; kmax = sim_transdeplist[d].size();
        while(k < kmax){
          if(sim_transdeplist[d][k] < 0){ rem_transdeplist(d,k); kmax--;}
          else k++;
        }
      }

      for(d = 0; d < ntransnotdepeq; d++){
        k = 0; kmax = sim_transnotdeplist[d].size();
        while(k < kmax){
          if(sim_transnotdeplist[d][k] < 0){ rem_transnotdeplist(d,k); kmax--;}
          else k++;
        }
      }
    }

    for(k = 0; k < nsourcetra[ts]; k++){          // Adds in source transitions
      tr = sourcetra[ts][k]; eq = tra[tr].eq;
      d = transdepref[eq];
      switch(transdep[eq]){
        case 1: sim_transdeplist[d].push_back(-1-tr); break;
        case 0: sim_transnotdeplist[d].push_back(-1-tr); break;
        default: emsg("Sim: EC6"); break;
      }
    }

    t = tmi;
    do{
		  if(simon == 1){ 
				pernew = long(100*(samp + t/tmax2)/nsamp); 
				if(pernew > per){ per = pernew; cout << "3|" << per << "|\n"; cout.flush();}
			}
			 
      R = 0;
      for(d = 0; d < ntransdepeq; d++){
				R += long(sim_transdeplist[d].size())*sim_Rtransdep[d]; Rst[d] = R;
			}
			
      for(d = 0; d < ntransnotdepeq; d++){ 
				R += long(sim_transnotdeplist[d].size())*sim_Rtransnotdep[d]; Rst[ntransdepeq+d] = R;
			}
		
			if(1 == 0){
				for(d = 0; d < ntransdepeq; d++){
					cout << sim_transdeplist[d].size() << " "<< eqnstr[transdepeq[d]]
								<< " " << sim_Rtransdep[d] << "  dep\n";
				}
      	for(d = 0; d < ntransnotdepeq; d++){
					cout << long(sim_transnotdeplist[d].size()) << " " 
								<< eqnstr[transnotdepeq[d]] << " " << sim_Rtransnotdep[d] << "  not dep\n";
				}
			}
			
      //if(checkon == 1) sim_check(R,t);

      if(R == 0) t = tma;
      else{
        t -= log(ran())/R;
        if(t > tma) t = tma;
      }
			
		  postr.clear();

      fl = 0;
      while(ti < ntdiv){        // Checks to see if a future event has occured
        while(index < tdiv[ti].size() && futev[tdiv[ti][index]].t < t){
          j = tdiv[ti][index]; index++;
          cl = futev[j].cl;
          tr = futev[j].tr;
          i = futev[j].i;
          c = sim_statind[i];
				  if(c != NOTALIVE){
            if(compval[sim_statind[i]][cl] == tra[tr].i){
              flag = 0; e = nindev_sim[i]-1; 
							while(e >= 0 && indev_sim[i][e].t > futev[j].tnow){ // Checks for intervening events
								if(tra[indev_sim[i][e].tr].cl == cl) flag = 1; e--;
							} 
              if(flag == 0){
                cf = c + (tra[tr].f-tra[tr].i)*classmult[cl];
                tr = compiftra[c][cf]; if(tr == -1) emsg("Sim: EC7");
                postr.push_back(tr); t = futev[j].t; fl = 1;
								break;
              }
            }
          }
        }
        if(fl == 1 || index < tdiv[ti].size() || tdivmin + (ti+1)*(tdivmax-tdivmin)/ntdiv > t) break;
        ti++; index = 0;
      };

      if(fl == 0 && t < tma){
        z = ran()*R; d = 0; while(d < ntransdepeq+ntransnotdepeq && z > Rst[d]) d++;

        if(d == ntransdepeq+ntransnotdepeq) emsg("Sim: EC8");

        if(d < ntransdepeq){                      // Selects individual
					dep = 1; i = sim_transdeplist[d][long(ran()*long(sim_transdeplist[d].size()))];
				}                  
        else{
					d -= ntransdepeq; dep = 0; 
					i = sim_transnotdeplist[d][long(ran()*long(sim_transnotdeplist[d].size()))];
				}

        if(i < 0){  // An individual is born
          tr = -(i+1);
          c = tra[tr].cf;
          i = nindtot_sim; nindtot_sim++; 
          if(nindtot_sim > indmax){
            stringstream ss; ss << "Maximum number of " << indmax << " individuals exceeded.";
            emsg(ss.str());
          }

					if(corona == 0){
						sim_depref.push_back(vector< vector<long> > ()); sim_depref[i].resize(ntransdepeq);
						sim_notdepref.push_back(vector< vector<long> > ()); sim_notdepref[i].resize(ntransnotdepeq);
					}
					else{
						sim_depref_d.push_back(vector<long> ());
						sim_depref_k.push_back(vector<long> ());
						sim_notdepref_d.push_back(vector<long> ());
						sim_notdepref_k.push_back(vector<long> ());
					}
					
          nindev_sim.push_back(0);
          indev_sim.push_back(vector<EV> ());
          sim_statind.push_back(NOTALIVE);

          ag = compval[c][agecl];
          if(ag == 0) indtbirth_sim.push_back(t-ran()*age[ag]); 
					else indtbirth_sim.push_back(t-(age[ag-1] + ran()*(age[ag]-age[ag-1])));

          postr.push_back(tr);
        }
        else{
          c = sim_statind[i]; if(c == NOTALIVE) emsg("Sim: EC9");
          if(dep == 1){ 
						for(k = 0; k < ncompleavedep[c]; k++){ 
							tr = compleavedep[c][k]; if(transdepref[tra[tr].eq] == d) postr.push_back(tr);
						}
					}
          else{
						for(k = 0; k < ncompleavenotdep[c]; k++){ 
							tr = compleavenotdep[c][k]; if(transdepref[tra[tr].eq] == d) postr.push_back(tr);
						}
					}
        }
      }

      if(t >= tma){ t = tma; break;}

			if(postr.size() == 0) emsg("Sim: EC10");
      j = long(ran()*postr.size());

      tr = postr[j];

      nindev_sim[i]++;
      EV ev; ev.t = t; ev.tr = tr;
      indev_sim[i].push_back(ev);

      ci = tra[tr].ci; cf = tra[tr].cf;
      if(sim_statind[i] != ci) emsg("Sim: EC11");
 
      if(ci != NOTALIVE){ 
				for(j = 0; j < ncomppopnumch[ci]; j++){
					sim_popnum[comppopnumch[ci][j]] -= comppopnumchweight[ci][j];
				}
			}
      if(cf != NOTALIVE){ 
				for(j = 0; j < ncomppopnumch[cf]; j++){
					sim_popnum[comppopnumch[cf][j]] += comppopnumchweight[cf][j];
				}
			}
      sim_statind[i] = cf;

      sim_indchange(i,ci,cf);
    
      if(cf != NOTALIVE){
        if(ci == NOTALIVE){ for(cl = 0; cl < nclass-1; cl++) addfuture(i,cf,cl,t);}
        else addfuture(i,cf,tra[tr].cl,t);
      }
    }while(1 == 1);

    ts++; tmi = tma;
  }

  for(i = 0; i < nindtot_sim; i++){        // Finds the composition in popluation at tmax
    c = tra[indev_sim[i][nindev_sim[i]-1].tr].cf;
    if(c != NOTALIVE){ 
			EV evend; evend.t = tmax2; evend.tr = traend+c; indev_sim[i].push_back(evend); nindev_sim[i]++;
		}
  }
	
	timeprop[SIMTIME] += clock();
}

void Chain::addfuture(long i, long c, long cl, double tstart)  // When simulating adds future events 
{
  long k, j, loop, loopmax = 100;
  long tr;
  double mean, shape, lam, kk, dt, tadd;

  for(k = 0; k < ncompclleave[c][cl]; k++){
    tr = compclleave[c][cl][k];
		switch(tra[tr].type){
      case EXP_TR: emsg("Sim: EC12"); break;

      case FIXED_TR:
        dt = ratecalc(tra[tr].eq,sim_popnum,param);
        break;

      case GAMMA_TR:
        mean = calculatenotdep(tra[tr].eq,param);
        shape = calculatenotdep(tra[tr].eqshape,param);
        dt = gammasamp(shape,shape/mean);
        break;

      case WEI_TR:
        lam = calculatenotdep(tra[tr].eq,param);
        kk = calculatenotdep(tra[tr].eqshape,param);
        dt = weibullsamp(lam,kk); loop++;
        break;

      case GROW_TR:
        if(cl != agecl) emsg("Sim: EC13");
        j = compval[c][agecl];
        if(j < nage) dt = indtbirth_sim[i]+age[j] - tstart-tiny;
        else dt = large;
        break;

      case SETTIME_TR:
			  if(cl != settimecl) emsg("Sim: EC14");
        j = compval[c][settimecl];
				
        if(j < settime.size()) dt = settime[j] - tstart-tiny;
        else dt = large;
        break;
    }

    if(dt < 0){
			emsg("Sim: EC15");
		}
		
    tadd = tstart + dt;
    if(tadd < tdivmin) emsg("Sim: EC16");
	   
    if(tadd < tdivmax) addfuturetra(i,tr,tstart,tadd);
  }
}

void Chain::addfuturetra(long i, long tr, double tstart, double tadd) 
{
  long  di, j;

  di = long(ntdiv*(tadd-tdivmin)/(tdivmax-tdivmin)); if(di < 0 || di >= ntdiv) emsg("Sim: EC17");
	j = tdiv[di].size()-1; while(j >= 0 && tadd < futev[tdiv[di][j]].t) j--;
	if(j == tdiv[di].size()-1) tdiv[di].push_back(futev.size());
	else tdiv[di].insert(tdiv[di].begin()+j+1,futev.size());	

  //j = 0; while(j < tdiv[di].size() && tadd > futev[tdiv[di][j]].t) j++;
  //tdiv[di].insert(tdiv[di].begin()+j,futev.size());

  FUTEV fev; 
	fev.tr = tr; fev.cl = tra[tr].cl; fev.t = tadd; fev.tnow = tstart; fev.i = i; futev.push_back(fev);
}

void Chain::sim_indchange(long i, long ci, long cf)    // When an individual changes in the simulation
{
  long k, kmax, kk, kkmax, j, jmax, d, dn, dd;
  long ref;
  EQCH eqch;
  NDEQCH ndeqch;

  kmax = transchref[ci][cf].size();  // Goes through basic changes to go from ci to cf
  for(k = 0; k < kmax; k++){
    ref = transchref[ci][cf][k];
    jmax = transdepeqch[ref].size();
    for(j = 0; j < jmax; j++){       // Lookup table is used to find the changes in the undelying dep eqs
      eqch = transdepeqch[ref][j]; d = eqch.d; dn = eqch.n;
      if(dn != 0){
        if(dn > 0){  // Adding an individual
          for(dd = 0; dd < dn; dd++) add_transdeplist(d,i);
        }
        else{        // Removing an individual
          for(dd = 0; dd < -dn; dd++){
						if(corona == 0){
							k = sim_depref[i][d][long(sim_depref[i][d].size())-1];							
							sim_depref[i][d].pop_back();
						}
						else{
							kk = 0; kkmax = sim_depref_d[i].size(); 
							while(kk < kkmax && sim_depref_d[i][kk] != d) kk++;
							if(kk == kkmax) emsg("Sim: EC17b");
							k = sim_depref_k[i][kk];
							
							if(kkmax > 1){
								sim_depref_d[i][kk] = sim_depref_d[i][kkmax-1];
								sim_depref_k[i][kk] = sim_depref_k[i][kkmax-1];
							}
							sim_depref_d[i].pop_back();
							sim_depref_k[i].pop_back();
						}
						
            if(sim_transdeplist[d][k] != i) emsg("Sim: EC18");
            rem_transdeplist(d,k);
          }
        }
      }

      if(eqch.valch == 1) sim_Rtransdep[d] = ratecalc(transdepeq[d],sim_popnum,param);
    }

    jmax = transnotdepeqch[ref].size();
    for(j = 0; j < jmax; j++){  // Lookup table is used to find the changes in the undelying not dep. eqs
      ndeqch = transnotdepeqch[ref][j];
      d = ndeqch.d; dn = ndeqch.n;
      if(dn != 0){
        if(dn > 0){ // Adding an individual
          for(dd = 0; dd < dn; dd++) add_transnotdeplist(d,i);
        }
        else{       // Removinging an individual
          for(dd = 0; dd < -dn; dd++){
						if(corona == 0){
							k = sim_notdepref[i][d][long(sim_notdepref[i][d].size())-1];
							sim_notdepref[i][d].pop_back();
						}
						else{
							kk = 0; kkmax = sim_notdepref_d[i].size(); 
							while(kk < kkmax && sim_notdepref_d[i][kk] != d) kk++;
							if(kk == kkmax) emsg("Sim: EC17b");
							k = sim_notdepref_k[i][kk];
							
							if(kkmax > 1){
								sim_notdepref_d[i][kk] = sim_notdepref_d[i][kkmax-1];
								sim_notdepref_k[i][kk] = sim_notdepref_k[i][kkmax-1];
							}
							sim_notdepref_d[i].pop_back();
							sim_notdepref_k[i].pop_back();
						}
						
            if(sim_transnotdeplist[d][k] != i) emsg("Sim: EC19");
            rem_transnotdeplist(d,k);
          }
        }
      }
    }
  }
}

void add_transnotdeplist(long d, long i)
{
	if(corona == 0){
		sim_notdepref[i][d].push_back(sim_transnotdeplist[d].size());
	}
	else{
		sim_notdepref_d[i].push_back(d);
		sim_notdepref_k[i].push_back(sim_transnotdeplist[d].size());
	}
	sim_transnotdeplist[d].push_back(i);
}

void rem_transnotdeplist(long d, long k)
{
  long kmax, kk, jj, jjmax;

  kmax = sim_transnotdeplist[d].size();
  if(k < kmax-1){ 
    sim_transnotdeplist[d][k] = sim_transnotdeplist[d][kmax-1];
    kk = sim_transnotdeplist[d][k]; 
    if(kk >= 0){
			if(corona == 0){
				jjmax = sim_notdepref[kk][d].size();
				for(jj = 0; jj < jjmax; jj++){
					if(sim_notdepref[kk][d][jj] == kmax-1){ sim_notdepref[kk][d][jj] = k; break;}
				}
				if(jj == jjmax) emsg("Sim: EC20");
			}
			else{
				jjmax = sim_notdepref_d[kk].size();
				for(jj = 0; jj < jjmax; jj++){
					if(sim_notdepref_d[kk][jj] == d && sim_notdepref_k[kk][jj] == kmax-1){
						sim_notdepref_k[kk][jj] = k; break;
					}
				}
				if(jj == jjmax) emsg("Sim: EC21");
			}
    }
  }
  sim_transnotdeplist[d].pop_back();
}

void add_transdeplist(long d, long i)
{
	if(corona == 0){
		sim_depref[i][d].push_back(sim_transdeplist[d].size());
	}
	else{
		sim_depref_d[i].push_back(d);
		sim_depref_k[i].push_back(sim_transdeplist[d].size());
	}
	sim_transdeplist[d].push_back(i);
}

void rem_transdeplist(long d, long k)
{
  long kmax, kk, jj, jjmax;

  kmax = sim_transdeplist[d].size();
  if(k < kmax-1){ 
    sim_transdeplist[d][k] = sim_transdeplist[d][kmax-1];
    kk = sim_transdeplist[d][k]; 
    if(kk >= 0){
			if(corona == 0){
				jjmax = sim_depref[kk][d].size();
				for(jj = 0; jj < jjmax; jj++){
					if(sim_depref[kk][d][jj] == kmax-1){ sim_depref[kk][d][jj] = k; break;}
				}
				if(jj == jjmax) emsg("Sim: EC21");
			}
			else{
				jjmax = sim_depref_d[kk].size();
				for(jj = 0; jj < jjmax; jj++){
					if(sim_depref_d[kk][jj] == d && sim_depref_k[kk][jj] == kmax-1){
						sim_depref_k[kk][jj] = k; break;
					}
				}
				if(jj == jjmax) emsg("Sim: EC21");
			}
    }
  }
  sim_transdeplist[d].pop_back();
}

void siminit()                                   // Initialises variables for use in simuations
{
  long c, p, i, k;
  long tr;

  for(c = 0; c <= ncomp; c++){
    ncomppopnumch.push_back(0);
    comppopnumch.push_back(vector<long>());
    comppopnumchweight.push_back(vector<double>());
  }
  sim_popnum = new double[npopnum];

  for(p = 0; p < npopnum; p++){
    for(i = 0; i < popnumterm[p].size(); i++){
      c = popnumterm[p][i];
      comppopnumch[c].push_back(p);
      comppopnumchweight[c].push_back(popnumtermweight[p][i]);
      ncomppopnumch[c]++;
    }
  }

  sourcetra.resize(nclassval[settimecl]); nsourcetra.resize(nclassval[settimecl]);
  for(k = 0; k < ncompleave[NOTALIVE]; k++){
    tr = compleave[NOTALIVE][k];
    c = tra[tr].cf;
    sourcetra[compval[c][settimecl]].push_back(tr);
  }
  for(k = 0; k < nclassval[settimecl]; k++) nsourcetra[k] = sourcetra[k].size();
}

void Chain::sim_check(double R, double t)   // Checks that the simulation is being performed correctly
{
  long c, ncompind[ncomp], p, d, j, jj, eq, ts, k;
  long i, tr;
  double num, dd, sum;

  if(indev_sim.size() != nindtot_sim) emsg("Sim: EC22");
  if(sim_statind.size() != nindtot_sim) emsg("Sim: EC23");

  // Checks popnum is ok
  for(c = 0; c < ncomp; c++) ncompind[c] = 0;

  for(i = 0; i < nindtot_sim; i++){
    c = sim_statind[i]; if(c != NOTALIVE) ncompind[c]++;
  }

  for(p = 0; p < npopnum; p++){
    num = 0; 
		for(i = 0; i < popnumterm[p].size(); i++){
			num += ncompind[popnumterm[p][i]]*popnumtermweight[p][i];
		}
    dd = num - sim_popnum[p];
    if(dd*dd > tiny) emsg("Sim: EC24");
  }

  // Checks transdeplist ok
  for(d = 0; d < ntransdepeq; d++){
    for(j = 0; j < sim_transdeplist[d].size(); j++){
      i = sim_transdeplist[d][j];
      if(i >= 0){
				if(corona == 0){
					for(jj = 0; jj < sim_depref[i][d].size(); jj++){
						if(sim_depref[i][d][jj] == j) break;
					}
					if(jj == sim_depref[i][d].size()) emsg("Sim: EC25");
				}
				else{
					for(jj = 0; jj < sim_depref_d[i].size(); jj++){
						if(sim_depref_d[i][jj] == d && sim_depref_k[i][jj] == j) break;
					}
					if(jj == sim_depref_d[i].size()) emsg("Sim: EC25b");
				}
      }
    }
  }

  for(d = 0; d < ntransnotdepeq; d++){
    for(j = 0; j < sim_transnotdeplist[d].size(); j++){
      i = sim_transnotdeplist[d][j];
			if(i >= 0){
				if(corona == 0){
					for(jj = 0; jj < sim_notdepref[i][d].size(); jj++){
						if(sim_notdepref[i][d][jj] == j) break;
					}
					if(jj == sim_notdepref[i][d].size()) emsg("Sim: EC26");
				}
				else{
					for(jj = 0; jj < sim_notdepref_d[i].size(); jj++){
						if(sim_notdepref_d[i][jj] == d && sim_notdepref_k[i][jj] == j) break;
					}
					if(jj == sim_notdepref_d[i].size()) emsg("Sim: EC26b");
				}
      }
    }
  }
	
  // Checks rates are OK
  for(d = 0; d < ntransdepeq; d++){
    dd = sim_Rtransdep[d] - ratecalc(transdepeq[d],sim_popnum,param);
    if(dd*dd > tiny) emsg("Sim: EC27");
  }

  sim_Rtransnotdep.resize(ntransnotdepeq);
  for(d = 0; d < ntransnotdepeq; d++){
    dd = sim_Rtransnotdep[d] - ratecalc(transnotdepeq[d],sim_popnum,param);
    if(dd*dd > tiny) emsg("Sim: EC28");
  }

  // Checks the overall rate is OK
  sum = 0;
  for(i = 0; i < sim_statind.size(); i++){
    c = sim_statind[i];
    if(c != NOTALIVE){
      for(j = 0; j < ncompleave[c]; j++){
        tr = compleave[c][j];
        if(tra[tr].type == EXP_TR){
          eq = tra[tr].eq; if(eq < 0) emsg("Sim: EC29");
          sum += ratecalc(eq,sim_popnum,param);
        }
      }
    }
  }

  ts = 0; while(ts < settime.size() && t >= settime[ts]) ts++; 
	if(ts >= nclassval[settimecl]) emsg("Sim: EC30");
	
  for(k = 0; k < nsourcetra[ts]; k++){           // Adds in source transitions
    tr = sourcetra[ts][k]; eq = tra[tr].eq; sum += ratecalc(eq,sim_popnum,param);
  }

  dd = R-sum;
  if(dd*dd > tiny) emsg("Sim: EC31");
}
