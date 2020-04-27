// The file generates output files

vector <long> traoutref;
vector <string> traoutname;
long ntraout;
vector <long> ntraev;
long evtot, nindtotst;

void traceinit()                                        // Initialises the trace plot
{
  long p, pp, c, cl, i, f, j, tr, d;

  ntraout = 0;
  traoutref.resize(ntra);

  for(tr = 0; tr < ntra; tr++){
    traoutref[tr] = -1;

    cl = tra[tr].cl;
    if(cl >= 0){
      i = tra[tr].i; f = tra[tr].f;
      if(i >= 0 && i < nclassval[cl] && f >= 0 && f < nclassval[cl]){
        stringstream ss; ss << classval[cl][i] << " → " << classval[cl][f];
        j = 0; while(j < ntraout && traoutname[j] != ss.str()) j++;
        if(j == ntraout){ traoutname.push_back(ss.str()); ntraout++;}
        traoutref[tr] = j;
      }
    }
  }
  ntraev.resize(ntraout);

  if(noout == 1) return;

  if(simon == 1) bout << "c|"; else bout << "1|";

  for(pp = 0; pp < nparam; pp++){
    p = paramorder[pp];
    bout << paramclass[p] << "|" << paramname[p] << "|Transition rate parameter " << paramname[p] << "|";
  }

  if(simon == 0){
    for(c = 0; c < ncompswa; c++){
      bout << "Init. Prob.|";
      for(cl = 0; cl < nclass-2; cl++){ 
				if(!(cl == agecl && nclassval[agecl] == 1)){
					if(cl != 0) bout << ","; 
					bout << classval[cl][compval[c][cl]];
				}
			}
      bout << "|";
      bout << "Probability of initially being in ";
      for(cl = 0; cl < nclass-2; cl++){ 
				if(!(cl == agecl && nclassval[agecl] == 1)){ 
					if(cl != 0) bout << ","; bout << classval[cl][compval[c][cl]];
				}
			}
      bout << "|";
    }
  }

  for(j = 0; j < ntraout; j++){
		bout << "Trans.|" << traoutname[j] << "|Number of " <<  traoutname[j] << " transitions|";
	}
	
  for(d = 0; d < nderive; d++){
    if(neq_popnum[derive[d]] == 0){
			bout << "Der.|" << derivename[d] << "|Derived quantity: " << derivename[d] << "|";
		}
  }

  if(simon == 0){
    if(checkon == 1){
      bout << "Misc.|Lev|Event likelihood|";
      bout << "Misc.|Lnm|Non Markovian event likelihood|";
      bout << "Misc.|Lob|Observation probability|";
      bout << "Misc.|Lcap|Capture probability|";
      bout << "Misc.|Levcap|Event capture probability|";
      bout << "Misc.|Linit|Initial condition probability|";
      bout << "Misc.|Lpop|Population probability|";
      bout << "Misc.|Lpop|Derived probability|";
      bout << "Misc.|Pr|Prior probability|";
      bout << "Misc.|Nerr|The number of observation|";
    }
    else{
      bout << "Misc.|log(π(θ,ξ¦D))|The log of the posterior probability|";
      bout << "Misc.|log(π(D¦ξ,θ))|The log of the observation probability|";
      bout << "Misc.|log(π(ξ¦θ))|The log of latent process likelihood|";
      bout << "Misc.|log(π(θ))|The log of the prior probability|";
    }
  }
  bout << "Misc.|# Ind.|Total number of individuals|";
  bout << "Misc.|# Transitions|Total number transitions|";
  bout << "\n";
  bout.flush();
}

void Chain::traceplot()                                   // Plots traces for the parameter values
{
  long p, pp, cl, k, d, c, e;
  long i, j, tr, evtotst, nindtotst;
  double Ltot, Lobtot;

  if(noout == 1) return;

  timeprop[TRACEPLOT] -= clock();

  evtotst = 0; 
  if(simon == 1){ nindtotst = nindtot_sim; for(i = 0; i < nindtot_sim; i++) evtotst += nindev_sim[i];}
  else{ nindtotst = nindtot; for(i = 0; i < nindtot; i++) evtotst += nindev[i];}

  if(simon == 1) bout << "d|"; else bout << "0|";

  for(pp = 0; pp < nparam; pp++){ p = paramorder[pp]; bout << param[p] << "|";}

  if(simon == 0){
    for(c = 0; c < ncompswa; c++) bout << probinit[c] << "|";
  }

  for(j = 0; j < ntraout; j++) bout << ntraev[j] << "|";

  for(d = 0; d < nderive; d++){
		if(neq_popnum[derive[d]] == 0) bout << calculatenotdep(derive[d],param) << "|";
	}

  if(simon == 0){
    Ltot = Lir-Liexp + Linm + Liinit;
    Lobtot = Lob_st+Lob_pd +Lev_pd + Lpop+Lder;
    if(checkon == 1){
			bout << Lir-Liexp << "|" << Linm << "|" << Lob_st << "|" << Lob_pd << "|" << Lev_pd 
						<< "|" << Liinit << "|"<<  Lpop << "|"<< Lder << "|" << Lpri << "|" << nwrong << "|";
		}
    else bout << Ltot+Lobtot+Lpri << "|" << Lobtot << "|" << Ltot << "|" << Lpri << "|";
  }
  bout << nindtotst << "|" << evtotst << "|" ;

  bout << "\n";
  bout.flush();

  timeprop[TRACEPLOT] += clock();
}

void oe(string st, vector<EV> ev)                 // Outputs an event sequence
{
  long e, ci, cf;
  long tr;
  double t;

  timeprop[EVENTPLOT] -= clock();

  if(ev.size() == 0) cout << "NA ";
  else{
    if(1 == 1){
      for(e = 0; e < ev.size(); e++){
        tr = ev[e].tr; t = ev[e].t;
        if(tra[tr].ci == NOTALIVE) cout << "NA";
        cout << " -> "; outcomp(tra[tr].cf); cout << " " << t << ",       "; 
        //outcomp(tra[tr].ci); cout << " -> "; outcomp(tra[tr].cf); cout << " " << t << ",       "; 
      }
    }
    else{
      outcomp(tra[ev[0].tr].ci); cout << ": ";
      for(e = 0; e < ev.size(); e++){
        tr = ev[e].tr; t = ev[e].t;
        if(tra[tr].type < 0){
          cout << "END";
        }
        else{
          ci = tra[tr].ci; cf = tra[tr].cf;
          cout << " -> ";
          if(ci == NOTALIVE || cf == NOTALIVE) outcomp(cf);
          else cout << classval[tra[tr].cl][tra[tr].f];
        }
        cout << " " << t << ",  ";
      }
    }
  }
  cout << st << "\n";

  timeprop[EVENTPLOT] += clock();
}

void oprop()           // Outputs a given proposal
{
  cout << "\n";
  oe("BEF",evrev);
  oe("AFT",evnew);
}

void outputmodel()      // Outputs the model (used in debugging)
{
  long cl, tr, d, p, j, c, eq, ob, cap, f;

  ofstream model("GP/model.txt");

  for(p = 0; p < npopm; p++){
    model << p << " t: " <<  popmt[p] << "   c: ";
    for(j = 0; j < popmcomp[p].size(); j++) model << popmcomp[p][j] << ",";
    model << "\n";
  }

  model << "Trasitions:\n";
  for(tr = 0; tr < ntra; tr++){
    cl = tra[tr].cl;
    if(cl < 0){
      c = tra[tr].ci; 
			if(c == NOTALIVE) model << "NA"; 
			else{ 
				for(cl = 0; cl < nclass; cl++){
					model << classval[cl][compval[c][cl]]; if(cl < nclass-1) model << ",";
				}
			}
      model << " -> ";
      c = tra[tr].cf;
			if(c == NOTALIVE) model << "NA"; 
			else{ 
				for(cl = 0; cl < nclass; cl++){
					model << classval[cl][compval[c][cl]]; if(cl < nclass-1) model << ",";
				}
			}
    }
    else model << classname[cl] << ": " << classval[cl][tra[tr].i] << " -> " << classval[cl][tra[tr].f];
    model << "     ";
    eq = tra[tr].eq; if(eq >= 0) model << eqnstr[tra[tr].eq];
    model << "  Capev: " << tra[tr].capev << "   Like: " << tra[tr].like << "\n";
  }

  for(d = 0; d < ntransdepeq; d++) model << eqnstr[transdepeq[d]] << " depeqn\n";
  for(d = 0; d < ntransnotdepeq; d++) model << eqnstr[transnotdepeq[d]] << " notdepeqn\n";

  for(p = 0; p < paramname.size(); p++){
    model << paramname[p] << "\n";
  }

  for(p = 0; p < popnumname.size(); p++){
    model << popnumname[p] << ": ";
    for(j = 0; j < popnumterm[p].size(); j++){
			model << popnumterm[p][j] << " " << popnumtermweight[p][j] << ",  ";
		}
    model << "\n"; 
  }

  for(cl = 0; cl < nclass; cl++){
    model << classname[cl] << ":\n";
    for(j = 0; j < nclassval[cl]; j++){
      model << classval[cl][j] << "\n";
    }
  }

  for(cap = 0; cap < ncap; cap++){
    model << capt[cap] << " " << capname[cap] << "; ";
    model << " cap\n\n";
  }
  model << "\n";

  for(ob = 0; ob < nobs; ob++){
    model << obst[ob] << " " << indid[obsi[ob]] << " ";
    if(obscap[ob] >= 0)  model << capname[obscap[ob]] << "; ";
    //for(c = 0; c < ncomp; c++) if(obsprobeqn[ob][c] != -1) model << compname[c] << ",";
    model << " obs\n\n";
  }

  for(eq = 0; eq < eqnstr.size(); eq++){
    model << "Eqn: " << eqnstr[eq] << "\n";
  }

  model << "Prior:\n";
  for(p = 0; p < priortype.size(); p++){
    model << p << " ";
    switch(priortype[p]){
      case FLAT: model << "FLAT" << priorminval[p] << " " <<  priormaxval[p] << "\n"; break;
      case GAMMA: model << "GAMMA" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case NORMAL: model << "NORMAL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case LOGNORMAL: 
				model << "LOGNORMAL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; 
				break;
      case EXPO: model << "EXPO" << eqnstr[prioreq1[p]] << "\n"; break;
      case BETA: model << "BETA" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case WEIBULL: model << "WEIBULL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;

      default: model << "error prior type\n"; break;
    }
  }

  model << "\n";
  for(f = 0; f <= ncapevtrange; f++) model << capevtrange[f] << ","; model << "capevtrange\n";
  model << "\n";

  model << "fixfl: " << fixfl << "\n";
  model << "sinkfl: " << sinkfl << "\n";
  model << "sourcefl: " << sourcefl << "\n";
  model << "gammafl: " << gammafl << "\n";
  model << "nmfl: " << nmfl << "\n";
  model << "capfl: " << capfl << "\n";
  model << "popfl: " << popfl << "\n";
  model << "tbirthfl: " << tbirthfl << "\n";
  model << "pairfl: " << pairfl << "\n";
  model << "singeventfl: " << singeventfl << "\n";
  model << "twothreefl: " << twothreefl << "\n";
  model << "gapfl: " << gapfl << "\n";
 
  model << "\n";
}

void outcomp(long c)                // Outputs the state of a compartment
{
  long cl;
  if(c == NOTALIVE) cout << "NA";
  else cout << compname[c];
}

void Chain::eventplot()            // Outputs an events sample
{
  long cl, k, loop, tr, cf, i, j, e;
  double fac;
  vector <long> list;
  vector <double> listt;

	if(corona == 1){
		if(nderive > 0){ 
			sim(tmax); 
			deriveplotcalc(); derivepl();
		}
		return;
	}
	
  for(loop = 0; loop < 2; loop++){
    if(simon == 1){ sim(0); bout << "b|";}
    else{
      switch(loop){
        case 0: sim(tmax); bout << "5|"; break;
        case 1: sim(0); bout << "a|"; break;
      }
    }

    fac = (tmaxactual-tminactual)/tmax2; 
    if(fac > 1.01 || fac < 0.99) emsg("Output: EC1");

		if(corona == 0){
			if(noout == 0){
				bout << nindtot_sim << "|";
				for(i = 0; i < nindtot_sim; i++){
					list.clear(); listt.clear();
					for(k = 0; k < nindev_sim[i]; k++){
						tr = indev_sim[i][k].tr;
						if(tra[tr].type >= 0){
							cf = tra[tr].cf;
							if(cf == ncomp) cf = classmult[capevcl]; else cf = cf%classmult[capevcl];
							list.push_back(cf); listt.push_back(tminactual + fac*indev_sim[i][k].t);
						}
					}
					bout << indtbirth_sim[i] << "|";
					bout << list.size() << "|"; 
					for(k = 0; k < list.size(); k++) bout << list[k] << "|" << listt[k] << "|";
				}
				bout << "\n";
			}

			if(noout == 2){
				for(i = 0; i < nindtot_sim; i++){
					bout << "ind" << i << ": ";
					for(k = 0; k < nindev_sim[i]; k++){
						tr = indev_sim[i][k].tr;
						cf = tra[tr].cf;  bout.flush();
						if(cf == ncomp) bout << classmult[capevcl];
						else bout << tra[tr].cf%classmult[capevcl];
						bout << " " << tminactual + fac*indev_sim[i][k].t << ",  ";
					}
					bout << "\n";
				}
			}
		}
		
    if(noout == 0 && loop == 0 && nderive > 0){ deriveplotcalc(); derivepl();}

    if(loop == 0){
      nindtotst = nindtot_sim;

      evtot = 0;
      for(j = 0; j < ntraout; j++) ntraev[j] = 0;

      for(i = 0; i < nindtot_sim; i++){
        for(e = 0; e < nindev_sim[i]; e++){
          tr = indev_sim[i][e].tr;
          if(traoutref[tr] >= 0) ntraev[traoutref[tr]]++;
          if(indev_sim[i][e].t != 0 && indev_sim[i][e].t != tmax && tra[tr].cl < nclass-3) evtot++;
        }
      } 
    }
    if(simon == 1) break;
  }
  bout.flush();
}

void Chain::diagnosticsfile()       // Outputs a file giving MCMC diagnostic information
{
	double f, nf, fmin, fmax, v, p;
	long tr;
	
	stringstream sdi; sdi << root << "/diagnostic" << simnum << ".txt";
	ofstream dia(sdi.str().c_str());

	for(p = 0; p < nparam; p++) if(ntr_param[p] >= 1) break;
  if(p < nparam){
    dia << "RANDOM WALK PARAMETER PROPOSALS\n";
    for(p = 0; p < nparam; p++){
      if(ntr_param[p] >= 1) dia << paramname[p] << " Acceptance probability: " << nac_param[p]/ntr_param[p] << "  Jumping size: " << jump_param[p] << "\n";
    }
  }

  for(p = 0; p < nparam; p++) if(ntr_paramnorm[p] >= 1) break;
  if(p < nparam){
		dia << "\n";
    dia << "NORMAL FITTED PARAMETER PROPOSAL\n";
    for(p = 0; p < nparam; p++){
      if(ntr_paramnorm[p] >= 1){
				dia << paramname[p] << " Acceptance probability: " << nac_paramnorm[p]/ntr_paramnorm[p] << "\n";
			}
    }
  }
	dia << "\n";
	
  dia << "MULTIMOVE\n";
	dia << "ac:";
	f =0; nf = 0; fmin = 100; fmax = 0;
	for(tr = 0; tr < ntra; tr++){
		if(ntr_multimove[tr] > 0){ 
			v = nac_multimove[tr]/ntr_multimove[tr]; 
			f += v; nf++; if(v > fmax) fmax = v; if(v < fmin) fmin = v;
		}
	}
	dia << f/nf << " (" << fmin << "-" << fmax << "),  ";
	 	
	dia << "fa:";
	f =0; nf = 0; fmin = 100; fmax = 0;
	for(tr = 0; tr < ntra; tr++){
		if(ntr_multimove[tr] > 0){ 
			v = nfa_multimove[tr]/ntr_multimove[tr]; 
			f += v; nf++; if(v > fmax) fmax = v; if(v < fmin) fmin = v;
		}
	}
	dia << f/nf << " (" << fmin << "-" << fmax << "),  ";
	dia << "\n";
	
	dia << "actot:";
	f =0; nf = 0; fmin = 100; fmax = 0;
	for(tr = 0; tr < ntra; tr++){
		if(ntr_multimove[tr] > 0){ 
			v = nac_multimovetot[tr]/ntr_multimovetot[tr]; 
			f += v; nf++; if(v > fmax) fmax = v; if(v < fmin) fmin = v;
		}
	}
	dia << f/nf << " (" << fmin << "-" << fmax << "),  ";
	 	
	dia << "jump:";
	f =0; nf = 0; fmin = 100; fmax = 0;
	for(tr = 0; tr < ntra; tr++){
		if(ntr_multimove[tr] > 0){ 
			v = jump_multimove[tr]; 
			f += v; nf++; if(v > fmax) fmax = v; if(v < fmin) fmin = v;
		}
	}
	dia << f/nf << " (" << fmin << "-" << fmax << "),  ";
	 	
	dia << "pr:";
	f =0; nf = 0; fmin = 100; fmax = 0;
	for(tr = 0; tr < ntra; tr++){
		if(ntr_multimove[tr] > 0){ 
			v = multimovepr[tr]; 
			f += v; nf++; if(v > fmax) fmax = v; if(v < fmin) fmin = v;
		}
	}
	dia << f/nf << " (" << fmin << "-" << fmax << ")\n";
	dia << "\n";
	
	dia << "ADDREM MULTINF\n" << nac_multinfadd/ntr_multinfadd << " rem ac " 
			<< nac_multinfrem/ntr_multinfrem <<"\n";

	dia << "\n";
	
	dia << "CPU TIME FOR PROPOSALS\n";
	dia << "Sec: " << long(100*timeprop[PARAM_PROP]/totaltime) << "%\n";
  dia << "Sampprob: " << long(100*timeprop[SAMPPROB]/totaltime) << "%\n";
  dia << "Paramameters: " << long(100*timeprop[PARAM_PROP]/totaltime) << "%\n";
	dia << "Add/rem inf: " << long(100*timeprop[ADDREMINF_PROP]/totaltime) << "%\n";
	dia << "Multimove: " << long(100*timeprop[MULTIMOVE_PROP]/totaltime) << "%\n";
	dia << "Diagnostic checks: " << long(100*timeprop[CHECK]/totaltime) << "%\n";
	dia << "Sim: " << long(100*timeprop[SIMTIME]/totaltime) << "%\n";
  
	dia << "Total CPU time used: " << totaltime/(60.0*CLOCKS_PER_SEC) << " minutes\n";
	dia << "Sample: " << samp << "\n";
}

void Chain::diagnosticschain()        // Outputs the success of MCMC proposals
{
  long p, e, cl, fl1=0, fl2=0, fl3=0, fl4=0, fl5=0, fl6=0, fl7=0, c, tr;
  long i;
  double nf, f, fmin, fmax, fav;

  if(noout == 1) return;

  stringstream ss; ss << setprecision(2);
  
  ss << "6|";
  for(p = 0; p < nparam; p++) if(ntr_param[p] >= 1) break;
  if(p < nparam){
    ss << "RANDOM WALK PARAMETER PROPOSALS|";
    for(p = 0; p < nparam; p++){
      if(ntr_param[p] >= 1){
				ss << paramname[p] << " Acceptance probability: " << nac_param[p]/ntr_param[p] 
						<< "  Jumping size: " << jump_param[p] << "|";
			}
    }
  }

  for(p = 0; p < nparam; p++) if(ntr_paramnorm[p] >= 1) break;
  if(p < nparam){
    ss << "|NORMAL FITTED PARAMETER PROPOSALS|";
    for(p = 0; p < nparam; p++){
      if(ntr_paramnorm[p] >= 1){
				ss << paramname[p] << " Acceptance probability: " << nac_paramnorm[p]/ntr_paramnorm[p] << "|";
			}
    }
  }
  
  for(p = 0; p < nsmooth; p++) if(ntr_smooth[p] >= 1) break;
  if(p < nsmooth){
    ss << "|RANDOM WALK SMOOTHING PRIOR PROPOSALS|";
    for(p = 0; p < nsmooth; p++){
      if(ntr_smooth[p] >= 1){
				ss << smoothname[p] << " Acceptance probability: " << nac_smooth[p]/ntr_smooth[p] 
						<< "  Jumping size: " << jump_smooth[p] << "|";
			}
    }
  }
	
	for(cl = 0; cl <= clall; cl++){
		ss << "|PARTICLE PROPOSALS IN ";
		if(cl == clall) ss << "ALL CLASSIFICATIONS";
		else ss << "'" << classname[cl] << "'";
		ss << "|";

		fmin = 1; fmax = 0; fav = 0; 
		for(i = 0; i < nind; i++){
			f = nac_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
		}
		ss << "Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")|";

		fmin = 1; fmax = 0; fav = 0; 
		for(i = 0; i < nind; i++){
			f = nfa_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
		}
		ss << "Failure: " << fav << " (" << fmin << "-" << fmax << ")|";

		fmin = 1; fmax = 0; fav = 0; 
		for(i = 0; i < nind; i++){ 
			f = nde_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
		}
		ss << "Degenerate: " << fav << " (" << fmin << "-" << fmax << ")|";

		fmin = npartmax; fmax = 0; fav = 0; 
		for(i = 0; i < nind; i++){ 
			f = nindpart[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
		}
		ss << "# Particles used: " << fav << " (" << fmin << "-" << fmax << ")|";
	}
	
	ss << "|ADDREM MULTINF|add ac " << nac_multinfadd/ntr_multinfadd
										<< " rem ac " << nac_multinfrem/ntr_multinfrem <<"|";

  if(ntr_indsimuo >= 1){
    ss << "|SIMULATED PROPOSALS (USED ON UNOBSERVED INDIVIDUALS)|"; 
    ss << "Acceptance probability: " <<  nac_indsimuo/ntr_indsimuo << "|";
    ss << "Failure: " <<  nfa_indsimuo/ntr_indsimuo << "      |";
  }

  for(cl = 0; cl < clall; cl++){
    ss << "|LOCAL CHANGES IN '" <<  classname[cl] << "'|";
    
		if(corona == 0){
			for(i = 0; i < nind; i++) if(ntr_move[i][cl] >= 1) break;
			if(i < nind){
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){
					f = nac_move[i][cl]/ntr_move[i][cl];
					if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "MOVE EVENT   Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")      |";
			}
		}
		
		fmin = 1; fmax = 0; fav = 0; double nfav = 0;
		for(c = 0; c < ncomp; c++){ 
			if(ntr_movecomp[c][cl] > 1){
				f = nac_movecomp[c][cl]/ntr_movecomp[c][cl]; 
				if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f; nfav++;
			}
		}
		
		if(corona == 0){
			ss << "MOVE EVENT   Acceptance probability: " << fav/nfav << " (" << fmin << "-" << fmax << ")      |";

			for(i = 0; i < nind; i++) if(ntr_sing[i][cl] >= 1) break;
			if(i < nind){
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nac_sing[i][cl]/ntr_sing[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "ADD/REMOVE EVENT   Acceptance probability: " << fav << " (" 
						<< fmin << "-" << fmax << ")      ";
		
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nfa_sing[i][cl]/ntr_sing[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "Failure: " << fav << " (" << fmin << "-" << fmax << ")      |";
			}

			for(i = 0; i < nind; i++) if(ntr_twothree[i][cl] >= 1) break;
			if(i < nind){
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nac_twothree[i][cl]/ntr_twothree[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "TWO/THREE EVENT    Acceptance probability: " 
						<< fav << " (" << fmin << "-" << fmax << ")      ";
		
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nfa_twothree[i][cl]/ntr_twothree[i][cl]; 
					if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "Failure: " << fav << " (" << fmin << "-" << fmax << ")      |";
			}

			for(i = 0; i < nind; i++) if(ntr_pair[i][cl] >= 1) break;
			if(i < nind){
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nac_pair[i][cl]/ntr_pair[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "PAIR EVENT   Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")      ";

				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nfa_pair[i][cl]/ntr_pair[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "Failure: " << fav << " (" << fmin << "-" << fmax << ")      |";
			}

			for(i = 0; i < nind; i++) if(ntr_gap[i][cl] >= 1) break;
			if(i < nind){
				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nac_gap[i][cl]/ntr_gap[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "EVENT GAP   Acceptance probability: " << fav << ", " << fmin << "-" << fmax << "      ";

				fmin = 1; fmax = 0; fav = 0; 
				for(i = 0; i < nind; i++){ 
					f = nfa_gap[i][cl]/ntr_gap[i][cl]; 
					if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
				}
				ss << "Failure: " << fav << " (" << fmin << "-" << fmax << ")      |";
			}
		}

		for(i = 0; i < nind; i++){
			if(ntr_tbirth[i] >= 1) fl1 = 1;
			if(ntr_tbirthent[i] >= 1) fl2 = 1;
			if(ntr_tent[i] >= 1) fl3 = 1;
			if(ntr_tlea[i] >= 1) fl4 = 1;
			if(ntr_entswitch[i] >= 1) fl5 = 1;
			if(ntr_birthentswitch[i] >= 1) fl6 = 1;
			if(ntr_leaswitch[i] >= 1) fl7 = 1;
		}
	}
	
  if(fl1+fl2+fl3+fl4+fl5+fl6+fl7 > 0){
    ss << "|LIFESPAN PROPOSALS FOR OBSERVED INDIVIDUALS|";

    if(fl1 == 1){
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_tbirth[i]/ntr_tbirth[i]; if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "BIRTH TIME   Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")      |";	
    }

    if(fl2 == 1){  
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_tbirthent[i]/ntr_tbirthent[i]; 
				if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "BIRTH/ENTER TIME   Acceptance probability: " 
					<< fav << " (" << fmin << "-" << fmax << ")      |";	
    }

    if(fl3 == 1){ 
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_tent[i]/ntr_tent[i]; if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "ENTER TIME   Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")      |";
    }

    if(fl4 == 1){ 
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_tlea[i]/ntr_tlea[i]; if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "LEAVE TIME   Acceptance probability: " << fav << " (" << fmin << "-" << fmax << ")      |";	
    }

    if(fl5 == 1){ 
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_entswitch[i]/ntr_entswitch[i]; 
				if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "ENTER/EXIST SWITCH   Acceptance probability: " 
					<< fav << " (" << fmin << "-" << fmax << ")      |";
    }

    if(fl6 == 1){ 
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_birthentswitch[i]/ntr_birthentswitch[i]; 
				if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "BIRTH/ENTER SWITCH   Acceptance probability: " 
					<< fav << " (" << fmin << "-" << fmax << ")      |";
    }

    if(fl7 == 1){
      fmin = 1; fmax = 0; fav = 0; 
			for(i = 0; i < nind; i++){ 
				f = nac_leaswitch[i]/ntr_leaswitch[i]; 
				if(f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;
			}
      ss << "LEAVE/EXIST SWITCH   Acceptance probability: " 
					<< fav << " (" << fmin << "-" << fmax << ")      |";
    }
  }
  
  if(ntr_tbirthuo >= 1 || ntr_tbirthentuo >= 1 || ntr_tentuo >= 1 || ntr_tleauo >= 1 ||
			ntr_entswitchuo >= 1 || ntr_birthentswitchuo >= 1 || ntr_leaswitchuo >= 1){
    ss << "|LIFESPAN PROPOSALS FOR UNOBSERVED INDIVIDUALS|";
    if(ntr_tbirthuo >= 1){
			ss << "BIRTH TIME  Acceptance probability: " << nac_tbirthuo/ntr_tbirthuo << "      |";
		}
    if(ntr_tbirthentuo >= 1){
			ss << "BIRTH/ENTER TIME  Acceptance probability: " << nac_tbirthentuo/ntr_tbirthentuo << "      |";
		}
    if(ntr_tentuo >= 1){
			ss << "ENTER TIME Acceptance probability: " << nac_tentuo/ntr_tentuo << "      |";
		}
    if(ntr_tleauo >= 1){
			ss << "LEAVE TIME   Acceptance probability: " <<  nac_tleauo/ntr_tleauo << "      |";
		}
    if(ntr_entswitchuo >= 1){
			ss << "ENTER/EXIST SWITCH  Acceptance probability: " << nac_entswitchuo/ntr_entswitchuo << "      |";
		}
    if(ntr_birthentswitchuo >= 1){
			ss << "BIRTH/ENTER SWITCH  Acceptance probability: "
					<<  nac_birthentswitchuo/ntr_birthentswitchuo << "      |";
		}
    if(ntr_leaswitchuo >= 1){
			ss << "LEAVE/EXIST SWITCH  Acceptance probability: "
					<< nac_leaswitchuo/ntr_leaswitchuo << "      |"; 
		}
  }
  
  if(ntr_add >= 1){
    ss << "|ADD INDIVIDUAL   Acceptance probability: " 
				<< nac_add/ntr_add << "     Failure: " << nfa_add/ntr_add << "   |";
    ss << "REMOVE INDIVIDUAL   Acceptance probability: " 
				<< nac_rem/ntr_rem << "    |";
  }
  
  ss << "|CPU TIME FOR PROPOSALS|";
  ss << "Total CPU time used: " << totaltime/(60.0*CLOCKS_PER_SEC) << " minutes|";
	ss << "Sec: " << long(100*timeprop[SEC]/totaltime) << "%|";
	ss << "Sampprob: " << long(100*timeprop[SAMPPROB]/totaltime) << "%|";
		
  ss << "Paramameters: " << long(100*timeprop[PARAM_PROP]/totaltime) << "%|";
  ss << "Particles: " << long(100*timeprop[PART_PROP]/totaltime) << "%|";
  ss << "Simulation of unobserved: " << long(100*timeprop[INDSIMUO_PROP]/totaltime) << "%|";
  ss << "Moving events: " << long(100*timeprop[MOVE_PROP]/totaltime) << "%|";
  ss << "Lifespan proposals: " << long(100*timeprop[LIFE_PROP]/totaltime) << "%|";
  ss << "Adding/removing indivduals: " << long(100*timeprop[ADDREM_PROP]/totaltime) << "%|";
  ss << "Adding/removing single event: " << long(100*timeprop[SINGEVENT_PROP]/totaltime) << "%|";
  ss << "Adding/removing pairs of events: " << long(100*timeprop[PAIR_PROP]/totaltime) << "%|";
  ss << "Local gap changes to events: " << long(100*timeprop[GAP_PROP]/totaltime) << "%|";
  ss << "Diagnostic checks: " << long(100*timeprop[CHECK]/totaltime) << "%|";
  ss << "Time indchange: " << long(100*timeprop[INDCHANGE]/totaltime) << "%|";
  ss << "Time trace: " << long(100*timeprop[TRACEPLOT]/totaltime) << "%|";
  ss << "Time event plot: " << long(100*timeprop[EVENTPLOT]/totaltime) << "%|";
  ss << "Time part sim: " << long(100*timeprop[PARTSIM]/totaltime) << "%|";
  ss << "Time part init: " << long(100*timeprop[PARTINIT]/totaltime) << "%|";
  ss << "Time init: " << long(100*timeprop[INIT]/totaltime) << "%|";
  ss << "Time initback: " << long(100*timeprop[INITBA]/totaltime) << "%|";

  bout << ss.str() << "\n";

  bout.flush();
}

double Chain::Lout()                // Outputs elements of the likelihood
{
   cout << "Lpri: " << Lpri << "   Liinit: " <<  Liinit << "   Lir: " << Lir 
				<< "   Liexp: " << Liexp << "    Linm: " << Linm
        << "   Lev_pd: " << Lev_pd << "   Lob_pd: " << Lob_pd 
				<< "    Lob_st: " << Lob_st << "   Lpop: " << Lpop << "   Lder: " << Lder << "\n";
}

string rep(string s, string s1, string s2)           // Replaces string s1 in s with s2
{
    long pos = s.find(s1);
    if(pos == std::string::npos) return s;
    return s.replace(pos,s1.length(),s2);
}

void tracefileinit()                                 // Initialises the trace file
{
	long r, a, p;
	string st;

	trace << "state\t";
	for(p = 0; p < nparam; p++){
		st = paramname[p];
		st = rep(st,"→","->");
		st = rep(st,"β","Beta");
		st = rep(st,"ρ","Rho");
		st = rep(st,"σ","Sigma");
		st = rep(st,"Φ","Phi");
		st = rep(st,"μ","Mu");
		trace << st << "\t";
	}
	
	for(r = 0; r < nregion; r++){
		for(a = 0; a < nag; a++){
			trace << "inf" << classval[1][r] << "_" << a << "\t";
		}
	}
	trace << "inftot\t" << "Li\t";
	
	for(r = 0; r < nregion; r++) trace << "multacf" << classval[1][r] << "\t";	 
	trace << "end\n";
	
	paramst.resize(nrun);
	for(r = 0; r < nrun; r++) paramst[r].resize(nparam);
}

void Chain::tracefile()               // Outputs the trace file
{
	long p, r, a, tot, num;
	
  trace << samp << "\t";
	for(p = 0; p < nparam; p++) trace << param[p] << "\t"; 
	
	tot = 0;
	for(r = 0; r < nregion; r++){
		for(a = 0; a < nag; a++){
			num = obsinflist[r][a].size()+inflist[r][a].size(); tot += num;
			trace << obsinflist[r][a].size()+inflist[r][a].size() << "\t";
		}
	}
	trace << tot << "\t" << L() << "\t";
	
	for(r = 0; r < nregion; r++) trace << multacf[r] << "\t"; trace << "0\n";
}
