 // This file considers making multiple event time proposals

void Chain::multimove()   // Each transition type is considered in turn and multiple event move are made
{
	double facup = 1.005, facdo = 0.995, facup2 = 1.05, facdo2 = 0.95;
	
	long ci, cf, k, i, e, tr, c, eq, loop, loopmax;
	double t, tst, List, gr, al, probif, probfi, Ltoti, Ltotf, pr, jump;
	vector <double> timest, timeoldst;
	vector <long> ist, est, mv;
	vector <vector <long> > movetri, movetre;
	
	movetri.clear(); movetre.clear();              // Finds out which events need moving
	movetri.resize(ntra); movetre.resize(ntra);
	
	for(i = 0; i < nind; i++){
		for(e = 0; e < nindev[i]; e++){
			tr = indev[i][e].tr;
			if(indev[i][e].t != 0 && indev[i][e].t != tmax && tra[tr].cl < agecl){
				movetri[tr].push_back(i);
				movetre[tr].push_back(e);
			}
		}
	}
	
	for(tr = 0; tr < ntra; tr++){
		loopmax = movetri[tr].size();
		if(loopmax > 0){
			jump = jump_multimove[tr]; pr = multimovepr[tr];
			
			timest.clear(); timeoldst.clear(); ist.clear(); est.clear(); mv.clear();
	
			Ltoti = L();
			probif = 0; probfi = 0;
			
			for(loop = 0; loop < loopmax; loop++){
				if(ran() < pr){
					i = movetri[tr][loop];
					e = movetre[tr][loop];
					c = indinit[i];

					tst = indev[i][e].t;
					if(checkiffixed(i,tst) == 0){ 
						t = normal(tst,jump);
			
						ntr_multimove[tr]++;
						if(t < 0 || t > tmax || 
							(e > 0 && t <= indev[i][e-1].t) || (e < nindev[i]-1 && t >= indev[i][e+1].t)){
							nfa_multimove[tr]++;
							if(samp < burnin) jump *= facdo;
						}
						else{
							al = exp(movecalcgrad(tr,tst)*(t-tst));
						
							eq = tra[tr].eq; if(transdep[eq] == 1) al *= transrate(eq,t)/transrate(eq,tst);  
							
							if(al > 1) al = 1;
							
							ist.push_back(i); est.push_back(e); timest.push_back(t); timeoldst.push_back(tst);
							
							if(ran() < al){
								nac_multimove[tr]++;
								if(samp < burnin) jump *= facup;
				
								probif += log(al);
								mv.push_back(1);
							}
							else{
								if(samp < burnin) jump *= facdo;
				
								probif += log(1-al);
								mv.push_back(0);
							}
						}
					}
				}
			}
			
			multich = 1;
			for(loop = 0; loop < ist.size(); loop++){   // Actually moves the events
				i = ist[loop]; e = est[loop]; t = timest[loop];
				if(mv[loop] == 1) eventmove(i,e,t);
			}
			multisec();
			
			for(loop = 0; loop < ist.size(); loop++){
				i = ist[loop]; e = est[loop]; t = timest[loop]; tst = timeoldst[loop];
				if(mv[loop] == 0){
					al = exp(movecalcgrad(tr,tst)*(t-tst));
					eq = tra[tr].eq; if(transdep[eq] == 1) al *= transrate(eq,t)/transrate(eq,tst);  
					if(al > 1) al = 1;
					probfi += log(1-al);
				}
				else{
					al = exp(movecalcgrad(tr,t)*(tst-t));
					eq = tra[tr].eq; if(transdep[eq] == 1) al *= transrate(eq,tst)/transrate(eq,t);  
						
					if(al > 1) al = 1;
					probfi += log(al);
				}
			}
			Ltotf = L();
			
			al = exp(Ltotf - Ltoti + probfi - probif);
		
			ntr_multimovetot[tr]++;
			if(ran() < al){
				nac_multimovetot[tr]++;
				if(samp < burnin){ pr *= facup2; if(pr > 0.3) pr = 0.3;}
			}
			else{
				if(samp < burnin) pr *= facdo2;
					
				multich = -1; multisec();
				for(loop = 0; loop < ist.size(); loop++){   // Moves events back
					if(mv[loop] == 1) eventmove(ist[loop],est[loop],timeoldst[loop]);
				}
			}
			
			jump_multimove[tr] = jump; multimovepr[tr] = pr;
			multisecend();
		}
	}
}

double Chain::movecalcgrad(long tr, double t)            // Works out the gradient in the loglikelihood
{
	long e, ee, dn, d, j, jmax, k, kmax, ref, ci, cf;
	double tt, grad = 0;
	NDEQCH ndeqch;
	
	ci = tra[tr].ci; cf = tra[tr].cf;
	
	kmax = transchref[ci][cf].size();   
	for(k = 0; k < kmax; k++){
    ref = transchref[ci][cf][k];
		
		// Assumes that the contribution coming from infection event is small (because few people become inf)
		if(1 == 0){
			EQCH eqch;
		
			jmax = transdepeqch[ref].size();   // This contribution is quiet small
			for(j = 0; j < jmax; j++){ 
				eqch = transdepeqch[ref][j];
				dn = eqch.n;
				if(dn != 0){
					d = eqch.d; 
					e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];      // finds where we are on the timeline
					do{
						ee = depeq_evnext[e]; tt = depeq_evt[ee]; if(tt > t) break;
						e = ee;
					}while(1 == 1);	
					
					grad += dn*depeq_evval[e];
				}
			}
		}
		
		jmax = transnotdepeqch[ref].size();
		for(j = 0; j < jmax; j++){        
			ndeqch = transnotdepeqch[ref][j];
			grad += transnotdepeq_val[ndeqch.d]*ndeqch.n;
		}
	}
	
	return grad;
}

void Chain::multimoveinit()              // Initialises multimove proposals
{
	long tr;
	
	for(tr = 0; tr < ntra; tr++){
		jump_multimove.push_back(1); multimovepr.push_back(0.1);
		ntr_multimove.push_back(0); nac_multimove.push_back(0);	nfa_multimove.push_back(0);
		ntr_multimovetot.push_back(0); nac_multimovetot.push_back(0);
	}
}

void addsettime()                       // Adds settimes to a constructed event sequence evnew
{
	long e, s, val, ci, cf, tr, ee, nev;
	double t;

	val = classmult[settimecl];
	e = 1; nev = evnew.size();
	for(s = 0; s < nsettime; s++){
		t = settime[s];
		while(e < nev && evnew[e].t < t) e++; if(e == nev) emsg("Multimove: EC1");
		
		ci = tra[evnew[e].tr].ci;
		cf = ci + val;
		EV ev; ev.tr = compiftra[ci][cf]; ev.t = t; evnew.insert(evnew.begin()+e,ev);	
		e++; nev++;
		for(ee = e; ee < nev; ee++){
			tr = evnew[ee].tr;
			ci = tra[tr].ci+val; cf = tra[tr].cf; 
			if(cf != NOTALIVE) tr = compiftra[ci][cf+val];
			else tr = traend+ci;
			evnew[ee].tr = tr;
		}
	}
}

void remsettime()                      // Removes settimes to a constructed event sequence evnew
{
	long e, s, val, ci, cf, tr, ee, nev;
	double t;
	
	val = classmult[settimecl];
	
	e = 1; nev = evnew.size();
	for(s = 0; s < nsettime; s++){
		t = settime[s];
		while(e < nev && evnew[e].t < t) e++; if(evnew[e].t != t) emsg("Multimove: EC2");
		
		evnew.erase(evnew.begin()+e);	
		for(ee = e; ee < nev; ee++){
			tr = evnew[ee].tr;
			ci = tra[tr].ci-val; cf = tra[tr].cf; 
			if(cf != NOTALIVE) tr = compiftra[ci][cf-val];
			else tr = traend+ci;
			evnew[ee].tr = tr;
		}
	}
}
