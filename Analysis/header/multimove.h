 // This file considers making multiple event time proposals

void Chain::multimove()   // Each transition type is considered in turn and multiple event move are made
{
	double facup = 1.005, facdo = 0.995, facup2 = 1.05, facdo2 = 0.95;
	
	long i, e, trr, trm, c, eq, loop, loopmax, ddti, ddtf, ddi, ddf;
	double t, tst, gr, al, probif, probfi, Ltoti, Ltotf, pr, jump;
	vector <double> timest, timeoldst, grst;
	vector <long> ist, est, mv, ddist, ddfst;
	vector <vector <long> > movetri, movetre;
	
	if(discon == 0) emsg("multimove: EC1");
	
	movetri.clear(); movetre.clear();              // Finds out which events need moving
	movetri.resize(ntra); movetre.resize(ntra);
	
	for(i = 0; i < nind; i++){
		for(e = 0; e < nindev[i]; e++){
			trm = tra[indev[i][e].tr].tramm;
			if(trm >= 0){
				movetri[trm].push_back(i);
				movetre[trm].push_back(e);
			}
		}
	}
	
	for(trm = 0; trm < ntramm; trm++){
		loopmax = movetri[trm].size();
		if(loopmax > 0){
			jump = jump_multimove[trm]; pr = multimovepr[trm];
			
			timest.clear(); timeoldst.clear(); ist.clear(); est.clear(); mv.clear(); grst.clear(); ddist.clear(); ddfst.clear();
	
			Ltoti = L();
			probif = 0; probfi = 0;
		
			for(loop = 0; loop < loopmax; loop++){
				if(ran() < pr){
					i = movetri[trm][loop];
					e = movetre[trm][loop];
					c = indinit[i];

					tst = indev[i][e].t;
					if(checkiffixed(i,tst) == 0){ 
						t = normal(tst,jump);
			
						ntr_multimove[trm]++;
						if(t < 0 || t > tmax || 
							(e > 0 && t <= indev[i][e-1].t) || (e < nindev[i]-1 && t >= indev[i][e+1].t)){
							nfa_multimove[trm]++;
							if(samp < burnin) jump *= facdo;
						}
						else{
							trr = indev[i][e].tr;
							gr = movecalcgrad(trr,tst);  // TO DO this assumes that the gradient is time invarient...
						
							ist.push_back(i); est.push_back(e); timest.push_back(t); timeoldst.push_back(tst); grst.push_back(gr);
						
							al = exp(gr*(t-tst));
							
							eq = tra[trr].eq;
							if(transdep[eq] == 1){       // TODO assumes both initial and final states are dependent
								ddti = DDconv[long(tst*DDfac)]; if(DDt[ddti+1] < tst) ddti++;
								if(tramm[trm].tra[DDsettime[ddti]] != trr) emsg("multimove: EC80");
								ddi = transdepref[eq]*nDD + ddti;  
							
								ddtf = DDconv[long(t*DDfac)]; if(DDt[ddtf+1] < t) ddtf++;
								trr = tramm[trm].tra[DDsettime[ddtf]];
								eq = tra[trr].eq;
								ddf = transdepref[eq]*nDD + ddtf; 

								if(DDcalc[ddf] == 0 || DDcalc[ddi] == 0) emsg("multimove: EC81");
								
								al *= depeq_disc_evval[ddf]/depeq_disc_evval[ddi];
								ddist.push_back(ddi); ddfst.push_back(ddf);
							}
							else{
								ddist.push_back(-1); ddfst.push_back(-1);
							}
							if(al > 1) al = 1;
							
							if(ran() < al){
								nac_multimove[trm]++;
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
				i = ist[loop]; e = est[loop]; t = timest[loop]; tst = timeoldst[loop]; gr = grst[loop]; ddi = ddist[loop]; ddf = ddfst[loop];		
				if(mv[loop] == 0){
					al = exp(gr*(t-tst));
					if(ddi >= 0) al *= depeq_disc_evval[ddf]/depeq_disc_evval[ddi];
					if(al > 1) al = 1;
					probfi += log(1-al);
				}
				else{
					al = exp(gr*(tst-t));
					if(ddi >= 0) al *= depeq_disc_evval[ddi]/depeq_disc_evval[ddf];
					if(al > 1) al = 1;
					probfi += log(al);
				}
			}
			Ltotf = L();
						
			al = exp(Ltotf - Ltoti + probfi - probif);
			
			ntr_multimovetot[trm]++;
			if(ran() < al){
				nac_multimovetot[trm]++;
				if(samp < burnin){ pr *= facup2; if(pr > 0.3) pr = 0.3;}
			}
			else{
				if(samp < burnin) pr *= facdo2;
					
				multich = -1; multisec();
				for(loop = 0; loop < ist.size(); loop++){   // Moves events back
					if(mv[loop] == 1) eventmove(ist[loop],est[loop],timeoldst[loop]);
				}
			}
			
			jump_multimove[trm] = jump; multimovepr[trm] = pr;
			multisecend();
		}
	}
}

/*
double Chain::multial(long tr, double t, double tst)     // Predicts the acceptance rate for a move
{
	long ci, cf, k, kmax, ref, j, jmax;
	double grad;
	NDEQCH ndeqch;
	
	ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
		
	ci = tra[tr].ci; cf = tra[tr].cf;
	
	grad = 0;
	kmax = transchref[ci][cf].size();   
	for(k = 0; k < kmax; k++){
    ref = transchref[ci][cf][k];
		
		jmax = transnotdepeqch[ref].size();
		for(j = 0; j < jmax; j++){        
			ndeqch = transnotdepeqch[ref][j];
			grad += transnotdepeq_val[ndeqch.d]*ndeqch.n;
		}
	}
	
	al = exp(grad*(t-tst));
	 
	eq = tra[tr].eq; 
	if(transdep[eq] == 1){
		
		al *= transrate(eq,tst)/transrate(eq,t);  
	
		if(DDcalc[d*nDD + ddti] == 0)
	}	

	if(al > 1) al = 1;
	return al;
}
*/

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
	long tr, tr2, ci, cf, s;
	
	ntramm = 0;                            // Multimove transitions collect together all all those from the top three classifications
	for(tr = 0; tr < ntra; tr++){	
		ci = tra[tr].ci; cf =  tra[tr].cf;
		
		if(tra[tr].cl >= agecl || tra[tr].cl < 0 || ci == NOTALIVE || cf == NOTALIVE) tr2 = -1;
		else{		 
			ci = ci%classmult[nclass-3];
			cf = cf%classmult[nclass-3];
			
			for(tr2 = 0; tr2 < ntramm; tr2++) if(ci == tramm[tr2].ci && cf == tramm[tr2].cf) break;
			if(tr2 == ntramm){
				TRAMM tram;
				tram.ci = ci; tram.cf = cf; for(s = 0; s <= nsettime; s++) tram.tra.push_back(-1);
				tramm.push_back(tram); 
				ntramm++; 
			}	
			tramm[tr2].tra[compval[tra[tr].ci][settimecl]] = tr;
		}
		tra[tr].tramm = tr2;
	} 
	
	if(1 == 0){
		for(tr = 0; tr < ntramm; tr++){
			cout << compname[tramm[tr].ci] << "->" << compname[tramm[tr].cf] << "  tramm\n";
			for(s = 0; s <= nsettime; s++) cout << tramm[tr].tra[s] << ",";
      cout << "tra\n";
		}
	}
	
	for(tr = 0; tr < ntramm; tr++){
		for(s = 0; s <= nsettime; s++) if(tramm[tr].tra[s] == -1) emsg("multimove: EC77");
	}
	
	for(tr = 0; tr < ntramm; tr++){
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
		while(e < nev && evnew[e].t < t) e++;
    if(e == nev) emsg("Multimove: EC1");
		
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
		while(e < nev && evnew[e].t < t) e++;
    if(evnew[e].t != t) emsg("Multimove: EC2");
		
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
