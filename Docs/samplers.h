// This file contains different event samplers
//vector <long> samp_type;

double Chain::simfromfixed(long i)  // Simulated around a single fixed event 
{
	long fi, c, cs, cinit, s, k;
	double tfix, t, dt, r, prob=0, stnext;
	vector <double> tlist;
	vector <long> cst;
	EV evm;
		
	fi = indfixev[i][0];
	tfix = fixev[fi].t;

	cinit = indinit[i];
	
	c = cinit + fixev[fi].traf - compval[cinit][0];  

	t = tfix;
	do{
		r = samprev_r[c]; if(r == 0) break;
		c = samprev_enter[c];
		dt = -log(ran())/r;
		
		t -= dt; 
		if(t < 0){ emsg("neg"); return -large;}
		prob += log(r) - dt*r;

		tlist.push_back(t);
		cst.push_back(c);
	}while(1 == 1);

	evnew.clear();	
	evm.tr = trabeg+cinit; evm.t = 0; evnew.push_back(evm);

	s = 0; stnext = settime[0];
	
	c = cinit;
	for(k = cst.size()-1; k >= 0; k--){
		t = tlist[k];
		while(stnext < t){ 
			cs = c+dcsett;
			evm.tr = compiftra[c][cs]; evm.t = stnext; evnew.push_back(evm); c = cs; s++; 
			stnext = settime[s];
		}
		cs = c + cst[k] - compval[c][0];
		evm.tr = compiftra[c][cs]; evm.t = t; evnew.push_back(evm);
		c = cs;
	}
	
	while(stnext < tfix){ 
		cs = c+dcsett; 
		evm.tr = compiftra[c][cs]; evm.t = stnext; evnew.push_back(evm); c = cs; s++; 
		stnext = settime[s];
	}
		
	cs = c + fixev[fi].traf - compval[c][0];
	evm.tr = compiftra[c][cs]; evm.t = tfix; evnew.push_back(evm);
	c = cs;
	
	prob += simfrom(c,tfix);
	
	return prob;
}

double Chain::probsimfromfixed(vector<EV> &vec, long i)
{
	long e, emax, fi, k, c;
	double t, tfix, r, dt, prob=0;
	vector <double> tlist;
	
	fi = indfixev[i][0];
	tfix = fixev[fi].t;
	
	e = 1; emax = long(vec.size())-1;
	do{
		if(tra[vec[e].tr].cl == 0){
			t = vec[e].t;
			tlist.push_back(t);
			if(t == tfix) break;
		}		
		e++;
	}while(e < emax);
	if(e == emax){
		cout << tfix << " " << nindfixev[i] << " nf\n";
		oe("he",vec);
		emsg("samplers: EC1");
	}
	
	c = tra[vec[e].tr].cf;
	
	for(k = tlist.size()-1; k >= 1; k--){
		r = samprev_r[c]; if(r == 0) emsg("samples: EC2");
		dt = tlist[k] - tlist[k-1];
		prob += log(r) - dt*r;
		c = samprev_enter[c];
	}
		
	prob += probsimfrom(vec,e);

	return prob;
}	

double Chain::simfromSE(long i, double t)  // From a given infection time simulated forwards in time
{
	long c, cs, s;
	EV evm;
	 
	evnew.clear();
	
	c = indinit[i];
	
	evm.tr = trabeg+c; evm.t = 0; evnew.push_back(evm);

	s = 0; 
	while(settime[s] < t){ 
		cs = c+dcsett; evm.tr = compiftra[c][cs]; evm.t = settime[s]; evnew.push_back(evm); c = cs; 
		s++; 
	}
	
	if(ncompleavedep[c] != 1) emsg("Samplers: EC8");
	cs = tra[compleavedep[c][0]].cf; evm.tr = compiftra[c][cs]; evm.t = t; evnew.push_back(evm); c = cs; 

	return simfrom(c,t);
}

double Chain::simfrom(long c, double t)
{	
	long cc, cs, s;
	double tt, r, pr, stnext, prob=0;
	EV evm;

	s = compval[c][settimecl]; stnext = settime[s];
	
	do{
		if(samp_nleave[c] == 0) break;
		
		r = samp_r[c];
	
		tt = t - log(ran())/r; 
	
		while(stnext < tt){
			cs = c+dcsett; cc += dcsett; 
			evm.tr = compiftra[c][cs]; evm.t = stnext; evnew.push_back(evm); c = cs; s++; stnext = settime[s];
		}

		if(tt < tmax){
			prob += log(r) - (tt-t)*r;
			t = tt;

			switch(samp_nleave[c]){
				case 1: cc = samp_leave[c][0]; break;
				case 2:
					pr = samp_prob[c][0];
					if(ran() < pr){ cc = samp_leave[c][0]; prob += log(pr);}
					else{ cc = samp_leave[c][1]; prob += log(1-pr);}
					break;
				default: emsg("Samples: EC3"); break;
			}

			evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
			c = cc;
		}
		else{
			prob += -(tmax-t)*r;
			break;
		}
	}while(1 == 1);
	
	while(stnext < tmax){
		cs = c+dcsett; cc += dcsett; 
		evm.tr = compiftra[c][cs]; evm.t = stnext; evnew.push_back(evm); c = cs; s++; stnext = settime[s];
	}
		
	evm.tr = traend+c; evm.t = tmax; evnew.push_back(evm);

	return prob;
}

double Chain::probsimfrom(vector<EV> &vec, long e)
{
	long c, cc, emax;
	double t, tt, r, pr, prob=0;
	
	c = tra[vec[e].tr].cf;
	
	t = vec[e].t;
	emax = long(vec.size())-1;
	do{
		if(samp_nleave[c] == 0) break;
		
		r = samp_r[c];
		e++; while(e < emax && tra[vec[e].tr].cl != 0) e++;
		tt = vec[e].t;

		if(tt < tmax){
			prob += log(r) - (tt-t)*r;
			t = tt;
	
			cc = tra[vec[e].tr].cf;
			switch(samp_nleave[c]){
				case 1: break;
				case 2:
					pr = samp_prob[c][0];
					if(compval[cc][0] == compval[samp_leave[c][0]][0]) prob += log(pr);
					else{
						prob += log(1-pr);
						
						if(compval[cc][0] !=  compval[samp_leave[c][1]][0])	emsg("Samplers: EC4");
					}
					break;
				default: emsg("Samples: EC3"); break;
			}	
			c = cc;
		}
		else{
			prob += -(tmax-t)*r;
			break;
		}
	}while(1 == 1);
	
	return prob;
}

void Chain::addsamp(long c, long leave, double r)
{
	samp_nleave[c] = 1; 
	samp_leave[c][0] = c+leave-compval[c][0];
	samp_r[c] = r;
	samp_prob[c][0] = 1;
}

void Chain::addsamp2(long c, long leave1, long leave2, double r, double prob)
{
	samp_nleave[c] = 2; 
	samp_leave[c][0] = c+leave1-compval[c][0];
	samp_leave[c][1] = c+leave2-compval[c][0];
	samp_r[c] = r;
	samp_prob[c][0] = prob;
	samp_prob[c][1] = 1-prob;
}

void Chain::addsamprev(long c, long enter, double r)
{ 
	samprev_enter[c] = c+enter-compval[c][0];
	samprev_r[c] = r;
}

double Chain::samplerinit()  // Generates the samplers used when adding individuals or during initialisation
{
	const long SS = 0, EE = 1, AA = 2, II = 3, HH = 4, RR = 5, DD = 6;
	long c, cc, cinit, s;
	double rEA, rAR, rAI, rIR, rIH, rID, rHR, rHD, fac;
	
	rEA = calculatenotdep(tra[compiftra[EE][AA]].eq,param);
	rAR = calculatenotdep(tra[compiftra[AA][RR]].eq,param);
	rAI = calculatenotdep(tra[compiftra[AA][II]].eq,param);
	rIR = calculatenotdep(tra[compiftra[II][RR]].eq,param);
	rIH = calculatenotdep(tra[compiftra[II][HH]].eq,param);
	rID = calculatenotdep(tra[compiftra[II][DD]].eq,param);
	rHR = calculatenotdep(tra[compiftra[HH][RR]].eq,param);
	rHD = calculatenotdep(tra[compiftra[HH][DD]].eq,param);

	samp_r.resize(ncomp); samp_nleave.resize(ncomp); samp_leave.resize(ncomp); samp_prob.resize(ncomp);
	for(c = 0; c < ncomp; c++){ samp_leave[c].resize(2); samp_prob[c].resize(2);}	
	samprev_r.resize(ncomp); samprev_enter.resize(ncomp);
	
	fac = (rIR+rID)/(rIR+rID+rIH);
	for(c = 0; c < ncomp; c++){
		switch(compval[c][0]){
		case SS: case RR: case DD: samp_nleave[c] = 0; samp_r[c] = 0; break;
		case EE: addsamp(c,AA,rEA); break;
		case AA: addsamp2(c,II,RR,rAI+rAR,fac*rAI/(rAI+rAR)); break;
		case II: addsamp2(c,RR,DD,rIR+rID+rIH,rIR/(rIR+rID)); break;
		case HH: addsamp2(c,RR,DD,rHR+rHD,rHR/(rHR+rHD)); break;
		}
		
		switch(compval[c][0]){
		case SS: case RR: case DD: case EE: samprev_enter[c] = -1; samprev_r[c] = 0; break;
		case AA: addsamprev(c,EE,rEA); break;
		case II: addsamprev(c,AA,rAR+rAI); break;
		case HH: addsamprev(c,II,rIH+rIR+rID); break;
		}
	}
	
	fracnotH = 1- (rAI/(rAI+rAR))*(rIH/(rIH+rIR+rID));
	
	
	evnothing.resize(c);
	for(cinit = 0; cinit < ncomp; cinit++){
		c = cinit;
		EV ev; ev.tr = trabeg+c; ev.t = 0; evnothing[cinit].push_back(ev);
		for(s = 0; s < nsettime; s++){
			cc = c+classmult[settimecl];
			EV evm; evm.tr = compiftra[c][cc]; evm.t = settime[s]; evnothing[cinit].push_back(evm);
			c = cc;
		}	
		EV eve; eve.tr = traend+c; eve.t = tmax; evnothing[cinit].push_back(eve);
	}
	
	dcsett = classmult[settimecl];
}
