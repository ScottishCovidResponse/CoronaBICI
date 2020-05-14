// This file considers proposals for simultaneously adding or removing infected individuals (Corona)

void Chain::addreminf()
{
	double upfac = 1.05, downfac = 0.975;
	long checkon = 0;
	long l, i, c, cc, tb, m, si, multac, a, r, s;
	double probif, probfi, Ltoti, Ltotf, al, dd, z;
	double t, tt;
	vector <long> ist;
	vector < vector <EV> > evrevst;	

	for(int cinit = 0; cinit < ncompswa; cinit++){
		if(notinflist[cinit].size() > 0){	
			evrevst.clear(); ist.clear();
			do{ multac = normal(multacf[cinit],multacf[cinit]/10.0);}while(multac <= 0);
			//multac  = 1;
			
			multich = 1;
	
			probif = 0; probfi = 0;
			if(ran() < 0.5){          // Adds a block of infected individuals
				Ltoti = Lwnm();

				getprob(cinit);
				for(m = 0; m < multac; m++){ 	
					tb = tbin-1; z = ran()*probsumtot; while(tb >= 0 && z > probsum[tb]) tb--;
					if(tb < 0) emsg("addreminf: EC2");
					probif += log(prob[tb]/probsumtot);
			
					si = notinflist[cinit].size(); if(si == 0) break;
					l = long(ran()*si); probif += log(1.0/si);
					i = notinflist[cinit][l];
					c = indinit[i];			
				
					t = (tb+ran())*tmax/tbin; probif += log(tbin/tmax);
					
					probif += simfromSE(i,t);
					
					notinflist[cinit][l] = notinflist[cinit][long(notinflist[cinit].size())-1];
					indinflist[notinflist[cinit][l]] = l;
					notinflist[cinit].pop_back();
					
					indinflist[i] = inflist[cinit].size();
					inflist[cinit].push_back(i);		
					
					probfi += log(1.0/(inflist[cinit].size()));	
				
					indcha(i); evrevst.push_back(evrev);
			
					ist.push_back(i);
				}
				multisec();
				Ltotf = Lwnm();
			
				if(m < multac) al = 0;
				else al = exp(Ltotf - Ltoti + probfi - probif);
				
				ntr_multinfadd++;
				if(ran() < al){
					nac_multinfadd++;
					if(samp < burnin) multacf[cinit] *= upfac;
					
					if(checkon == 1 && ran() < 0.1){  // CHECK
						multisecend(); multich = 1;
						
						double probifst = probif, probfist = probfi;
						probif = 0; probfi = 0;
						evrevst.clear();
						
						for(m = 0; m < multac; m++){ 
							si = inflist[cinit].size();
							i = ist[m]; l = indinflist[i];
						
							probif += log(1.0/si);	
						
							evnew = evnothing[cinit];
							indcha(i); evrevst.push_back(evrev);
							
							inflist[cinit][l] = inflist[cinit][long(inflist[cinit].size())-1];
							indinflist[inflist[cinit][l]] = l;
							inflist[cinit].pop_back();
								
							indinflist[i] = notinflist[cinit].size();
							notinflist[cinit].push_back(i);		
								
							probfi += log(1.0/long(notinflist[cinit].size()));
						}

						multisec();
						Ltotf = Lwnm();

						if(m < multac) emsg("prob");
						else{
							getprob(cinit);
							for(m = 0; m < multac; m++){
								size_t e = 1;
								size_t emax = evrevst[m].size();
								while(e < emax && tra[evrevst[m][e].tr].cl != 0)
									e++;
								if(e == emax) emsg("Addfreminf: EC 8");
								
								t = evrevst[m][e].t;
								tb = long(tbin*t/tmax);
								
								probfi += log(prob[tb]/probsumtot);
								probfi += log(tbin/tmax);
								probfi += probsimfrom(evrevst[m],e);
							}	
						}
						dd = probif - probfist; if(dd*dd > tiny) emsg("pr");
						dd = probfi - probifst; if(dd*dd > tiny)emsg("pr2");
					}						
					//CHECK
				}
				else{
					if(samp < burnin){ multacf[cinit] *= downfac; if(multacf[cinit] < 1) multacf[cinit] = 1;}
					
					multich = -1; multisec();
					for(m = long(ist.size())-1; m >= 0; m--){
						i = ist[m];			
						
						evrev = evrevst[m];
						indrev(i);
						
						l = indinflist[i];
						inflist[cinit][l] = inflist[cinit][long(inflist[cinit].size())-1];
						indinflist[inflist[cinit][l]] = l;

						inflist[cinit].pop_back();
						indinflist[i] = notinflist[cinit].size();
						notinflist[cinit].push_back(i);
					}
				}
			}
			else{              // Removes a block of infected individuals
				Ltoti = Lwnm();

				for(m = 0; m < multac; m++){ 
					si = inflist[cinit].size(); if(si == 0) break;
					l = long(ran()*si);
					i = inflist[cinit][l];
					probif += log(1.0/si);	
					evnew = evnothing[cinit];
					indcha(i); evrevst.push_back(evrev);
					ist.push_back(i);
				
					inflist[cinit][l] = inflist[cinit][long(inflist[cinit].size())-1];
					indinflist[inflist[cinit][l]] = l;
					inflist[cinit].pop_back();
						
					indinflist[i] = notinflist[cinit].size();
					notinflist[cinit].push_back(i);		
						
					probfi += log(1.0/long(notinflist[cinit].size()));
				}

				multisec();
				Ltotf = Lwnm();

				if(m < multac) al = 0;
				else{
					getprob(cinit);
					for(m = 0; m < multac; m++){
						int e = 1; 
						int emax = evrevst[m].size();
						while(e < emax && tra[evrevst[m][e].tr].cl != 0) e++;
						if(e == emax) emsg("Addfreminf: EC 8");
						
						t = evrevst[m][e].t;
						tb = long(tbin*t/tmax);
						
						probfi += log(prob[tb]/probsumtot);
						probfi += log(tbin/tmax);
						probfi += probsimfrom(evrevst[m],e);
					}	
					
					al = exp(Ltotf - Ltoti + probfi - probif);
				}
					
				ntr_multinfrem++;
				if(ran() < al){
					nac_multinfrem++;
				}
				else{
					multich = -1; multisec();
					for(m = long(ist.size())-1; m >= 0; m--){
						i = ist[m];			
						evrev = evrevst[m];
						indrev(i);
					
						l = indinflist[i];
						notinflist[cinit][l] = notinflist[cinit][long(notinflist[cinit].size())-1];
						indinflist[notinflist[cinit][l]] = l;

						notinflist[cinit].pop_back();
						indinflist[i] = inflist[cinit].size();
						inflist[cinit].push_back(i);
					}
				}
			}	
			multisecend();
		}
	}
}

void Chain::getprob(long c)      // Calculates probability distribution for samping infection times
{
	long e, ee, tb, d, ci, cf, a, dst[nsettime+1], s, ddti;
	double sum, t, tt;
	
	timeprop[SAMPPROB] -= clock();

	if(ncompleavedep[c] != 1) emsg("addreminf: EC58");
	long tr = compleavedep[c][0];
	ci = tra[tr].ci; cf = tra[tr].cf;
	
	for(s = 0; s <= nsettime; s++){
		dst[s] = transdepref[tra[compiftra[ci+s*dcsett][cf+s*dcsett]].eq]; 
	}
	
	sum = 0;
	for(tb = tbin-1; tb >= 0; tb--){
		t = (tb+0.5)*tmax/tbin;
		
		d = dst[tbins[tb]];
		
		if(discon == 1){
			ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
			prob[tb] = depeq_disc_evval[d*nDD + ddti];
			if(DDcalc[d*nDD + ddti] == 0) emsg("addreinf: EC88");
		}
		else{
			e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];      // finds where we are on the timeline
			do{
				ee = depeq_evnext[e]; tt = depeq_evt[ee]; if(tt > t) break;
				e = ee;
			}while(1 == 1);
			
			prob[tb] = depeq_evval[e];
		}
		sum += prob[tb];
		probsum[tb] = sum; 
	}
	probsumtot = sum;
	
	timeprop[SAMPPROB] += clock();
}

void Chain::addreminfinit()     // Initialises the procedure
{
	long tb, s, c;
	double t;

	for(c = 0; c < nclassval[0]; c++){ if(classval[0][c] == "H") break;}
	if(c == nclassval[0]) emsg("addreminf: EC56");
	HH = c;

	s = 0;
	tbins.resize(tbin);
	for(tb = 0; tb < tbin; tb++){
		t = (tb+0.5)*tmax/tbin;
		while(s < nsettime && settime[s] < t) s++; 
		tbins[tb] = s;
	}
	
	prob.resize(tbin); probsum.resize(tbin);
	
	multacf.resize(ncomp);
	for(c = 0; c < ncomp; c++) multacf[c] = 10;
}

void Chain::addstartuoinf()                                  // Adds infected unobserved individuals in the initial state
{
	long div, ndiv = 20, c, cc, tr, e, i, num, numadd, l;
	double probH, fac, t, tt;
	vector <vector <long> > bin;
	
	bin.resize(ncomp);
	for(c = 0; c < ncomp; c++){
		bin[c].resize(ndiv);
		for(div = 0; div < ndiv; div++) bin[c][div] = 0;
	}
			
  for(i = 0; i < nindtot; i++){
		for(e = 1; e < nindev[i]-1; e++){
			tr = indev[i][e].tr; 
			if(compval[tra[tr].ci][0] != compval[tra[tr].cf][0]){
				div = long(indev[i][e].t*ndiv/tmax); if(div < 0 || div >= ndiv) emsg("addreminf:EC77");
				bin[indinit[i]][div]++;
			}
		}	
	}
}
