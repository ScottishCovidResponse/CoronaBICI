// This file considers proposals for simultaneously adding or removing infected individuals (Corona)

void Chain::addreminf()
{
	double upfac = 1.05, downfac = 0.975;
	long checkon = 0;
	long l, i, c, cc, tb, m, si, multac, a, r, s;
	double probif, probfi, Ltoti, Ltotf, al, dd, z;
	double t, tt;
	vector <long> ist;
	vector <short> cst;
	vector < vector <EV> > evrevst;	

	for(r = 0; r < nregion; r++){
		evrevst.clear(); ist.clear(); cst.clear();

		do{ multac = normal(multacf[r],multacf[r]/10.0);}while(multac <= 0);
		
		multich = 1;
		
		probif = 0; probfi = 0;
		if(ran() < 0.5){          // Adds a block of infected individuals
			Ltoti = Lwnm();

			getprob(r);
			for(m = 0; m < multac; m++){ 
			  a = 0; z = ran()*probsumatot; while(a < nag && z > probsuma[a]) a++;
				if(a == nag) emsg("addreminf: EC1");
				
				tb = tbin-1; z = ran()*probsum[a][0]; while(tb >= 0 && z > probsum[a][tb]) tb--;
				if(tb < 0) emsg("addreminf: EC2");
					
				probif += log(prob[a][tb]/probsumatot);
		
				si = notinflist[r][a].size(); if(si == 0) break;
				l = long(ran()*si); probif += log(1.0/si);
				i = notinflist[r][a][l];
				c = indinit[i];			
			
				t = (tb+ran())*tmax/tbin; probif += log(tbin/tmax);
				evnew.clear();
				EV ev; ev.tr = trabeg+c; ev.t = 0; evnew.push_back(ev);
				cc = c+EE-SS;
				EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
				c = cc;
			
				tt = t - log(ran())/rEA;
				if(tt >= tmax) probif += -rEA*(tmax-t);
				else{
					probif += log(rEA) - rEA*(tt-t); 
					t = tt;
					cc = c + AA - EE;
					EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
					c = cc;
					
					tt = t - log(ran())/(rAR+rAIcor);
					if(tt >= tmax) probif += -(rAR+rAIcor)*(tmax-t);
					else{
						if(ran()*(rAR+rAIcor) < rAR){
							probif += log(rAR) - (rAR+rAIcor)*(tt-t);
							t = tt;
							cc = c + RR-AA;
							EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
							c = cc;
						}
						else{
							probif += log(rAIcor) - (rAR+rAIcor)*(tt-t);
							t = tt;
							cc = c + II-AA;
							EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
							c = cc;
							
							tt = t -log(ran())/(rIR[a]+rID[a]);
							if(tt >= tmax) probif += -(rIR[a]+rID[a])*(tmax-t);
							else{
								if(ran()*(rIR[a]+rID[a]) < rIR[a]){
									probif += log(rIR[a]) - (rIR[a]+rID[a])*(tt-t);
									t = tt;
									cc = c + RR-II;
									EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
									c = cc;
								}
								else{
									probif += log(rID[a]) - (rIR[a]+rID[a])*(tt-t);
									t = tt;
									cc = c + DD-II;
									EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
									c = cc;
								}
							}
						}
					}
				}
				EV eve; eve.tr = traend+c; eve.t = tmax; evnew.push_back(eve);	
				if(nsettime > 0) addsettime();
			
				notinflist[r][a][l] = notinflist[r][a][long(notinflist[r][a].size())-1];
				indinflist[notinflist[r][a][l]] = l;
				notinflist[r][a].pop_back();
				
				indinflist[i] = inflist[r][a].size();
				inflist[r][a].push_back(i);		
				
				probfi += log(1.0/(inflist[r][a].size()));	
			
				indcha(i); evrevst.push_back(evrev);
		
				ist.push_back(i);
			}
			multisec();
			Ltotf = Lwnm();
		
			if(m < multac) al = 0;
			else{
				getprob(r);
				for(m = 0; m < multac; m++){ 
					i = ist[m];
					a = (indinit[i]/(nclassval[0]*nclassval[1]))%nag;
					probfi += log(proba[a]/probsumatot);
				}
				al = exp(Ltotf - Ltoti + probfi - probif);
			}
			
			ntr_multinfadd++;
			if(ran() < al){
				nac_multinfadd++;
				if(samp < burnin) multacf[r] *= upfac;
				
				if(checkon == 1){  
					// CHECK Used for checking proposals are correct
					double probifst = probif, probfist = probfi;
					probif = 0; probfi = 0;
					vector < vector <EV> > evrevst2;
				
					multich = -1;
					getprob(r);
					for(m = 0; m < multac; m++){ 
						i = ist[m];
						a = (indinit[i]/(nclassval[0]*nclassval[1]))%nag;
						l = indinflist[i]; 
					
						probif += log(proba[a]/probsumatot);
						probif += log(1.0/(inflist[r][a].size()));	
					
						c = indinit[i]; cst.push_back(c);
					
						evnew.clear();
						EV ev; ev.tr = trabeg+c; ev.t = 0; evnew.push_back(ev);
						EV eve; eve.tr = traend+c; eve.t = tmax; evnew.push_back(eve);
						
						indcha(i); evrevst2.push_back(evrev);
						ist.push_back(i);
					
						inflist[r][a][l] = inflist[r][a][long(inflist[r][a].size())-1];
						indinflist[inflist[r][a][l]] = l;
						inflist[r][a].pop_back();
							
						indinflist[i] = notinflist[r][a].size();
						notinflist[r][a].push_back(i);		
							
						probfi += log(1.0/long(notinflist[r][a].size()));
					}
					multisec();

					getprob(r);
					for(m = 0; m < multac; m++){
						evnew = evrevst2[m];
						if(nsettime > 0) remsettime();
						
						c = cst[m];
						a = (c/(nclassval[0]*nclassval[1]))%nag;
						
						t = evnew[1].t;
						tb = long(tbin*t/tmax);
					
						probfi += log(prob[a][tb]/probsumatot);
						probfi += log(tbin/tmax);
					
						tt = evnew[2].t;
						probfi += -rEA*(tt-t);
						if(tt < tmax){
							probfi += log(rEA);
							t = tt;
							tt = evnew[3].t;
							probfi += -(rAR+rAIcor)*(tt-t);
							if(tt < tmax){
								t = tt;
								if(evnew[3].tr == compiftra[c+AA][c+RR]) probfi += log(rAR);
								else{
									if(evnew[3].tr != compiftra[c+AA][c+II]) emsg("addreminf: EC3");
									probfi += log(rAIcor);
									tt = evnew[4].t;
									probfi +=  -(rIR[a]+rID[a])*(tt-t);
									if(tt < tmax){
										if(evnew[4].tr == compiftra[c+II][c+RR]) probfi += log(rIR[a]);
										else{
											if(evnew[4].tr != compiftra[c+II][c+DD]) emsg("addreminf: EC3");
											probfi += log(rID[a]);
										}
									}
								}
							}
						}
					}
					dd = probif - probfist; if(dd*dd > 0.0000001) emsg("addreminf: EC4");
					dd = probfi - probifst; if(dd*dd > 0.0000001) emsg("addreminf: EC5");
					
					//cout << probif << " "<< probfi << " " << probifst <<" " << probfist << " Check\n";
					// CHECKEND
				}
			}
			else{
				if(samp < burnin){ multacf[r] *= downfac; if(multacf[r] < 1) multacf[r] = 1;}
				
				multich = -1; multisec();
				for(m = long(ist.size())-1; m >= 0; m--){
					i = ist[m];			
					
					evrev = evrevst[m];
					indrev(i);
					
					l = indinflist[i];
					a = (indinit[i]/(nclassval[0]*nclassval[1]))%nag;
		
					inflist[r][a][l] = inflist[r][a][long(inflist[r][a].size())-1];
					indinflist[inflist[r][a][l]] = l;

					inflist[r][a].pop_back();
					indinflist[i] = notinflist[r][a].size();
					notinflist[r][a].push_back(i);
				}
			}
		}
		else{              // Removes a block of infected individuals
			Ltoti = Lwnm();

			getprob(r);
			for(m = 0; m < multac; m++){ 
				a = 0; z = ran()*probsumatot; while(a < nag && z > probsuma[a]) a++;
				if(a == nag) emsg("addreminf: EC6");
				
				si = inflist[r][a].size(); if(si == 0) break;
				l = long(ran()*si);
				i = inflist[r][a][l];
			
				probif += log(proba[a]/probsumatot);
				probif += log(1.0/si);	
			
				c = indinit[i]; cst.push_back(c);
			
				evnew.clear();
				EV ev; ev.tr = trabeg+c; ev.t = 0; evnew.push_back(ev);
				for(s = 0; s < nsettime; s++){
					cc = c+classmult[settimecl];
					ev.tr =  compiftra[c][cc]; ev.t = settime[s]; evnew.push_back(ev);
					c = cc; 
				}
				EV eve; eve.tr = traend+c; eve.t = tmax; evnew.push_back(eve);
				
				indcha(i); evrevst.push_back(evrev);
				ist.push_back(i);
			
				inflist[r][a][l] = inflist[r][a][long(inflist[r][a].size())-1];
				indinflist[inflist[r][a][l]] = l;
				inflist[r][a].pop_back();
					
				indinflist[i] = notinflist[r][a].size();
				notinflist[r][a].push_back(i);		
					
				probfi += log(1.0/long(notinflist[r][a].size()));
			}

			multisec();
			Ltotf = Lwnm();

			if(m < multac) al = 0;
			else{
				getprob(r);
				for(m = 0; m < multac; m++){
					evnew = evrevst[m];
					if(nsettime > 0) remsettime();
					
					c = cst[m];
					a = (c/(nclassval[0]*nclassval[1]))%nag;
					
					t = evnew[1].t;
					tb = long(tbin*t/tmax);
					
					probfi += log(prob[a][tb]/probsumatot);
					probfi += log(tbin/tmax);
					
					tt = evnew[2].t;
					probfi += -rEA*(tt-t);
					if(tt < tmax){
						probfi += log(rEA);
						t = tt;
						tt = evnew[3].t;
						probfi += -(rAR+rAIcor)*(tt-t);
						if(tt < tmax){
							t = tt;
							if(evnew[3].tr == compiftra[c+AA][c+RR]) probfi += log(rAR);
							else{
								if(evnew[3].tr != compiftra[c+AA][c+II]) emsg("addreminf: EC7");
								probfi += log(rAIcor);
								tt = evnew[4].t;
								probfi +=  -(rIR[a]+rID[a])*(tt-t);
								if(tt < tmax){
									if(evnew[4].tr == compiftra[c+II][c+RR]) probfi += log(rIR[a]);
									else{
										if(evnew[4].tr != compiftra[c+II][c+DD]) emsg("addreminf: EC8");
										probfi += log(rID[a]);
									}
								}
							}
						}
					}
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
					a = (indinit[i]/(nclassval[0]*nclassval[1]))%nag;
					
					notinflist[r][a][l] = notinflist[r][a][long(notinflist[r][a].size())-1];
					indinflist[notinflist[r][a][l]] = l;

					notinflist[r][a].pop_back();
					indinflist[i] = inflist[r][a].size();
					inflist[r][a].push_back(i);
				}
			}
		}	
		multisecend();
	}
}

void Chain::getprob(long r)      // Calculates probability distribution for samping infection times
{
	long e, ee, tb, d, ci, cf, a, dst[nsettime+1], s, ddti;
	double sum, suma, t, tt;
	
	timeprop[SAMPPROB] -= clock();
	suma = 0;
	for(a = 0; a < nag; a++){
		ci = SS*classmult[0] + r*classmult[1] + a*classmult[2];
		cf = EE*classmult[0] + r*classmult[1] + a*classmult[2];
		
		for(s = 0; s <= nsettime; s++){
			dst[s] = transdepref[tra[compiftra[ci+s*classmult[settimecl]][cf+s*classmult[settimecl]]].eq]; 
		}
		
		sum = 0;
		for(tb = tbin-1; tb >= 0; tb--){
			t = (tb+0.5)*tmax/tbin;
			
			d = dst[tbins[tb]];
			
			if(discon == 1){
				ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
				prob[a][tb] = depeq_disc_evval[d*nDD + ddti];
				if(DDcalc[d*nDD + ddti] == 0) emsg("addreinf: EC88");
			}
			else{
				e = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];      // finds where we are on the timeline
				do{
					ee = depeq_evnext[e]; tt = depeq_evt[ee]; if(tt > t) break;
					e = ee;
				}while(1 == 1);
				
				prob[a][tb] = depeq_evval[e];
			}
			sum += prob[a][tb];
			probsum[a][tb] = sum; 
		}
		
		suma += sum;
		proba[a] = sum;
		probsuma[a] = suma;
	}
	probsumatot = suma;
	
	timeprop[SAMPPROB] += clock();
}

void Chain::addreminfinit()     // Initialises the procedure
{
	long tb, r, a, s;
	double t;

	nregion = nclassval[1];
	nag = nclassval[2];
	
	rEA = calculatenotdep(tra[compiftra[EE][AA]].eq,param);
	rAR = calculatenotdep(tra[compiftra[AA][RR]].eq,param);
	rAI = calculatenotdep(tra[compiftra[AA][II]].eq,param);
	for(a = 0; a < nag; a++){
		rIR[a] = calculatenotdep(tra[compiftra[II+a*classmult[2]][RR+a*classmult[2]]].eq,param);
		rIH[a] = calculatenotdep(tra[compiftra[II+a*classmult[2]][HH+a*classmult[2]]].eq,param);
		rID[a] = calculatenotdep(tra[compiftra[II+a*classmult[2]][DD+a*classmult[2]]].eq,param);
		rHR[a] = calculatenotdep(tra[compiftra[HH+a*classmult[2]][RR+a*classmult[2]]].eq,param);
		rHD[a] = calculatenotdep(tra[compiftra[HH+a*classmult[2]][DD+a*classmult[2]]].eq,param);
	}

	rAIcor = rAI*(1-(rIH[0]/(rIH[0]+rIR[0]+rID[0]))); 

	s = 0;
	tbins.resize(tbin);
	for(tb = 0; tb < tbin; tb++){
		t = (tb+0.5)*tmax/tbin;
		while(s < nsettime && settime[s] < t) s++; 
		tbins[tb] = s;
	}
	
	prob.resize(nag); probsum.resize(nag); proba.resize(nag); probsuma.resize(nag);
	for(a = 0; a < nag; a++){
		prob[a].resize(tbin); probsum[a].resize(tbin);
	}
	
	multacf.resize(nregion);
	for(r = 0; r < nregion; r++) multacf[r] = 10;
}

void Chain::addstartuoinf()                                  // Adds infected unobserved individuals in the initial state
{
	long r, div, ndiv = 20, c, cc, tr, e, i, num, numadd, a, l;
	double probH, fac, t, tt;
	vector <vector <long> > bin;
	
	bin.resize(nregion);
	for(r = 0; r < nregion; r++){
		bin[r].resize(ndiv);
		for(div = 0; div < ndiv; div++) bin[r][div] = 0;
	}
			
  for(i = 0; i < nindtot; i++){
		c = indinit[i];
		r = compval[c][1];
		for(e = 1; e < nindev[i]-1; e++){
			tr = indev[i][e].tr; 
			if(compval[tra[tr].ci][0] == SS && compval[tra[tr].cf][0] == EE){
				div = long(indev[i][e].t*ndiv/tmax); if(div < 0 || div >= ndiv) emsg("addreminf:EC77");
				bin[r][div]++;
			}
		}	
	}
	
	a = 0; if(nag != 1) emsg("addreminf: EC66");

	probH = (rAI/(rAR+rAI))*(rIH[0]/(rIH[0]+rID[0]+rIR[0]));
	fac = (1-probH)/probH;
	
	for(r = 0; r < nregion; r++){
		for(div = 0; div < ndiv; div++){
			numadd = long((bin[r][div]*fac) + 0.5);	
			
			for(num = 0; num < numadd; num++){
				t = (div+ran())*tmax/ndiv;
				if(notinflist[r][a].size()==0) emsg("addreminf: EC67");
				l = long(ran()*notinflist[r][a].size());
				i = notinflist[r][a][l];
				
				notinflist[r][a][l] = notinflist[r][a][long(notinflist[r][a].size())-1];
				indinflist[notinflist[r][a][l]] = l;
				notinflist[r][a].pop_back();
				
				indinflist[i] = inflist[r][a].size();
				inflist[r][a].push_back(i);		
				
				c = indinit[i];
				
				evnew.clear();
				EV ev; ev.tr = trabeg+c; ev.t = 0; evnew.push_back(ev);
				cc = c+EE-SS;
				EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
				c = cc;
			
				tt = t - log(ran())/rEA;
				if(tt < tmax){
					t = tt;
					cc = c + AA - EE;
					EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
					c = cc;
					
					tt = t - log(ran())/(rAR+rAIcor);
					if(tt < tmax){
						if(ran()*(rAR+rAIcor) < rAR){
							t = tt;
							cc = c + RR-AA;
							EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
							c = cc;
						}
						else{
							t = tt;
							cc = c + II-AA;
							EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
							c = cc;
							
							tt = t -log(ran())/(rIR[a]+rID[a]);
							if(tt < tmax){
								if(ran()*(rIR[a]+rID[a]) < rIR[a]){
									t = tt;
									cc = c + RR-II;
									EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
									c = cc;
								}
								else{
									t = tt;
									cc = c + DD-II;
									EV evm; evm.tr = compiftra[c][cc]; evm.t = t; evnew.push_back(evm);
									c = cc;
								}
							}
						}
					}
				}
				EV eve; eve.tr = traend+c; eve.t = tmax; evnew.push_back(eve);	
				
				if(nsettime > 0) addsettime();
		
				indcha(i);
			}
		}	
	}
}
