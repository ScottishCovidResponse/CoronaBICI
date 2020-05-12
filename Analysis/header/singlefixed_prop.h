// This file gives proposals for individuals with a single fixed event specified

void Chain::singlefixed_prop()
{
	long cinit, nob, k, i;
	double al, probif, probfi, Ltoti, Ltotf;
	
	for(cinit = 0; cinit < ncompswa; cinit++){
		nob = obsinflist[cinit].size();
		if(nob > 0){	
			for(k = 0; k < nob; k++){
				if(ran() < 0.01){
					i = obsinflist[cinit][k];
					
					probfi = probsimfromfixed(indev[i],i);
					probif = simfromfixed(i);

					Ltoti = Lwnm();
					indcha(i);
					Ltotf = Lwnm();
					
					al = exp(Ltotf - Ltoti + probfi - probif);
					
					ntr_singfix[cinit]++;
					if(ran() < al){
						nac_singfix[cinit]++;
					}
					else{
						indrev(i);
					}
				}
			}
		}
	}
}

void Chain::singlefixed_prop2()
{
	long cinit, nob, k, kmax, i, d, ddti, e, emax, ddi, ddf;
	double al, almid, pif, pfi, probif, probfi, Ltoti, Ltotf, t;
	vector <long> ist, change, ddist, ddfst;
	vector < vector <EV> > evrevst;	

	for(cinit = 0; cinit < ncompswa; cinit++){
		nob = obsinflist[cinit].size();
		if(nob > 0){
			evrevst.clear(); ist.clear(); change.clear(); ddist.clear(); ddfst.clear();
			probif = 0; probfi = 0;
			
			multich = 1;
			Ltoti = Lwnm();	
		
			kmax = long(sificha[cinit]);
			for(k = 0; k < kmax; k++){
				i = obsinflist[cinit][long(nob*ran())];		
				e = 0; emax = nindev[i]-1; while(e < emax && tra[indev[i][e].tr].cl != 0) e++;
				if(e == emax) emsg("Singledixed: EC1");
					
				t = indev[i][e].t;	
				d = transdepref[tra[indev[i][e].tr].eq];
				ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
				ddi = d*nDD + ddti; if(DDcalc[ddi] == 0) emsg("addreinf: EC88");

				pfi = probsimfromfixed(indev[i],i);
				pif = simfromfixed(i);
				
				e = 0; emax = long(evnew.size())-1; while(e < emax && tra[evnew[e].tr].cl != 0) e++;
				if(e == emax) emsg("Singledixed: EC2");
					
				t = evnew[e].t;	
				d = transdepref[tra[evnew[e].tr].eq];
				ddti = DDconv[long(t*DDfac)]; if(DDt[ddti+1] < t) ddti++;
				ddf = d*nDD + ddti; if(DDcalc[ddf] == 0) emsg("addreinf: EC88");

				ist.push_back(i); ddist.push_back(ddi); ddfst.push_back(ddf);
			
				almid = depeq_disc_evval[ddf]/depeq_disc_evval[ddi]; if(almid > 1) almid = 1;		
				if(ran() < almid){
					probif += log(almid);
					probif += pif; probfi += pfi;						
					change.push_back(1);
					
					indcha(i);
					evrevst.push_back(evrev);
				}
				else{
					probif += log(1-almid);
					change.push_back(0);	
					evrevst.push_back(vector <EV> ());				
				}
			}
			multisec();
			Ltotf = Lwnm();
			
			for(k = 0; k < ist.size(); k++){
				if(change[k] == 1){
					almid = depeq_disc_evval[ddist[k]]/depeq_disc_evval[ddfst[k]]; if(almid > 1) almid = 1;
					probfi += log(almid);
				}
				else{
					almid = depeq_disc_evval[ddfst[k]]/depeq_disc_evval[ddist[k]]; if(almid > 1) almid = 1;
					probfi += log(1-almid);
				}
			}
			
			al = exp(Ltotf - Ltoti + probfi - probif);	
	
			ntr_singfix[cinit]++;
			if(ran() < al){
				nac_singfix[cinit]++;
				if(samp < nsamp){
					sificha[cinit] *= 1.025; 
					if(sificha[cinit] > 20) sificha[cinit] = 20; 
					if(sificha[cinit] > nob/4) sificha[cinit] = nob/4;
				}					
			}
			else{
				if(samp < nsamp) sificha[cinit] *= 0.95;
				
				multich = -1; multisec();
				for(k = ist.size()-1; k >= 0; k--){
					if(change[k] == 1){
						evrev = evrevst[k];
						indrev(ist[k]);
					}
				}
			}

			multisecend();
		}
	}
}
