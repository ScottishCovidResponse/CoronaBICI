// This file deals with plotting derived quantities

void Chain::deriveplotinit()         // Initialises plots graphs for derived quantities 
{
  long d;

  nderivetemporalplot = 0;
  derivetemporalplot.resize(nderivetemporal);
  for(d = 0; d < nderivetemporal; d++){
    derivetemporalplot[d].resize(DERIVEX);
  }
  bout << "g|" << DERIVEX << "|" << nderivetemporal << "|"; 
  for(d = 0; d < nderivetemporal; d++) bout << derivename[derivetemporal[d]] << "|";
  bout << "\n";
}

void Chain::derivepl()                // Plots graphs for derived quantities (read by BICI GUI)
{
  long dt, i, j, n, v;
  double t, av, CImin, CImax, f;
  vector <double> list;

  bout << "h|";
  for(dt = 0; dt < nderivetemporal; dt++){
    for(i = 0; i < DERIVEX; i++){
      t = onefac*(double(i)/(DERIVEX-1))*tmax2;

      list.clear();
      for(j = long(nderivetemporalplot*0.2); j < nderivetemporalplot; j++){
				list.push_back(derivetemporalplot[dt][i][j]);
			}
      sort(list.begin(),list.end());
      n = list.size(); av = 0; for(j = 0; j < n; j++) av += list[j]; av /= n;
      if (list.size() <= 1) {
        bout << tminactual + t << "|" << av << "|" << av << "|" << av << "|";
        continue;
      }
			
      v = long((n-1)*0.025); f = (n-1)*0.025 - v;
      if ((v + 1) >= list.size()) {
        std::cout << "Resize too big v" << std::endl;
        v = list.size() - 2;
      }
      CImin = list[v]*(1-f) +  list[v+1]*f;

      v = long((n-1)*0.975);
      if ((v + 1) >= list.size()) {
        std::cout << "Resize too big v" << std::endl;
        v = list.size() - 2;
      }
      f = (n-1)*0.975 - v;
      CImax = list[v]*(1-f) +  list[v+1]*f;

      bout << tminactual + t << "|" << CImin << "|" << av << "|" << CImax << "|";
    }
  }
  bout << "\n";
}

void Chain::deriveplotcalc()                       // Calculates data for plots of derived values
{
  long dt, i, j, k, e, tr, p, c;
  double t, popnum[npopnum], sum;
  vector <EV> ev;
  vector <long> nc;

  for(i = 0; i < nindtot_sim; i++){
    for(e = 0; e < nindev_sim[i]; e++) ev.push_back(indev_sim[i][e]);
  }
  sort(ev.begin(),ev.end(),compareev);

	short ccc; for(ccc = 0; ccc < ncomp; ccc++) if(compname[ccc] == "S_S08000030_T<50") break;
	
  nc.resize(ncomp+1); for(c = 0; c < ncomp; c++) nc[c] = 0; nc[NOTALIVE] = nindtot_sim;
  j = 0;
  for(i = 0; i < DERIVEX; i++){
    t = onefac*(double(i)/(DERIVEX-1))*tmax2;
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
		
		for(dt = 0; dt < nderivetemporal; dt++){
			derivetemporalplot[dt][i].push_back(calculate(derive[derivetemporal[dt]],popnum,param));
    }
  }
  nderivetemporalplot++;
}

void deriveinit()                                // Initialises update of derived quantities
{
  long d, eq, j, c, k, m, p;

  dermcomp.resize(nderm);
  for(j = 0; j < nderm; j++){
    dermcomp[j].resize(ncomp); for(c = 0; c < ncomp; c++) dermcomp[j][c] = 0;
    d = derm[j]; eq = derive[d];
    for(k = 0; k < neq_popnum[eq]; k++){
      p = eq_popnum[eq][k];
      for(m = 0; m < popnumterm[p].size(); m++){
        c = popnumterm[p][m]; if(c < 0 || c >= ncomp) emsg("Derive: EC1");
        dermcomp[j][c] = 1;
      }
    }
  }

  nderive_dpop.resize(nderive); derive_dpop.resize(nderive); derive_dpopweight.resize(nderive);
  for(d = 0; d < nderive; d++){
    nderive_dpop[d].resize(ncomp); derive_dpop[d].resize(ncomp); derive_dpopweight[d].resize(ncomp);
    eq = derive[d];
    for(k = 0; k < neq_popnum[eq]; k++){
      p = eq_popnum[eq][k];
      for(m = 0; m < popnumterm[p].size(); m++){
        c = popnumterm[p][m]; if(c < 0 || c >= ncomp) emsg("Derive: EC2");
        derive_dpop[d][c].push_back(p);
        derive_dpopweight[d][c].push_back(popnumtermweight[p][m]);
      }
    }
    for(c = 0; c < ncomp; c++) nderive_dpop[d][c] = derive_dpop[d][c].size();
  }
}
