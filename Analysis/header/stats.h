// Combines information from diffirent chaims to derived statistics

void Chain::addparam()
{
	long p;

	for(p = 0; p < nparam; p++){
		paramst[simnum][p].push_back(param[p]);
	}
}

void stats()                                          // Calculates diagnostic statistics
{
	long r, i, p, N, n, d;
	double mu[nrun], vari[nrun], W, B, valav, varr, ESS, GR, CImin, CImax, a, cor, f, muav, sum;
	vector <double> temp;
	
	N = paramst[simnum][0].size();
	
	#ifdef MPI
	double *block;
	block = new double[N*nparam];
	
	if(simnum == 0){
		for(r = 1; r < nrun; r++){
		  MPI_Recv(block,nparam*N,MPI_DOUBLE,r,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(p = 0; p < nparam; p++){ for(i = 0; i < N; i++) paramst[r][p].push_back(block[p*N+i]);} 
		}
	}
	else{
		for(p = 0; p < nparam; p++){ for(i = 0; i < N; i++) block[p*N+i] = paramst[simnum][p][i];} 
	  MPI_Send(block, nparam*N,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	}
	#endif
	
	if(simnum == 0){
		stringstream sdi; sdi << root << "/stats" << ".txt";
		ofstream stats(sdi.str().c_str());

		stats << "Parameter\tMean\t90% Credible Interval\tESS\tGR\n";
		
		for(p = 0; p < nparam; p++){
			valav = 0; varr = 0;
			for(r = 0; r < nrun; r++){ for(i = 0; i < N; i++) valav += paramst[r][p][i]/(N*nrun);}
			for(r = 0; r < nrun; r++){ 
				for(i = 0; i < N; i++) varr += (paramst[r][p][i]-valav)*(paramst[r][p][i]-valav)/(N*nrun-1);
			}
			if(varr < 0.0000000000000001){
				ESS = -1; GR = -1;
			}
			else{
				temp.clear();
				for(r = 0; r < nrun; r++){
					for(i = 0; i < N; i++) temp.push_back((paramst[r][p][i]-valav)/sqrt(varr));
				}
				
				n = temp.size();
				sum = 1;
				for(d = 0; d < n/2; d++){             // calculates the effective sample size
					a = 0; for(i = 0; i < n-d; i++) a += temp[i]*temp[i+d]; 
					cor = a/(n-d);
					if(cor < 0) break;
					sum += 0.5*cor;			
				}
				ESS = n/sum;
				
				if(nrun == 1) GR = -1;
				else{
					muav = 0;
					for(r = 0; r < nrun; r++){
						valav = 0; 
						for(i = 0; i < N; i++) valav += paramst[r][p][i]/N;
						varr = 0; 
						for(i = 0; i < N; i++) varr += (paramst[r][p][i]-valav)*( paramst[r][p][i]-valav)/(N-1);
						mu[r] = valav;
						vari[r] = varr;
						muav += mu[r]/nrun;
					}
					W = 0; for(r = 0; r < nrun; r++) W += vari[r]/nrun;
					B = 0; for(r = 0; r < nrun; r++) B += (mu[r]-muav)*(mu[r]-muav)*N/(nrun-1);
					GR = sqrt(((1-1.0/N)*W+B/N)/W);
				}
			}
			
			temp.clear();
			for(r = 0; r < nrun; r++){
				for(i = 0; i < N; i++) temp.push_back(paramst[r][p][i]);
			}
			sort(temp.begin(),temp.end());
			
			n = nrun*N;
			i = long((n-1)*0.05); f = (n-1)*0.05 - i;
			CImin = temp[i]*(1-f) + temp[i+1]*f;
			
			i = long((n-1)*0.95); f = (n-1)*0.95 - i;
			CImax = temp[i]*(1-f) + temp[i+1]*f;
			
			stats << paramname[p] << "\t" << valav << "\t" << CImin << " - " << CImax << "\t";
			if(ESS == -1) stats << "---\t"; else stats	<< ESS << "\t";
			if(GR == -1) stats << "---\n"; else stats << GR << "\n";
		}
	}
}
