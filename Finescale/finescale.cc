// Compile with g++ fine.cc -o fine -O3
// Run with ./fine 1 10000 0

// This shows how the number of exposed "E" (infected but not infectious) infected "I" and recovered "R"
// individuals change as a function of time 

// The first number give the approach used (0=old approach  1=new approach);
// The second number gives the number of regions
// The last number changes the random seed


using namespace std;

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

short level;

const short EItrans = 0, IRtrans = 1;

const double tiny = 0.00000001;
const long popsize = 5000000;
long nhouse = 1000;

const short checkon = 0;

const double finegridsize = 0.02;  // The rate over which the best definition is done
//const double finegridsize = 0.05;  // The rate over which the best definition is done
const double d0 = 0.05;
	
const double tmax = 100;

struct HOUSE {
	double x, y;
	vector <long> ind;
};

vector <HOUSE> house;

struct MAPPOINT {
	double val;
	long colref;
};

vector <vector <vector <MAPPOINT> > > map;
vector < vector <vector <long> > > mapadd; 

bool compX(long lhs, long rhs) { return house[lhs].x < house[rhs].x; }
bool compY(long lhs, long rhs) { return house[lhs].y < house[rhs].y; }
	
struct COL {  // collection of houses
	vector <long> houseref;
	long parent;
	vector <long> child;
	double add;
	vector <long> fine;
	double R;
	long pop;
	double x, y;
	short done;
};
 
struct LEVEL {
	vector <COL> collist;
	vector <long> donelist;
	vector <double> add;
};

long Cfine;

vector <LEVEL> lev;

vector <long> indc;
vector <long> indstat;

long tdnext, tdfnext;

vector <double> ffine;
vector <vector <long> > popfine;

struct FEV {                                    // Future events
  long type; long ind; double t; short done; 
};

const long nsettime = 10;

double rEI = 1.0/5, rIR = 1.0/5;

long NE=0, NI=0, NR=0; 

double settime[nsettime] = {10,20,30,40,50,60,70,80,90,100}; 
double R0[nsettime] = {3.6,3.6,3.6,3.6,2.1,1.8,1.5,1,0.5,0.5};
double beta[nsettime];

const short tdiv = 10000;
const long fediv = tdiv;
vector < vector <FEV> > fev;

double ran(){ return double(0.999999999*rand())/RAND_MAX;}
void emsg(string msg){ cout << msg << "\n"; exit (EXIT_FAILURE);}
void addinfc(long c, double t);
void addinfection(long i, double t);
void addfev(double t, long type, long i);
long nextinfection();
void init();
void tauleaping();
void gillespie();
void dofe();

int main(int argc, char** argv)
{
	double t;
	long loop, time = 0, simnum, type;

	if(argc != 4) emsg("wrong num\n");
	
	type = atoi(argv[1]);
	nhouse = atoi(argv[2]);
	simnum = atoi(argv[3]); srand(simnum);// for(loop = 0; loop < simnum*1000; loop++) ran();
	
	cout << "Initialising....\n";
	
	init();	
	
	cout << "Running....\n";
	
	tdnext = tdiv; for(loop = 0; loop < 3; loop++) addinfc(long(ran()*Cfine),0);
	
	time = -clock();

	if(type == 0) tauleaping();
	else gillespie();
	
	time += clock();
	
	cout << double(time)/(CLOCKS_PER_SEC) << " Total time\n";
}

void gillespie()
{
	long td, j, s, c, NIfine[Cfine];
	double t, tpl, R, tnextinf, tnextfe, tnextst;
	
	t = 0; tpl = 1;

	s = 0; tnextst = settime[s];
	do{
		if(tdnext < tdiv) tnextfe = fev[tdnext][tdfnext].t; else tnextfe = tmax;
	
		R = lev[0].collist[0].R;
		if(R < tiny) tnextinf = tmax; 
		else tnextinf = t - log(ran())/(beta[s]*R);
		
		while(t > tpl){ 
			cout  << "Time: " << tpl << "  E: " << NE << "  I: "<< NI << " R: " << NR << "\n"; tpl++;
		}
	
		if(tnextst < tnextfe && tnextst < tnextinf){
			t = tnextst;
			s++; tnextst = settime[s];
		}
		else{
			if(tnextinf < tnextst && tnextinf < tnextfe){
				c = nextinfection();
				addinfc(c,tnextinf);	
				t = tnextinf;
			}
			else{             // transition
				if(tnextfe == tmax) break;
				
				dofe();
				t = tnextfe;
			}
		}
	}while(t < tmax);
	
	for(td = 0; td < tdiv; td++){
		for(j = 0; j < fev[td].size(); j++){
			if(fev[td][j].done == 0) emsg("tt");
		}
	}
}

void dofe()
{
	long i, c, cc, ccc, j, jmax, k, kmax, l, ll, pop;
	double fac, val, num, dd;
	
	i = fev[tdnext][tdfnext].ind;
	fev[tdnext][tdfnext].done = 1;
	c = indc[i];
				
	fac = 0; 			
	switch(fev[tdnext][tdfnext].type){
		case EItrans:
			NE--; NI++; fac = 1;
			break;
			
		case IRtrans:
			NI--; NR++; fac = -1;
			break;
	}

	tdfnext++;
	if(tdfnext == fev[tdnext].size()){
		tdnext++; tdfnext = 0; 
		while(tdnext < tdiv && fev[tdnext].size() == 0) tdnext++;
	}
		
	if(fac == 0) return;
	
	for(l = level-1; l >= 0; l--){
		kmax = map[c][l].size();
		for(k = 0; k < kmax; k++){
			cc = map[c][l][k].colref;
			val = fac*map[c][l][k].val*lev[l].collist[cc].pop;
			lev[l].add[cc] = val;
			lev[l].collist[cc].R += val;
		}
		
		kmax = mapadd[c][l].size();
		for(k = 0; k < kmax; k++){
			cc = mapadd[c][l][k];
			jmax = lev[l].collist[cc].child.size();
			val = 0; for(j = 0; j < jmax; j++) val += lev[l+1].add[lev[l].collist[cc].child[j]];
			lev[l].add[cc] = val;
			lev[l].collist[cc].R += val;
		}

		if(l < level-1){
			kmax = map[c][l].size();
			for(k = 0; k < kmax; k++){
				cc = map[c][l][k].colref;
				val = fac*map[c][l][k].val;
				lev[l].collist[cc].add += val;
			}
		}
	}
	
	if(checkon == 1){
		for(l = level-1; l >= 0; l--){
			kmax = map[c][l].size();
			for(k = 0; k < kmax; k++){
				cc = map[c][l][k].colref;
				val = fac*map[c][l][k].val;
			
				pop = lev[l].collist[cc].pop;
				num = val*pop;
				jmax = lev[l].collist[cc].fine.size();
				for(j = 0; j < jmax; j++){
					ccc = lev[l].collist[cc].fine[j];
					ffine[lev[l].collist[cc].fine[j]] += val;
				}
			}
		}
		
		for(l = level-1; l >= 0; l--){
			kmax = map[c][l].size(); for(k = 0; k < kmax; k++) lev[l].add[map[c][l][k].colref] = 0;
			kmax = mapadd[c][l].size(); for(k = 0; k < kmax; k++) lev[l].add[mapadd[c][l][k]] = 0;
		}
		
		double sum=0;
		for(cc = 0; cc < Cfine; cc++) sum += ffine[cc]*popfine[cc].size();
		dd = lev[0].collist[0].R - sum;
		if(dd < -tiny || dd > tiny){ cout << lev[0].collist[0].R << " " << dd << " dd\n";  emsg("pro2");}
	}
}

long nextinfection()
{
	long l, c, cc, j, jmax;
	double z, sum, sumst[4], val, dd;
	
	l = 0; c = 0;
	do{
		val = lev[l].collist[c].add; lev[l].collist[c].add = 0;
	
		jmax = lev[l].collist[c].child.size();
		sum = 0;
		for(j = 0; j < jmax; j++){
			cc = lev[l].collist[c].child[j];
			
			if(val != 0){
				lev[l+1].collist[cc].R += val*lev[l+1].collist[cc].pop;
				lev[l+1].collist[cc].add += val;
			}	
			sum += lev[l+1].collist[cc].R;
		
			sumst[j] = sum;
		}
		
		z = ran()*sum; j = 0; while(j < jmax && z > sumst[j]) j++;
		if(j == jmax) emsg("j");
		
		c = lev[l].collist[c].child[j];
		l++;
	}while(l < level-1);
	
	if(checkon == 1){
		dd = ffine[c]*popfine[c].size() - lev[l].collist[c].R;
		if(dd < -tiny || dd > tiny) emsg("pro1"); 
	}
	
	return c;
}

long poisson(double lam)
{
  double L, p;
  long k;

  L = exp(-lam); k = 0; p = 1;
  do{
    k++;
    p *= ran();
  }while(p > L);
  return k-1;
}

void tauleaping()
{
	long s, i, j, jmax, c, cc, f, l, h, k, kmax, n, td, pop, NIfine[Cfine];
	double t, dt, fac, val;
	
	for(c = 0; c < Cfine; c++) NIfine[c] = 0;
		
	s = 0;
	ofstream plot("plot.txt");
	dt = tmax/tdiv;
	for(td = 0; td < tdiv; td++){
		t = td*tmax/tdiv; while(settime[s] < t) s++;
		
		if(td%40 == 0) cout  << t << " " << NE << " "<< NI <<" " << NR << "\n";
			
		plot <<  t << " " << NE << " "<< NI <<" " << NR << "\n";
	
		for(c = 0; c < Cfine; c++){
			pop = popfine[c].size();
			n = poisson(beta[s]*ffine[c]*pop*dt);

			for(j = 0; j < n; j++) addinfc(c,(td+ran())*tmax/tdiv);
		}
	
		for(f = 0; f < fev[td].size(); f++){
			i = fev[td][f].ind;
			c = indc[i];
			
			fac = 0; 			
			switch(fev[td][f].type){
				case EItrans:
					NIfine[c]++;
					NE--; NI++; fac = 1;
					break;
					
				case IRtrans:
					NI--; NR++; fac = -1;
					break;
			}

			if(fac != 0){
				for(l = level-1; l >= 0; l--){
					kmax = map[c][l].size();
					for(k = 0; k < kmax; k++){
						cc = map[c][l][k].colref;
						val = fac*map[c][l][k].val;
	 
						jmax = lev[l].collist[cc].fine.size();
						for(j = 0; j < jmax; j++) ffine[lev[l].collist[cc].fine[j]] += val;
					}
				}
			}
		}
	}
	
	for(c = 0; c < Cfine; c++){
		//if(NIfine[c] != 0) cout << c << " " << NIfine[c] << "\n";
	}
}

void init()
{
	long h, l, i, imax, j, jmax, ii, jj, k, kmax, m, mmax, num, c, cmax, cc, ccc, par, s, pop, L, n, flag;
	double x, y, xx, yy, grsi, dd, sum;
	HOUSE ho;
	COL col;
	MAPPOINT mp;
	vector <long> housex1, housex2;
	vector <double> grsize;
	vector <long> grL;
	vector< vector< vector <long> > > grid;
	
	// Randomly distributes houses
	
	for(h = 0; h < nhouse; h++){
		ho.x = ran();
		ho.y = ran();
		house.push_back(ho);
	}

	for(i = 0; i < popsize; i++){
		h = long(ran()*nhouse);
		house[h].ind.push_back(i);
	}
	
	lev.push_back(LEVEL ());
	for(h = 0; h < nhouse; h++) col.houseref.push_back(h);
	col.parent = -1;
	lev[0].collist.push_back(col);

	l = 0;
	do{
		lev.push_back(LEVEL ());
		
		flag = 0;
		
		for(c = 0; c < lev[l].collist.size(); c++){
			sort(lev[l].collist[c].houseref.begin(),lev[l].collist[c].houseref.end(),compX);
			num = lev[l].collist[c].houseref.size();
			
			housex1.clear(); housex2.clear();
			for(j = 0; j < num/2; j++) housex1.push_back(lev[l].collist[c].houseref[j]);
			for(j = num/2; j < num; j++) housex2.push_back(lev[l].collist[c].houseref[j]);
		
			sort(housex1.begin(),housex1.end(),compY);
			num = housex1.size(); if(num > 2) flag = 1;
			
			col.houseref.clear(); col.child.clear(); 
			for(j = 0; j < num/2; j++) col.houseref.push_back(housex1[j]);

			if(col.houseref.size() > 0){
				col.parent = c;
				lev[l].collist[c].child.push_back(lev[l+1].collist.size());	
				lev[l+1].collist.push_back(col);
			}

			col.houseref.clear(); col.child.clear();
			for(j = num/2; j < num; j++) col.houseref.push_back(housex1[j]);
			if(col.houseref.size() > 0){
				col.parent = c;
				lev[l].collist[c].child.push_back(lev[l+1].collist.size());	
				lev[l+1].collist.push_back(col);
			}
			
			sort(housex2.begin(),housex2.end(),compY);
			num = housex2.size(); if(num > 2) flag = 1;
			
			col.houseref.clear(); col.child.clear(); 
			for(j = 0; j < num/2; j++) col.houseref.push_back(housex2[j]);
			if(col.houseref.size() > 0){
				col.parent = c;
				lev[l].collist[c].child.push_back(lev[l+1].collist.size());	
				lev[l+1].collist.push_back(col);
			}

			col.houseref.clear(); col.child.clear();
			for(j = num/2; j < num; j++) col.houseref.push_back(housex2[j]);
			if(col.houseref.size() > 0){
				col.parent = c;
				lev[l].collist[c].child.push_back(lev[l+1].collist.size());	
				lev[l+1].collist.push_back(col);
			}
		}
		l++;
	}while(flag == 1);

	level = l+1; cout << level << " Number of levels\n";
	
	Cfine = lev[l].collist.size();
	
	l = level-1;
	cmax = lev[l].collist.size();
	for(c = 0; c < cmax; c++) lev[l].collist[c].fine.push_back(c);
	
	for(l = level-2; l >= 0; l--){
		cmax = lev[l].collist.size();
		for(c = 0; c < cmax; c++){
			jmax = lev[l].collist[c].child.size();
			for(j = 0; j < jmax; j++){
				cc = lev[l].collist[c].child[j];
				kmax = lev[l+1].collist[cc].fine.size();
				for(k = 0; k < kmax; k++) lev[l].collist[c].fine.push_back(lev[l+1].collist[cc].fine[k]);
			}
		}
	}
	
	/*
	for(l = 0; l < level; l++){
		cout << l << ": ";
			cmax = lev[l].collist.size();
		for(c = 0; c < cmax; c++) cout << lev[l].collist[c].fine.size() << ",";
		cout << "\n";
	}
	*/
	
	l = level-1;
	for(c = 0; c < Cfine; c++){
		pop = 0; x = 0; y = 0;
		
		imax = lev[l].collist[c].houseref.size();
		for(i = 0; i < imax; i++){
			h = lev[l].collist[c].houseref[i];
			pop += house[h].ind.size();
			x += house[h].x;
			y += house[h].y;
		}
		lev[l].collist[c].pop = pop;
		lev[l].collist[c].x = x/imax;
		lev[l].collist[c].y = y/imax;
		lev[l].collist[c].done = 0;
		lev[l].collist[c].R = 0;
		lev[l].collist[c].add = 0;
	}
	
	for(l = level-2; l >= 0; l--){
		for(c = 0; c < lev[l].collist.size(); c++){
			pop = 0; x = 0; y = 0;
			jmax = lev[l].collist[c].child.size();
			for(j = 0; j < jmax; j++){
				cc = lev[l].collist[c].child[j];
				pop += lev[l+1].collist[cc].pop;
				x += lev[l+1].collist[cc].x;
				y += lev[l+1].collist[cc].y;
			}
			lev[l].collist[c].pop = pop;
			lev[l].collist[c].x = x/jmax;
			lev[l].collist[c].y = y/jmax;
			lev[l].collist[c].done = 0;
			lev[l].collist[c].R = 0;
			lev[l].collist[c].add = 0;
		}
	}
	
	grsi = finegridsize;
	grsize.resize(level); grL.resize(level); grid.resize(level);
	for(l = level-1; l >= 0; l--){
		L = long(1.0/(2*grsi)); if(L < 1) L = 1;
		grsize[l] = grsi;
		grL[l] = L;
		cmax = lev[l].collist.size();
		grid[l].resize(L*L);
		for(c = 0; c < cmax; c++){
			i = long(lev[l].collist[c].x*L);
			j = long(lev[l].collist[c].x*L);	
			grid[l][j*L+i].push_back(c);
		}
		grsi *= 2;
	}
	
	map.resize(Cfine); mapadd.resize(Cfine);
	for(c = 0; c < Cfine; c++){
		map[c].resize(level); mapadd[c].resize(level);
		
		if(c%1000 == 0) cout << c << " "  << Cfine << "\n";
		l = level-1;
		
		x = lev[l].collist[c].x;
		y = lev[l].collist[c].y;
	
		for(l = level-1; l >= 0; l--){
			flag = 0;
			
			if(l < level-1){
				kmax = lev[l+1].donelist.size();
				for(k = 0; k < kmax; k++){
					cc = lev[l+1].donelist[k];
					par = lev[l+1].collist[cc].parent; if(par == -1) emsg("prob");
					if(lev[l].collist[par].done == 0){
						lev[l].collist[par].done = 1;
						lev[l].donelist.push_back(par);
										
						mmax = lev[l].collist[par].child.size();
						for(m = 0; m < mmax; m++){
							ccc = lev[l].collist[par].child[m];
							if(lev[l+1].collist[ccc].done == 0){
								xx = lev[l+1].collist[ccc].x;
								yy = lev[l+1].collist[ccc].y;
								dd = ((xx-x)*(xx-x) + (yy-y)*(yy-y))/(d0*d0); if(dd < 0.2) dd = 0.2;
														
								mp.val = 1.0/dd;
								mp.colref = ccc;
								map[c][l+1].push_back(mp);
								flag = 1;
							}
						}
					}
				}
			}
			
			L = grL[l];	grsi = grsize[l];
			i = long(L*x+0.5)-1;
			j = long(L*x+0.5)-1;

			for(ii = i; ii <= i+1; ii++){
				for(jj = j; jj <= j+1; jj++){
					if(ii >= 0 && ii < L && jj >= 0 && jj < L){
						n = jj*L+ii;
						kmax = grid[l][n].size();
						for(k = 0; k < kmax; k++){
							cc = grid[l][n][k];
							if(lev[l].collist[cc].done == 0){
								xx = lev[l].collist[cc].x;
								yy = lev[l].collist[cc].y;
		
								dd = (xx-x)*(xx-x) + (yy-y)*(yy-y);
								if(dd < grsi*grsi){
									dd /= d0*d0; if(dd < 0.2) dd = 0.2;
					
									mp.val = 1.0/dd;
									mp.colref = cc;
									map[c][l].push_back(mp);
									flag = 1;
									
									lev[l].collist[cc].done = 1;
									lev[l].donelist.push_back(cc);
								}
							}
						}
					}
				}
			}
		}
	
		sum = 0;
		for(l = 0; l < level; l++){
			kmax = map[c][l].size(); 
			for(k = 0; k < kmax; k++) sum += map[c][l][k].val*lev[l].collist[map[c][l][k].colref].pop;
		}
		for(l = 0; l < level; l++){
			kmax = map[c][l].size(); for(k = 0; k < kmax; k++) map[c][l][k].val /= sum;
		}
		
		pop = 0;
		for(l = 0; l < level; l++){
			kmax = map[c][l].size(); 
			for(k = 0; k < kmax; k++) pop += lev[l].collist[map[c][l][k].colref].pop;
		}
		if(pop != popsize) emsg("popsize wrong");
		
		for(l = 0; l < level; l++){
			kmax = map[c][l].size(); 
			for(k = 0; k < kmax; k++) lev[l].collist[map[c][l][k].colref].done = 0;
			
			jmax = lev[l].donelist.size(); 
			for(j = 0; j < jmax; j++){
				cc = lev[l].donelist[j];
				if(lev[l].collist[cc].done == 1) mapadd[c][l].push_back(cc);
			}
		}			
		
		for(l = 0; l < level; l++){
			kmax = lev[l].donelist.size(); for(k = 0; k < kmax; k++) lev[l].collist[lev[l].donelist[k]].done = 0;
			lev[l].donelist.clear();
		}
		
		for(l = 0; l < level; l++){
			cmax = lev[l].collist.size();
			for(cc = 0; cc < cmax; cc++) if(lev[l].collist[cc].done != 0) emsg("done prob");
		}
	}
	
	ffine.resize(Cfine); popfine.resize(Cfine); indc.resize(popsize);
	for(c = 0; c < Cfine; c++){
		ffine[c] = 0;
		if(lev[level-1].collist[c].houseref.size() != 1) emsg("P");		
		h = lev[level-1].collist[c].houseref[0];
		for(i = 0; i < house[h].ind.size(); i++){
			popfine[c].push_back(house[h].ind[i]);
			indc[house[h].ind[i]] = c;
		}
	}

	fev.resize(fediv);

	for(s = 0; s < nsettime; s++) beta[s] = R0[s]*rIR;	
	
	for(l = 0; l < level; l++){
		cmax = lev[l].collist.size();
		for(c = 0; c < cmax; c++) lev[l].add.push_back(0); 
	}
}

void addinfc(long c, double t)
{
	long pop, l, i, cc;
	double dR;
	
	pop = popfine[c].size(); if(pop == 0){ cout << lev[level-1].collist[c].R << "R\n"; emsg("no pop");}
	l = long(ran()*pop);
	i = popfine[c][l];
	
	popfine[c][l] = popfine[c][pop-1];
	popfine[c].pop_back();
	
	l = level-1; cc = c;
	dR = -lev[l].collist[c].R/pop;
	do{
		lev[l].collist[cc].R += dR;
		lev[l].collist[cc].pop--;
		cc = lev[l].collist[cc].parent; l--;
	}while(l >= 0);
	
	addinfection(i,t);
}
	
void addinfection(long i, double t)
{
	double tI, tR;

	NE++;
	tI = t - log(ran())/rEI;
	tR = tI - log(ran())/rIR;

	addfev(tI,EItrans,i);
	addfev(tR,IRtrans,i);
}

void addfev(double t, long type, long i)
{
	long d, j;
	
	if(t >= tmax) return;
	
	d = long((t/tmax)*fediv);
	j = 0; while(j < fev[d].size() && t > fev[d][j].t) j++;
	
	if(d  == tdnext){ if(j < tdfnext) tdfnext = j;}
	if(d < tdnext){ tdnext = d; tdfnext = j;}
	
	FEV fe; fe.t = t; fe.type = type; fe.ind = i; fe.done = 0;
	fev[d].insert(fev[d].begin()+j,fe);
}
