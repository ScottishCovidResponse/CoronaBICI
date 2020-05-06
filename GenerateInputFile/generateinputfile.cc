// This code takes the model defined in BICI and adds data to it based on case and death data
// The output is an XML file which can the be read into BICI for analysis
 
// compile using  g++ generateinputfile.cc -o gen
// To create model run using ./gen Scotland model    
// To create input file run using ./gen Scotland input 
// Note, "Scotland" can also set to "UK", "England" or "Wales"

#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

double tEA = 2.5, tAI = 3, tAR = 3.5, tIR = 16.7;  // Estimates of transition times from literature
double tIH = 3.1, tID = 12.9, tHD = 8, tHR = 13.1;

short nage = 1;                                    // Different age classifications
string age[1]={"allage"};
double N[1]={100};                                 // Percentage in different age categories
double ph[1]={1};                                  // Probability of hostipalisation

double tmax;                                       // The time over which observations are made

void readdata(string cont);                        // Functions
void loadcensus();
void loadpostcodes();
void createbici(string geo);
void generateM();
void generatexml(string geo);
void displaymap();

double ran(){ return double(rand())/RAND_MAX;};

struct REG {                                       // Provides information about a region
	string name;         
	string code;
	double x, y;
	long pop;
	long popscot;
	long numcase;
	vector <long> inf; 
	vector <double> pcx, pcy;	
};

long nregion;
vector <REG> region;
	
vector< vector<double> > mat;                      // Matrix giving intermixing between regions
	
int main(int argc, char *argv[])
{
	long r, totpop;
	
	if(argc != 3){ cout << "Wrong number of arguments!\n"; return 0;}           
	string geo(argv[1]);
	
	readdata(geo);                            // Reads in the case data
	
	loadcensus();                             // Uses census data to get population information

	string st(argv[2]);
	if(st == "model") createbici(geo);   // This creates a BICI model file
	else generatexml(geo);                    // This adds data to a model specification
	
	cout << nregion << " regions:\n";
	totpop = 0;
	for(r = 0; r < nregion; r++){
		cout << region[r].name << "  Code: " << region[r].code << "  # Cases: " <<  region[r].numcase << " " << "   Population: " << region[r].pop << " " << region[r].x << "," << region[r].y << "\n";
		totpop += region[r].pop;
	}
	cout << "Total population: " << totpop << "\n";
}

void readdata(string cont)                               // Reads in the case data
{
	string date, country, code, name, line, nums;
	long num, r, i;
	double t;
	
	nregion = 0;

	ifstream input("Data/covid-19-cases-uk-02.05.txt");
	getline(input,line);
	tmax = 0;
  do{
		getline(input,line);
		if(input.eof()) break;
		
		stringstream ss(line);		
		getline(ss,date,'\t');

		// converts the date to a time
		if(date.substr(4,1) == "1") t = atoi(line.substr(0,2).c_str())-1; 
		else{
			if(date.substr(4,1) == "2") t = 31+atoi(line.substr(0,2).c_str())-1; 
			else{
				if(date.substr(4,1) == "3") t = 31+29+atoi(line.substr(0,2).c_str())-1; 
				else{
					if(line.substr(4,1) == "4") t = 31+29+31+atoi(line.substr(0,2).c_str())-1;
					else{
						if(line.substr(4,1) == "5") t = 31+29+31+30+atoi(line.substr(0,2).c_str())-1;
						else cout << line.substr(4,1) << "prob\n";
					}
				}
			}
		}
					
		getline(ss,country,'\t');
		getline(ss,code,'\t');
		getline(ss,name,'\t');
		getline(ss,nums,'\t'); num = atoi(nums.c_str());
		
		if((country == cont || cont == "UK" || (cont == "Single" && code == "S08000030")) && code != ""){
			for(r = 0; r < nregion; r++) if(region[r].code == code) break;
			if(r == nregion){
				REG reg; reg.name = name; reg.code = code; reg.numcase = 0; reg.pop = 0; reg.popscot = 0;
				region.push_back(reg); nregion++;
			}
		
			if(num > region[r].numcase){
				for(i = region[r].numcase; i < num; i++){
					region[r].inf.push_back(t);
					if(t > tmax) tmax = t;
				}
				region[r].numcase = num;
			}
		}
	}while(1 == 1);
	tmax++;
}

void loadcensus()                      // Census data is loaded to get the populations in regions
{
	string code, line, pop, x, y;
	long r;
	
	ifstream input("Data/censustable.txt");
	getline(input,line);
	do{
		getline(input,line);
		if(input.eof()) break;
		
		stringstream ss(line);	
		getline(ss,code,'\t');
		getline(ss,pop,'\t');
		getline(ss,x,'\t');
		getline(ss,y,'\t');

		for(r = 0; r < nregion; r++){ if(region[r].code == code) break;}
		if(r < nregion){
			region[r].pop = atoi(pop.c_str()); region[r].x = atof(x.c_str()); region[r].y = atof(y.c_str());
		}
	}while(1 == 1);
}

void createbici(string geo)                    // Creates a model file which can be loaded into BICI
{
	long r, rr, a, aa, flag;
	double fac;
	
	generateM();

	stringstream ss; ss << geo << "_model.txt";
	ofstream outp(ss.str().c_str());
		
	outp << "compartment classification='DS' name='S' x='0' y='0' color='#44ff44'\n";
	outp << "compartment classification='DS' name='E' x='1' y='0' color='#ffff44'\n";
	outp << "compartment classification='DS' name='A' x='2' y='0' color='#ff7777'\n";
	outp << "compartment classification='DS' name='I' x='3' y='0' color='#aa0000'\n";
	outp << "compartment classification='DS' name='H' x='4' y='0' color='#006600'\n";
	outp << "compartment classification='DS' name='R' x='3' y='-1' color='#4444ff'\n";
	outp << "compartment classification='DS' name='D' x='3' y='1' color='#444444'\n";
		 
	outp << "timetransitions values='83'\n";
		 
	outp << fixed << showpoint << setprecision(10);
	
	for(r = 0; r < nregion; r++){
		outp << "compartment classification='loc' name='" << region[r].code << "' x='" << region[r].x << "' y='" << -region[r].y << "' color='#444444'\n";
	}
	
	for(a = 0; a < nage; a++){
		outp << "compartment classification='age' name='" << age[a] << "' x='" << a << "' y='0' color='#44ffff'\n";
	}
	
	outp << "transition from='S' to='E' type='markovian' rate='\n";
	for(r = 0; r < nregion; r++){
		if(nage > 1){
			for(a = 0; a < nage; a++){
				outp << region[r].code << "," << age[a] << ":[Φ] + [β_Time]*(";
				outp << "(1-[d])*(";
				for(aa = 0; aa < nage; aa++){
					if(aa != 0) outp << "+";
					fac = mat[r][r]/region[r].pop;
					outp << 0.2*fac << "*{A," << region[r].code << "," << age[aa] << "}+" << fac << "*{I," << region[r].code << "," << age[aa] << "}";
				
					flag = 0;
					for(rr = 0; rr < nregion; rr++){
						if(rr != r){
							if(mat[r][rr] > 0.001){
								fac = mat[r][rr]/region[rr].pop;
								if(flag == 0){ outp << "+[d]*("; flag = 1;}
					      else outp << "+";
								outp << 0.2*fac << "*{A," << region[rr].code << "," << age[aa] << "}+" << fac << "*{I," << region[rr].code << "," << age[aa] << "}";
							}
						}
					}
					if(flag == 1) outp << ")";
					outp << ")";
				}
				outp << ")\n";
			}
		}
		else{
			outp << region[r].code << ":[Φ] + exp([G_" << region[r].code << "])*[β_Time]*(";
			outp << "(" << 1.0/region[r].pop << "-" << (1-mat[r][r])/region[r].pop << "*[d])*(";
			outp << 0.2 << "*{A," << region[r].code << "}+1*" << "{I," << region[r].code << "}";
			outp << ")";
			
			flag = 0;
			for(rr = 0; rr < nregion; rr++){
				if(rr != r){
					if(mat[r][rr] > 0.001){
						fac = mat[r][rr]/region[rr].pop;
						if(flag == 0){ outp << "+[d]*("; flag = 1;}
				    else outp << "+";
						outp << 0.2*fac << "*{A," << region[rr].code << "}+" << fac << "*{I," << region[rr].code << "}";
					}
				}
			}
			if(flag == 1) outp << ")";
			outp << ")\n";
		}
	}
	outp << "'\n";
	
	outp << "transition from='E' to='A' type='markovian' rate='[ρ(E→A)]'\n";
	outp << "transition from='A' to='I' type='markovian' rate='[ρ(A→I)]'\n";
	outp << "transition from='A' to='R' type='markovian' rate='[ρ(A→R)]'\n";
	outp << "transition from='I' to='R' type='markovian' rate='[ρ(I→R)]'\n";
	outp << "transition from='I' to='H' type='markovian' rate='[ρ(I→H)]'\n";
	outp << "transition from='I' to='D' type='markovian' rate='[ρ(I→D)]'\n";
	outp << "transition from='H' to='D' type='markovian' rate='[ρ(H→D)]'\n";
	outp << "transition from='H' to='R' type='markovian' rate='[ρ(H→R)]'\n";
	
	outp << "setprior param='Φ' prior='flat' min='0' max='0.000001'\n";
	outp << "setprior param='β_Time' prior='flat' min='0' max='3'\n";
	outp << "setprior param='d' prior='Fix' value='1'\n";
	outp << "setprior param='ρ(E→A)' prior='Fix' value='" << 1.0/tEA << "'\n";
	outp << "setprior param='ρ(A→I)' prior='Fix' value='" << 1.0/tAI << "'\n";
	outp << "setprior param='ρ(A→R)' prior='Fix' value='" << 1.0/tAR << "'\n";
	outp << "setprior param='ρ(I→R)' prior='Fix' value='" << 1.0/tIR << "'\n";
	outp << "setprior param='ρ(I→H)' prior='Fix' value='" << 1.0/tIH << "'\n";
	outp << "setprior param='ρ(I→D)' prior='Fix' value='" << 1.0/tID << "'\n";
	outp << "setprior param='ρ(H→D)' prior='Fix' value='" << 1.0/tHD << "'\n";
	outp << "setprior param='ρ(H→R)' prior='Fix' value='" << 1.0/tHR << "'\n";
	outp << "setdistribution param='G_loc' prior='Normal' mean='[μ]' sd='[σ]'\n";
	outp << "setprior param='σ' prior='flat' min='0.01' max='0.5'\n";
	outp << "setprior param='μ' prior='Fix' value='0'\n";
	outp << "inftimerange min='0' max='" << tmax << "'\n";
	
	for(r = 0; r < nregion; r++){
		outp << "addderived param='" << region[r].code << "-S' expression='{" << region[r].code << ",S}'\n";
		outp << "addderived param='" << region[r].code << "-E' expression='{" << region[r].code << ",E}'\n";
		outp << "addderived param='" << region[r].code << "-A' expression='{" << region[r].code << ",A}'\n";
		outp << "addderived param='" << region[r].code << "-I' expression='{" << region[r].code << ",I}'\n";
		outp << "addderived param='" << region[r].code << "-H' expression='{" << region[r].code << ",H}'\n";
		outp << "addderived param='" << region[r].code << "-R' expression='{" << region[r].code << ",R}'\n";
		outp << "addderived param='" << region[r].code << "-D' expression='{" << region[r].code << ",D}'\n";
	}
}

void generateM()                // Calculates the matrix which describes the mixing between regions
{
	long r, rr;
	double sum, d = 10, dd, d0;
	long BX, BY, P, mx, my, k, ii, jj, kk, p, kst;
	double x, y, xmin, xmax, ymin, ymax, max;
	vector< vector <double> > px;
	vector< vector <double> > py;
	vector< vector <short> > pr;
	
	mat.resize(nregion); 
	for(r = 0; r < nregion; r++){ for(rr = 0; rr < nregion; rr++) mat[r].push_back(0);}
	
	d0 = 30;
	for(r = 0; r < nregion; r++){
		for(rr = 0; rr < nregion; rr++){
			dd = (region[r].x - region[rr].x)*(region[r].x - region[rr].x) + 
					 (region[r].y - region[rr].y)*(region[r].y - region[rr].y); 
			if(dd < d0) dd = d0;
			
			if(r == rr) mat[r][rr] = 1;
			else mat[r][rr] = 0.1*exp(-dd/(d0*d0));
		}
	}
	
	vector <short> done; done.resize(nregion);
			
	for(r = 0; r < nregion; r++){
		for(rr = 0; rr < nregion; rr++) done[rr] = 0;
		
		for(rr = 0; rr < nregion; rr++){
			max = 0;
			for(k = 0; k < nregion; k++){ 
				if(done[k] == 0 && mat[r][k] > max){ max = mat[r][k]; kst = k;}
			}
			if(rr < 3 && max < 0.002) mat[r][kst] = 0.002;
			if(rr > 6) mat[r][kst] = 0;
			done[kst] = 1;
		}
	}
	
	for(r = 0; r < nregion; r++){                          // Normalises
		sum = 0; for(rr = 0; rr < nregion; rr++) sum += mat[r][rr];
		for(rr = 0; rr < nregion; rr++) mat[r][rr] /= sum;
	}
}

void generatexml(string geo)
{
	long a, r, pp, pop, jj, i;
	string line;
	double sumst[nage], Nsum[nage], sum, z;

	stringstream ss; ss << geo << "_bici_inputsmall.xml";
	ofstream bfile(ss.str().c_str());
	
	stringstream ss2; ss2 << geo << "_bici_model.xml";
	ifstream bhead(ss2.str().c_str());
	do{
		getline(bhead,line);
		if(bhead.eof()){ cout << ss2.str().c_str() << ": File not found\n"; break;}
		bfile << line << "\n";
		if(line == "</model>") break;
	}while(1 == 1);
	
	bfile << "<data>\n";
	
	bfile << "<capev name='capev0' type='trans' from='I' to='H' tmin='0' tmax='" << tmax << "' pd='1'/>\n";
	
	sum = 0; for(a = 0; a < nage; a++){ sum += N[a]/100.0; Nsum[a] = sum;}
	sum = 0; for(a = 0; a < nage; a++){ sum += N[a]*ph[a]; sumst[a] = sum;} 
  for(a = 0; a < nage; a++) sumst[a] /= sum; 
 
	pop = 0;
	for(r = 0; r < nregion; r++){
		cout << region[r].name << "   Total population:" << pop << " r\n";
		pp = region[r].pop; pp /= 10000;

		vector <vector <double> > agedist;
		agedist.resize(nage);
		for(jj = 0; jj < region[r].inf.size(); jj++){
			z = ran(); a = 0; while(a < nage && z > sumst[a]) a++;
			if(a == nage) cout << "pr\n";
			if(ran() < 0.001) agedist[a].push_back(region[r].inf[jj]);
			//agedist[a].push_back(region[r].inf[jj]);
		}
		//for(a = 0; a < nage; a++) cout << a << " " << agedist[a].size() << " agedist\n";
	
		if(nage > 1){
			/*
			for(a = 0; a < nage; a++){ 		
				bfile << "<setallinit id='Ind' min='"<< pop << "' max='" << pop + long(sumst[a]*region[r].inf.size()) << " DS='S' loc='" << region[r].code << "' age='" << age[a] << "' Age='A0+' Time='All' />\n";
			}
			*/
		}
		else{
			bfile << "<setallinit id='Ind' min='"<< pop << "' max='" << pop+pp-1 << "' DS='S' loc='" << region[r].code << "' age='" << age[a] << "' Age='A0+' Time='T<83' />\n";
		}
		
		jj = 0;
		for(i = 0; i < pp; i++){
			stringstream ss; ss << "Ind" << pop;
	
			z = i/pp; a = 0; while(a < nage && z > Nsum[a]) a++;
			if(a == nage) cout << "prob\n";
			
			if(agedist[a].size() > 0){
				bfile << "<transition id='" << ss.str() << "' time='" << agedist[a][agedist[a].size()-1] << "' capev='capev0' from='I' to='H'/>\n";
				agedist[a].pop_back();
			}
			pop++;
		}
	}
	
	bfile << "</data>\n";
	bfile << "<inference tmin='0' tmax='" << tmax << "' indmax='70000000'/>\n";
	bfile << "</BICI>\n";

	cout << pop << " " << "Done!\n";	
}

