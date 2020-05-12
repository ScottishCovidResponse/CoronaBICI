/*
Load mpi: module load mpi/openmpi-x86_64
Compile using:    mpic++ bici.cc header/tinyxml2.cc -O3 -o bici

Run using:        mpirun -n 4 ./bici Scotland_bici_input.xml 1000000

-n 4 gives the number nodes (i.e. MCMC chains)
"Scotland_bici_input.xml" is the model/data file to be analysed
1000000 is the number of MCMC iterations

Once run the software will initialise and show the samples incrementaly increase
The following outputs are placed into the Outputs/<model name> directory:


a) "Trace[X].txt" (where X is the number of the chain) gives trace plots of the model parameters. 
   These can be loaded into other software to visually check that MCMC is mixing well (e.g. Tracer).

b) "Diagnostic[X].txt" gives information about the MCMC proposals to help identify inefficiencies in the 
   code.

c) "Bici[X].txt" is an output file which can be read into the BICI GUI for visualisation.

d) "Stats.txt" calculates the posterior means, 90% credible intervals, effective sample sizes and 
   Gelman-Rubin diagnostic statics (these last two are used to confirm that MCMC is well mixed).
*/

// In order to use MPI, compile with -DUSE_MPI

const long noout = 0;                 // Suppresses the output when error checking
const long checkon = 0;               // Determines if checking is done
long simnum;                          // Store the MCMC chain number

#include "header/tinyxml2.h"

using namespace tinyxml2;

using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

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
                                  // Most of the code is placed into header files:
#include "header/var.h" 					// Information about the global variables and constants                  
#include "header/equation.h"			// Converts strings representing equations into calculations
#include "header/readinput.h"			// Reads in the input XML file and sets up the model
#include "header/incompgamma.h"		// Calculates the incomplete gamma function
#include "header/dist.h"					// Samples from or calcualtes probabilities for distributions
#include "header/sim.h"						// Simlates from the model
#include "header/output.h"        // Generates output files
#include "header/derive.h"				// Deals with plotting derived quantities
#include "header/init.h"          // Initialises the MCMC procedure
#include "header/startseq.h"			// Generates a pausible initial event sequence consistent with data
#include "header/prior.h"         // Calculates prior probabilities
#include "header/likelihoodinit.h"// Initialises fast calculation of the likelihood
#include "header/likelihood.h"    // Calculates changes to likelihood under a change to an individual
#include "header/check.h"         // Used to check every is working correctly
#include "header/ratecalc.h"      // Calculates the rates of transitions
#include "header/observation.h"   // Deals with state observations on the model
#include "header/param_prop.h"    // Makes parameter MCMC proposals 
#include "header/part_propnew.h"  // Makes particle proposals for individual event sequences
#include "header/local_prop.h"		// Makes local changes to individual event sequences
#include "header/move_prop.h"			// Considers moving individual event time proposals
#include "header/multimove.h"     // Considers making multiple event time proposals
#include "header/life_prop.h"     // Proposals affecting the lifespan of an individual
#include "header/addrem_prop.h"		// Proposals to add and remove individuals from the system
#include "header/indsim_prop.h"		// Proposals which simulate event sequences for individuals
#include "header/addreminf.h"     // Proposals for adding and removing infected individuals (Corona)
#include "header/samplers.h"      // Different event samplers used in proposals
#include "header/singlefixed_prop.h"   // Proposals for individuals with a single fixed event specified

#include "header/stats.h"         // Generate statistics to be output at the end of execution 

int main(int argc, char** argv)
{
  string file;
  vector <double> v;
 
	#ifdef USE_MPI
	int process_Rank, size_Of_Cluster;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
	simnum = process_Rank;
	nrun = size_Of_Cluster;
  #endif
 
  #ifndef USE_MPI
	simnum = 0;
	nrun = 1;
  #endif

	if (simnum == 0) {
		cout << "CoronaBICI running ";
#ifdef USE_MPI
		cout << "on " << size_Of_Cluster << " MPI processes" << endl;
#else
		cout << "in serial" << endl;
#endif
	}
 
	if(argc != 3) emsg("Not the right number of arguments");
	file = argv[1]; 
	nsamp = atoi(argv[2]); burnin = nsamp/3;
	
	stringstream sd; sd << "Output/" << file.substr(0,file.length()-10);
	root = sd.str();

	if (simnum == 0) {
		ensuredirectory("Output");
		ensuredirectory(root);
	}
	
	// Ensure directories are created before other processes proceed
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	stringstream st; st << root << "/trace" << simnum << ".txt";
	trace.open(st.str().c_str());
	
	stringstream sb; sb << root << "/bici" << simnum << ".txt";
	bout.open(sb.str().c_str());
	
	srand(simnum);
	
  bout << "3|0|\n"; bout.flush();

	if(simnum == 0) cout << "Reading input file...\n";
  
  readinput(file); bout << "3|10|\n"; bout.flush();

	if(simnum == 0) cout << "Initialising...\n";
	
  init(); bout << "3|15|\n"; bout.flush();

  likelihoodinit(); bout << "3|20|\n"; bout.flush();

  if(derfl == 1) deriveinit();
  priorinit(); bout << "3|30|\n"; bout.flush();
  observationinit(); bout << "3|40|\n"; bout.flush();
  paraminit();
  partinit();

  localpropinit();
  lifepropinit(); bout << "3|50|\n"; bout.flush();
  moveinit();

  if(simon == 0) startseqinit();

  bout << "3|60|\n"; bout.flush();

  ch[0] = new Chain;
  nchain = 1;

  ch[0]->invT_pd = 1;
  ch[0]->invTLi = 1;
  ch[0]->invT_st = 1;
  ch[0]->invT_pop = 20;
  ch[0]->invT_der = 10;

  siminit();

  ch[0]->initparamsamp();

  ch[0]->initchain(); 

  if(simon == 0){
		if(corona == 1) ch[0]->addreminfinit();
		if(corona == 1) ch[0]->samplerinit();
		ch[0]->start();
	}
 
  bout << "3|90|\n"; bout.flush();
  //
	if(corona == 0) ch[0]->addreminit(); 
	
  bout << "3|95|\n"; bout.flush();

  if(nderive > 0) ch[0]->deriveplotinit();

  bout << "3|100|\n"; bout.flush();

  traceinit();
	
	/*
	ch[0]->eventplot();
  ch[0]->traceplot();
	return 0;
	*/
	
  if(simon == 1){
		bout << "8|\n"; bout << "3|0|\n"; bout.flush();
	 for(samp = 0; samp < nsamp; samp++){ 
      ch[0]->eventplot();
		  ch[0]->traceplot();
    }
		
		if(corona == 1){
			ch[0]->outputcoronadata();
			ch[0]->savesimulated();
		}
		
		#ifdef USE_MPI
		MPI_Finalize();
		#endif
    return 0;
  }
	
	if(corona == 1)	ch[0]->initinflist();
	ch[0]->multimoveinit();
	
	//if(corona == 1) ch[0]->addstartuoinf();
	ch[0]->paramstart(); 
	
  ch[0]->check(-1);

	tracefileinit();

  for(samp = 0; samp < nsamp; samp++){
		ch[0]->update();

    altertemperature();

    if(samp%100 == 0 && simnum == 0) cout << "Sample: " << samp << " / " << nsamp << "\n";
	
    if(samp%(nsamp/1000) == 0){ ch[0]->traceplot(); ch[0]->tracefile();}
		if(samp >= burnin && samp%(nsamp/20) == 0) ch[0]->eventplot();
		if(samp >= burnin && samp%(nsamp/1000) == 0) ch[0]->addparam();
    
    if(samp != 0 && samp%(nsamp/10) == 0){
			if(corona == 0) ch[0]->diagnosticschain();
			ch[0]->diagnosticsfile();
			ch[0]->check(0);
		}
	}
	if(simnum == 0) cout << "\nSee files in '" << root << "' for outputs.\n";
	
	ch[0]->diagnosticsfile();
	stats();
	
	#ifdef USE_MPI
	MPI_Finalize();
	#endif
}

void Chain::update()                                                // Performs an MCMC update
{
  totaltime -= clock();
	 
  timeprop[PARAM_PROP] -= clock();
	param_prop();
	//for(loop = 0; loop < paramloop; loop++) paramnorm_prop();
  timeprop[PARAM_PROP] += clock();
  ntimeprop[PARAM_PROP]++;
	
  timeprop[ADDREMINF_PROP] -= clock();
 	addreminf();
	timeprop[ADDREMINF_PROP] += clock();

	timeprop[MULTIMOVE_PROP] -= clock();
	//if(ran() < 0.02) multimove();
	timeprop[MULTIMOVE_PROP] += clock();
	
	timeprop[SINGFIX_PROP] -= clock();
	//singlefixed_prop();
	singlefixed_prop2();
	timeprop[SINGFIX_PROP] += clock();
	
	// These MCMC proposals are currently turned off as not needed for analysing Corona data (yet)
	/*  
	 long chon = 0, paramloop, loop, loopmax, i, imax, r;
  double pr;
	
  //if(samp == 99) chon = 1;

  if(samp < burnin && samp%5 == 0) optimiseproposals();
  if(samp < burnin && samp%10 == 0) simsumcalc();

  if(totaltime == 0) paramloop = 4;
	
  else{ paramloop = long(5-8*timeprop[PARAM_PROP]/totaltime); if(paramloop < 1) paramloop = 1;}
 
	if(ran() < part_prob){
    timeprop[PART_PROP] -= clock();
    for(i = 0; i < nind; i++) part_prop(i);
    timeprop[PART_PROP] += clock();
    ntimeprop[PART_PROP]++;
  }

  //if(chon == 1) check(2);

  if(nindtot > nind){
    if(ran() < induo_prob){
      timeprop[INDSIMUO_PROP] -= clock();
      for(i = nind; i < nindtot; i++) indsim_prop(i);
      timeprop[INDSIMUO_PROP] += clock();
      ntimeprop[INDSIMUO_PROP]++;
    }
  }
	
  if(ran() < move_prob){
    timeprop[MOVE_PROP] -= clock();
    for(i = 0; i < nind; i++){
			if(i%1 == 0) cout << i << "i\n";
			move_prop(i);
		}
    timeprop[MOVE_PROP] += clock();
    ntimeprop[MOVE_PROP]++;
  }

  if(chon == 1) check(4);

  if(singeventfl == 1){
    if(ran() < sing_prob){
      timeprop[SINGEVENT_PROP] -= clock();
      for(i = 0; i < nind; i++) sing_prop(i);
      timeprop[SINGEVENT_PROP] += clock();
      ntimeprop[SINGEVENT_PROP]++;
    }
  }

  if(chon == 1) check(5);

  if(twothreefl == 1){
    if(ran() < twothree_prob){
      timeprop[TWOTHREE_PROP] -= clock();
      for(i = 0; i < nind; i++) twothree_prop(i);
      timeprop[TWOTHREE_PROP] += clock();
      ntimeprop[TWOTHREE_PROP]++;
    }
  }

  if(chon == 1) check(6);
	
  if(pairfl == 1){
    if(ran() < pair_prob){
      timeprop[PAIR_PROP] -= clock();
      for(i = 0; i < nind; i++) pair_prop(i);
      timeprop[PAIR_PROP] += clock();
      ntimeprop[PAIR_PROP]++;
    }
  }

  if(chon == 1) check(7);

  if(gapfl == 1){
    if(ran() < gap_prob){
      timeprop[GAP_PROP] -= clock();
      for(i = 0; i < nind; i++) gap_prop(i);
      timeprop[GAP_PROP] += clock();
      ntimeprop[GAP_PROP]++;
    }
  }

  if(chon == 1) check(8);

  if(tbirthfl == 1 || sourcefl == 1 || sinkfl == 1){
    if(ran() < life_prob){
      timeprop[LIFE_PROP] -= clock();
      for(i = 0; i < nindtot; i++) life_prop(i);
      timeprop[LIFE_PROP] += clock();
      ntimeprop[LIFE_PROP]++;
    }
  }

  if(chon == 1) check(9);

  if(addremfl == 1){
    if(ran() < addrem_prob){
      timeprop[ADDREM_PROP] -= clock();

      imax = (nind+nindtot)/4; if(imax == 0) imax = 1;
      for(i = 0; i < imax; i++) addrem_prop();
      timeprop[ADDREM_PROP] += clock();
      ntimeprop[ADDREM_PROP]++;
    }
  }

  if(chon == 1) check(10);
	*/

  totaltime += clock();
}

