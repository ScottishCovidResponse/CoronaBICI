// This file gives information about all the global variables and constants

long corona = 1;                        // This is set to 1 for analysing the corona virus data

const double gammalim = 0.05;                 // Limit on gamma / beta distribution shape parameter

const long discon = 1;                        // Set to 1 if dependent equations are discretised in time
const long DX = 400;                          // Number of discretisations 
double dfac;	                                // Factor converting from t to discretisation

ofstream bout, trace;                         // Used for outputting bici and trace files
string root;                                  // The directory for output
vector< vector < vector <double> > > paramst; // Store parameter (for use in stats.h)

long samp, nsamp;                             // The MCMC sample number
long burnin;                                  // The MCMC number of samples for burnin 

// CORONA SPECIFIC
long tbin = 100;                              // Temporal discretation when adding/removing inf ind
vector <double> multacf;                      // Number to add/rem from a region
vector <short> indinit; 											// This sets the initial state of an individual

const long SS = 0, EE = 1, AA = 2, II = 3;    // Used to specify states
const long HH = 4, RR = 5, DD = 6;

vector <vector <double> > prob;               // Used to sample infection time
vector <vector <double> > probsum;
vector <double> proba;
vector <double> probsuma;
double probsumatot;
double rEA, rAR, rAI, rIR[7], rIH[7], rID[7], rHR[7], rHD[7];
long nregion;                                 // Number of regions
long nag;                                     // Number of age groups in model
long nrun;                                    // The number of MCMC runs

// END CORONA SPECIFIC

long plfl = 0;                                // Determines if diagnostic information is output

// Codes used form timings of different procedures:
const long PARAM_PROP = 0, PART_PROP = 1, MOVE_PROP = 2, LIFE_PROP = 5, ADDREM_PROP = 6;
const long SINGEVENT_PROP = 7, TWOTHREE_PROP = 8, INDSIMUO_PROP = 10, PAIR_PROP = 11, GAP_PROP = 12;
const long CHECK = 13, INDCHANGE = 14, TRACEPLOT = 15, EVENTPLOT = 16, PARTSIM = 17, PARTINIT = 18;
const long INITBA = 19, INIT = 20, ADDREMINF_PROP = 21, MULTIMOVE_PROP = 22, SEC = 23, SAMPPROB = 24;
const long SIMTIME = 25;

// Different types of transitions:
const long EXP_TR = 0, FIXED_TR = 1, GAMMA_TR = 2, WEI_TR = 3, GROW_TR = 4, SETTIME_TR = 5;
const long CAPEVTRANGE_TR = 6, BEG_TR = 7, END_TR = 8, NULL_TR = 9;

// Different types of priors:
const long FLAT = 0, GAMMA = 1, NORMAL = 2, LOGNORMAL = 3, EXPO = 4, BETA = 5, WEIBULL = 6, FIX = 7;
const long SMOOTH = 0, LOGSMOOTH = 1;

// Different types of observation
const long CAP = 0, POP = 1, DER = 2;

long DERIVEX = 400;                           // Temporal discretisation of derived quantities

long indmax;                                  // The maximum number of individuals
const double large = 10000000;                // Used to refer to a large quantity
const double mlarge = -large;                 // A large negative quantity
const double tiny = 0.000000001;              // A tiny quantity
const double ttiny = 0.000000000001;          // A very tiny quantity

const double evdtmin = 0.00000001;            // The smallest seperation time between events
long logsummax = 10000;                       // The size of the lookup table for logsum

const long depeqdivmax = 1000;                // Temporal discretisation of lookup for dependent equations
const long ncompcapdiv = 1000;								// Temporal discretisation of lookup for captures
const double ratesmall = 0.0000000001;        // The smallest rate size
const double obssmall = 0.00001;     				 	// Give a minimum limit on the obsevation probability
const double notobsdL = -1000;           			// Give the penalty for a non-observation
const double onefac = 0.99999999;             // A factor just less than 1
const long npartmax = 40;                     // The maximum particle number
const long npartstart = 10;                   // The number of particles to start with

double muest;                                 // Estimate for mortality rate 
double capprobest;                            // Estimate for the capture probability

long loopstartmax = 1000;                     // Max loop when finding a valid initial starting seqeucne 

double totaltime = 0;                         // Total execution time
double timeprop[30], ntimeprop[30];           // The CPU timings for the various proposals

long fixfl=0;                                 // Set to one if any events are obseved in the data
long sinkfixfl=0;                             // Set to one if any fixed sink transitions in the model
long sourcefixfl=0;                           // Set to one if any ficed source transition in the model
long movefl=0;                                // Set to one if the data contains move, sourcemove, sinkmove
long sinkfl=0, sourcefl=0;                    // Set to one if the model contains sinks or sources
long gammafl=0;                               // Set to one if any gamma transitions in model
long nmfl=0;                                  // Set to one if any non-Markovian transitions in model
long capevfl = 0;                             // Set to one if there is a probability of captures
long capfl=0;                                 // Set to one if there are any captures
long popfl=0;                                 // Set to one if there are any population measurements
long derfl=0;                                 // Set to one if there are any derived measurements
long tbirthfl=0;                              // Set to one if there are more than one age classification
long pairfl=0;                                // Set to one if zero/two event proposals used
long singeventfl=0;                           // Set to one if single events can be inserted and deleted
long twothreefl=0;                            // Set to one if interchange between two and three events
long gapfl=0;                                 // Set to one if single events can be inserted and deleted
long addremfl=0;                              // Set to one if there are unobserved individuals 

long simon=0;                                 // Set to one if a simulation in run

// Model parameters

long agecl, nage;                              // Characterises the age time line
vector<double> age;

long settimecl, nsettime;                      // Characterises the Time time line
vector<double> settime;

long nclass;                                   // # of classifications
vector<string> classname;                      // Names of classifications
vector<long> nclassval;                        // # of compartments in a classification
vector< vector<string> > classval;             // Names of compartments
vector<long> classmult;

vector<string> compname;                       // Overall names of compartments 

vector<long> paramorder;                       // Orders the parameters for output

struct TT {                                    // Structure used for transitions
  long type; long ci; long cf; long dc; long cl; long i; long f; long eq; long eqshape;
  long capev; long like; long trans; long ntraend; vector<long> traend; long nm;
  long nnonexptra; 
	vector<long> nonexptra;// The nonexponential transitions activated when individual undergoes transition
};

long ntra;                                     // Total number of transitions going between compartments
vector<TT> tra;                                // Stores all possible transitions in the model
long moventra=0;                               // Denotes which transitions are used for move events
long trabeg;                                   // Denotes which transitions are use to begin the timeline
long traend;                                   // Denotes which transitions are use to end the timeline
long tranull;                                  // Denotes a null transition

long ncomp;                                    // The total number of compartents
long ncomps;                                   // Does not include age, time and fix timelines
long ncompswa;                                 // Does not include time and fix timelines
long ncompswf;                                 // Does not include fix timelines
vector< vector <long> > compiftra;             // [ci][cf] Gets the tra number
vector< vector <long> > compval;               // [c][cl] Gets the class value for a particular compartment

vector <long> ncompleave;                      // Gives the tra leaving a compartment
vector < vector <long> > compleave;

vector< vector <long> > ncompclleave;          // Gives the tra leaving a compartment with a particular cl
vector< vector < vector <long> > > compclleave;

vector < vector < vector <long> > > compleavesimcldep;    // Dep tras leaving a compartment for sim
vector < vector < vector <long> > > compleavesimclnotdep; // Not dep tras leaving a compartment for sim

vector <long> ncompleavedep;                   // Gives dep tra leaving a compartment 
vector < vector <long> > compleavedep;
vector <long> ncompleavenotdep;                // Gives notdep tra leaving a compartment 
vector < vector <long> > compleavenotdep;

// Prior

long nprior;                                   // The number of priors
vector<long> priortype;                        // The type of prior
vector<long> priorparam;                       // The parameters on which the prior is placed
vector<double> priorminval;                    // The min and max values for a flat prior
vector<double> priormaxval;
vector<long> prioreq1;                         // Eqs can be used to define priors defining quantitives
vector<long> prioreq2;

vector<double> agprior;                        // The prior used for the birth time

vector <double> dirchal;                       // The Dirchlet prior placed on the initial composition

// Data

double tmax;                                   // The range over which inference is performed
double tmax2;                                  // Additional time used for furture prediction
double tminactual, tmaxactual;                 // The total range over which outputs are given

long nind;                                     // The total number of observed individuals
vector<string> indid;                          // Individual ids
vector<short> nindobs;                         // Individual observations
vector< vector<short> > indobs;
vector< vector<double> > indobst;
vector<long> nindfixev;                        // All fixed events associated with an individual
vector< vector<long> > indfixev;
vector<double> indfixtenter;                   // The entry time of the individual (if set otherwise large)
vector<double> indfixtleave;                   // The leave time of the individual (if set otherwise large)
vector<short> indlikeenter;                    // Determines if the enter likelihood is evaluated
vector<short> indlikeleave;                    // Determines if the leave likelihood is evaluated 
vector<short> indfixentercapev;                // The capev reference when entering
vector<short> indfixleavecapev;                // The capev reference when leaving

long nobs;                                     // The state observations
vector<double> obst;
vector<long> obsi;
vector<long> obscap;                           // The capture corrspending to a particular observation
vector< vector<long> > obsprobeqn;             // [#obs][c] Prob of observation in each compartment

long ncope;                                    // This gives a list of all the eqations used in obsprobeqn
vector<long> cope;
vector< vector<long> > obsprob_cope;

struct FEV {                                   // Structure used for fixed events
	long i; double t; short cl; short trai; short traf; short dc; short capev; short type; short like;
};
bool compareFEV(FEV lhs, FEV rhs) { return lhs.t < rhs.t; }

long nfixev;                                   // Stores all the fixed events given in the data
vector<FEV> fixev;

long ncapev;                                   // Event captures
vector<long> capevall;                         // Set if all all events are captured
vector<long> capevprobeqn;                     // The equation for the capture event probability

long capevcl, ncapevtrange;                    // Characterises the capture events time line
vector<double> capevtrange;
vector<string> capevname;
vector<double> capevtmin;
vector<double> capevtmax;
vector<long> capevtype;
vector< vector<long> > capevfilt;
const short CESOURCE = 0, CESINK = 1, CETRANS = 2;

long ncap;                                     // The capture campaigns
vector<string> capname;
vector<double> capt;
vector< vector<long> > capprobeqn;             // [#cap][c] Prob for capturing in each of the compartments

long ncpe;                                     // Gives a list of equations used in capprobeqn
vector<long> cpe;
vector< vector<long> > capprob_cpe;

struct CC { long val; long ty; double t;};     // Structure used for capture campaigns
bool comparecc(CC lhs, CC rhs) { return lhs.t < rhs.t; }
vector<long> ncompcap;                         // Gives the captures which occur in a given compartment
vector< vector<CC> > compcap;
vector< vector<long> > compcapdiv;

vector< vector<short> > capindob;              // Obs number for a given individual in a given capture

long npopm;                                    // Population measurements
vector<double> popmt;                          // The times of the measurements
vector<double> popmal;                         // Measurements gamma distributed with an alpha and beta
vector<double> popmbe;
vector< vector<long> > popmcomp;               // The compartments that make up the population
vector <long> popmap;                          // Used to make changes in the observation probability 
vector <long> poplist;

long nderm;                                    // Derived quantity measurements
vector<long> derm;                             // Which quantity
vector<double> dermt;                          // The times of the measurements
vector<double> dermav;                         // Measurements are gamma distributed with alpha and beta
vector<double> dermvar;
vector< vector<long> > dermcomp;               // The compartments which affect derivation
vector <long> dermap;                          // Used to make changes in obs prob for derived measurements
vector <long> derlist;

long NOTALIVE;                                 // Represents the not alive compartment

long nparam;                                   // Number of parameters in the model
vector<string> paramname;                      // The name of the parameter
vector<string> paramclass;                     // A reference names used in BIP to classify traces
vector<long> paramprior;                       // The prior used for the prior
vector<long> nparampriordep;                   // the priors which depend on that parameter
vector< vector<long> > parampriordep;

long nderive;                                  // The number of derived quantities
vector<string> derivename;
vector<long> derive;
vector< vector< long> > nderive_dpop;
vector< vector< vector<long> > > derive_dpop;  // Works out how populations change as compartments change
vector< vector< vector<double> > > derive_dpopweight;

long nsmooth;                                  // Time smoothing priors
vector<string> smoothname;
vector<long> smoothtype;
vector<string> smoothref;
vector<long> nsmoothparam;
vector< vector<double> > smoothval;
vector< vector<long> > smoothparam;

vector< vector<long> > paramsmooth;
vector< vector<long> > paramsmoothi;

long nderivetemporal;                          // Lists all derived quantities that are temporal
vector <long> derivetemporal;

vector<long> nparam_dep;                       // Which dep equations a parameter is in
vector< vector<long> > param_dep;

vector<long> nparam_notdep;                    // Which not dep equations a parameter is in
vector< vector<long> > param_notdep;

vector<long> nparam_cope;                      // Which observation probabilities a parameter is in
vector< vector<long> > param_cope;

vector<long> nparam_cpe;                       // Which detection probabilities a parameter is in
vector< vector<long> > param_cpe;

vector<long> nparam_capev;                     // Which capture events a parameter is in
vector< vector<long> > param_capev;

vector<long> nparam_derm;                      // Which derived quantities a parameter is in
vector< vector<long> > param_derm;

vector<long> param_nmfl;                       // Determines if change in param requires the nm update 

vector <double> firstobst, lastobst;           // The first and last time an individual is observed

// Uses in MCMC proposals

double problifeprop[7] = {0.25,0.25,0.1,0.1,0.1,0.1,0.1}; // Determines the frequecy of life proposals
double problifepropsum[7];

long nclnm;                                    // Quantities related to non-Markovian transitions
vector <long> clnm;
vector <long> cstnm;
vector <double> tstnm;

struct EV { double t; short tr;};              // Structure for events
bool compareev(EV lhs, EV rhs) { return lhs.t < rhs.t; }

struct PART {                                  // Structure for particles
	vector <EV> ev; long cend; double w; double wsum; long pback; long cstart; vector <long> futerefst;
};
vector<PART> part;                             // Stores particles

vector<long> neq_popnum;                       // List of all the populations involved in an equation
vector< vector<long> > eq_popnum; 
vector<long> neq_param;                        // Llist of all the parameters involved in an equation
vector< vector<long> > eq_param;

long npopnum;                                  // The number of populations which are kept track of
vector<string> popnumname;                     // Name of population
vector< vector<long> > popnumterm;             // The compartment populations which are added 
vector< vector<float> > popnumtermweight;      // The weightings given to compartments
long popnummax;                                // Maximum number of populatins an equation can depend on

long nnmeq;                                    // List of equation used in non-Markovian transitions
vector<long> nmeq;

vector<long> transdep;                         // [eqn]   -1 = not an exponential transitions  
// 0 = does not depend on population number   1 = depends on population number

vector<long> transdepref;                      // [eqn]  references which transdepeq
long ntransdepeq, ntransnotdepeq;
vector<long> transdepeq;                       // List of all the equations dependent on population
vector<long> transnotdepeq;                    // List of all the equations independent of population

vector< vector <long> > transdepeqrecalc;      // Determines if an eq needs to be recalculated
vector< vector <double> > popnumtake;

vector< vector< vector<long> > > transchref;



struct EQCH {                                  // Structure giving the change in a dep eqn
	long d; long n; long valch; double* popnum;
};
vector< vector<EQCH> > transdepeqch;

struct REFL { double t; double ref; long j; short sign;};
bool comparerefl(REFL lhs, REFL rhs) { return lhs.t < rhs.t; }
vector <long> depmap;                          // Used for making multiple changes to individuals
vector <long> depmaplist;
double *ch_popn;
vector <vector <REFL> >  depmaprefl;
long multich = 0;

struct NDEQCH { long d; long n;};
vector< vector<NDEQCH> > transnotdepeqch;

struct TRANS{ long cl; long i; long f;};       // Store all unique transitions in move (used for jump move)
long ntrans;                                   // Unique transitions
vector<TRANS> trans;

// Used for particles

const long CHC = 0, BEG = 2, END = 3, OBS = 4, FIXEV = 5, OTHCL = 6;
struct EVPART { double t; long ty; long ch;};

long nevp;                                     // When simulating from particles this shows the fixed times
vector <EVPART> evp;
vector<long> evpref;                           // Gives a reference to cexist

long irev;                                     // Used to keep track of event on the exisiting sequence
long erev, nevrev;
vector <EV> evrev;

vector <long> psel;

double wsum;
vector <double> wsumst;

long pbegst, pendst;                           // Temporarily stores variables

long nexist;
vector<long> cexist;                           // Keeps track of changes to c in the existing event seq
vector<double> texist;

long nRdep, nRnotdep;                          // Stores transition rates when performing simulation
double Rdep,  Rnotdep;
vector <double>  Rdepst, Rnotdepst, rdep, rnotdep;

bool evpcomp(EVPART lhs, EVPART rhs) { return lhs.t < rhs.t; }

long nfuteref;                                 // Stores future events
vector <long> futeref;

long npart;                                    // The number of particles
long partcl;  // Determines if transitions are performed along a particular classification or all
long clall;   // The classification used to represent all

struct EVPSEC{ long evpbeg; long evpend; double tbeg; double tend;};

long nevpsec;
vector<EVPSEC> evpsec ;                        // Divides up into sections which are simulated

double partdt;                                 // The m time step used when doing simulation

double tent, tlea;                             // Stores the entry and leaving time for an individual

long ctop;                                     // This age time and fix part of the compartment on entry

vector <vector <vector <vector< long> > > > fixnmmap;

// Used in adding and removing individuals

const long naddremdiv = 10;                    // Number of time divisions used when inserting individual

vector <double> induoagesumtot;
vector <vector <double> > loginduoage, induoage;// The age distribution for added individuals
vector <vector <double> > induoagesum;         // The age distribution for added individuals

double addremdivsumtot;
vector <double> logaddremdiv, addremdiv;
vector <double> addremdivsum;

vector <double> induolifetime;                 // The expected lifetime for added individuals

vector <double> logsum;                        // Stores the log of a number factorial
vector <EV> evnew;                             // Use to generate a proposed event sequence

short notimerange = 0;                         // Set to one if the data takes up no time

double pmean, pvar, pgrad;                     // Used when calculating the grapdient in the likelihood
vector<double> pvel;

vector <string> warning;                       // Keeps a list of any warning messages

// Temporary variables 

long ist, nwrong, n_popnum, sampill;
double Linmst, uocapmin, uocapmax;
string capwarn;

class Chain                                    // Store all quantities and funcion for an MCMC chain
{
  public:
    double Lir;                                // Likelihood of the events                
    double Liexp;                              // Likelihood of the gaps between events
    double Linm;                               // Likelihood for non-Markovian transitions
    double Liinit;                             // Likelihood of initial state
    double Lev_pd;                             // Likelihood from probability of detection
    double Lob_pd;                             // Observation probability from probability of detection
    double Lob_st;                             // Observation probability from observed states
    double Lpop;                               // Likelihood for the population measurements
    double Lder;                               // Likelihood for the derived measurements
    double Lpri;                               // The prior

    vector <long> pop;                         // Current population at an observe population point
    vector <double> popL;             				 // The likelihood associated with measurement

    vector<double*> dermpopnum;                // Populations used to calculate derived values
    vector <double> derL;                      // The likelihood associated with derived values

    double invT_pd;                            // Inverse temperature applied to the likelihoods
    double invT_st;
    double invTLi;
    double invT_pop;
    double invT_der;

    vector< vector <float> > nindpart;         // Number of particles used for an individual

    long nindtot;                              // Information about the individuals
    vector <double> indtbirth;
    vector<short> nindev;
    vector <vector <EV> > indev;

    double *param;                             // The parameter values

    long nderivetemporalplot;                  // Used for plotting the derived quantities
    vector< vector< vector<double> > > derivetemporalplot;
 
    vector< vector <double> > pairac;          // Used to optimise zero two pairs

		// Diagostic variables used to assess acceptance probabities for different proposals:
		double ntr_multinfadd, nac_multinfadd, ntr_multinfrem, nac_multinfrem;
		
    vector <float> ntr_smooth, nac_smooth, jump_smooth;
    vector <float> ntr_param, nac_param, jump_param;
    vector <float> ntr_paramsamp, nac_paramsamp;
    vector <float> ntr_paramnorm, nac_paramnorm;
    vector < vector <float> > ntr_part, nac_part, nfa_part, nde_part;

    vector <float> ntr_indsim, nac_indsim, nfa_indsim;
    float ntr_indsimuo, nac_indsimuo, nfa_indsimuo;

    vector <float> ntr_tbirth, nac_tbirth, jump_tbirth; 
    vector <float> ntr_tbirthent, nac_tbirthent, jump_tbirthent;
    vector <float> ntr_tent, nac_tent, jump_tent;
    vector <float> ntr_tlea, nac_tlea, jump_tlea;
    vector <float> ntr_entswitch, nac_entswitch;
    vector <float> ntr_birthentswitch, nac_birthentswitch;
    vector <float> ntr_leaswitch, nac_leaswitch;

    float ntr_tbirthuo, nac_tbirthuo, jump_tbirthuo;  
    float ntr_tbirthentuo, nac_tbirthentuo, jump_tbirthentuo;
    float ntr_tentuo, nac_tentuo, jump_tentuo;
    float ntr_tleauo, nac_tleauo, jump_tleauo;
    float ntr_entswitchuo, nac_entswitchuo;
    float ntr_birthentswitchuo, nac_birthentswitchuo;
    float ntr_leaswitchuo, nac_leaswitchuo;

    float ntr_add, nac_add, nfa_add;
    float ntr_rem, nac_rem;

    vector< vector <float> > ntr_move, nac_move;
    vector <vector <float> > jump_move;

		vector< vector <float> > ntr_movecomp, nac_movecomp;
    vector <vector <float> > jump_movecomp;
		
    vector< vector <float> > ntr_sing, nac_sing, nfa_sing;
    vector< vector <float> > ntr_twothree, nac_twothree, nfa_twothree;
    vector< vector <float> > ntr_pair, nac_pair, nfa_pair;
    vector< vector <float> > ntr_gap, nac_gap, nfa_gap;

    long nindinit;                            // Number of initial individuals
    vector <long> ncompinit;                  // Number of inds entering population in compartments
    vector <double> probinit;                 // Prob of an initial individual having a given compartment

		vector<long> depeq_disc_evn;              // The discretised time lines for dependent equations
		vector<double> depeq_disc_evdt;
		vector<double> depeq_disc_evval;
		vector<double*> depeq_disc_evpopnum;  
		
    vector< vector<long> > depeqdiv;          // Hash table in time

    long ndepeq_ev;         			            // Unorder list for time variation in the dependent equations
    vector<long> depeq_evn;             			// The number of times this equation appears in R
    vector< double* > depeq_evpopnum;   		 	// The population numbers for that particular
    vector<double> depeq_evt;            			// The time of the change in the dep eq 
    vector<long> depeq_evnext;           			// The next member on the list
    vector<long> depeq_evback;           			// The last member on the list
    vector<long> ndepeq_ev_ev;           			// The cases where the actual event happens
    vector< vector<double> > depeq_ev_evt;
    vector<double> depeq_evval;          			// The value for the equation

    long ndepeq_evstack;                      // Stores unused values to be dynamically added
    vector<long> depeq_evstack;

    vector<long> transnotdepeq_num;           // Quantities for the not dep equations
    vector<double> transnotdepeq_dt;
    vector<double> transnotdepeq_val;

    vector<double> nmeq_val;

    vector<long> cope_num;                    // Quantities for the observation probability 
    vector<double> cope_val;
    vector<double> cope_op;

    vector<long> cpe_num;                     // Quantities for the detection probability
    vector<long> cpe_oneminus;
    vector<double> cpe_pd;
    vector<double> cpe_val;
    vector<double> cpe_valoneminus;

    vector<double> capev_val;                 // Qunatities for capture events
    vector<double> capev_valoneminus;                      
    vector<long> capev_num;                    
    vector<long> capev_oneminus;          

    // Used for sampling the initial compartment within a classification
    vector< vector <double> > simcl_initcsumtot;
		vector< vector< vector <double> > > logsimcl_initc, simcl_initc, simcl_initcsum;

    // Used for sampling the enter compartment within a classification
    vector< vector <double> > simcl_entercsumtot;
		vector< vector< vector <double> > > logsimcl_enterc, simcl_enterc, simcl_entercsum;  

    // Used for sampling the initial compartment when simulating
    vector <double> sim_initcsumtot;
		vector< vector <double> > logsim_initc, sim_initc, sim_initcsum;           
    
		// Used for sampling the enter compartment when simulating
    vector <double> sim_entercsumtot;
		vector< vector <double> > logsim_enterc,sim_enterc, sim_entercsum;         
  
		// Used for sampling the initial compartment when simulating for the unobserved individuals
		vector <double> simuo_initcsumtot;
    vector< vector <double> > logsimuo_initc, simuo_initc, simuo_initcsum;           
    
		// Used for sampling the birth compartment when simulating for the unobserved individuals
		vector <double> simuo_entercsumtot;
    vector< vector <double> > logsimuo_enterc, simuo_enterc, simuo_entercsum;         
		
		// Used for addreminf proposal
		vector <long> indinflist;
		vector <vector <vector <long> > > obsinflist; // List of observed individuals
		vector <vector <vector <long> > > inflist;    // List of unobserved infected individuals
		vector <vector <vector <long> > > notinflist; // List of unobserved uninfected individuals
	
		// Used for multimove
		vector <float> jump_multimove;
		vector <float> multimovepr;
		vector <float> ntr_multimove, nac_multimove, nfa_multimove;
		vector <float> ntr_multimovetot, nac_multimovetot;
			
  public:   // Functions used within a chain
  void initchain();
  long start(); 
  void addemptyind(double tbirth); 
  void rememptyind(long i);  
  double priorprob(long pr);  
  double priorsmooth(long s, long j);           
  double priorcalc();                
  void indchange(long i, long estart, double t, double tend);
  void indcha(long i);
  void indrev(long i);
  long eventmove(long i, long e, double tnew);  
  void secchange(long i, long ci, long cf, double ti, double tf); 
	void secchange2(long d, long dn, double *ch_popnum, long valch, double ti, double tf);

  void eqtimelineinit();              
  long insertevent(long d, long e, double t, long eq, long dn, double *dpopnum, long sign);    
  void trydeleteevent(long d, long e);             
  void check(long num);                          
	void checkparamgrad();                        
  void checklike(long num);                       
  void checkderive();                          
  double calcobsprob(long eq);                    
  void indremevent(long i, long e);               
  void indaddevent(long i, long k, EV ev);     
  void chainobsinit();                      
  void param_prop();                             
  void paramnorm_prop();                       
  void paramgrad(long p);
  void param_gradprop();                             
  void param_sampprop();                                       
  void changeparam(long dir, long p, double param_new);      
  void diagnosticschain();                          
	void diagnosticsfile(); 
	void part_prop(long i);                          
  void insertexisting(vector<EV> &ev);           
  void addexistingfute(PART &pa, long tr, double tbeg, double tend);

  void getdep(vector <long> &vec, long co, double t);
  void getnotdep(vector <long> &vec);               
  void chainpartinit();                        
  double transra(long eq, double t, long cnow);
  double transradep(long eq, double t, long cnow);  
  void move_prop(long i);                           
  void addfuture(long i, long c,long cl, double tstart);  
  void addfuturetra(long i, long tr, double tstart, double tadd);
  void sim(double tmin);                            
  void sim_indchange(long i, long ci, long cf); 
  void eventplot();                              
  void sim_check(double R, double t);               
  double likenm(long i);                         
  void part_simsec(long i, long evpbeg, long evpend, PART &pa);   
  void part_sim(long i);       
  double indsimulate(long i);   
  double indsim(long i, long c);  
  double probindsimulate(long i);                                 
  double probindsim(vector<EV> &ev);    
  long selectstart(long i, double t);                           
  double probstart(long i, double t, long c);        
  void indsim_init(long i, double tbirth, long incobs);  
  double part_simprob(long i, long evpbeg, long evpend, PART &pa);
  long getlifespan(long i);    
  void setnmeq();            
  void update();                        
  void traceplot();                                
  void tracefile();
	void nmchange(long tr, double dt, double dt2);    
  void nmupchange(long tr, double dt, double dt2);  
  double L();                                      
  double Lwnm();                                 
  void life_prop(long i);                      
  long lifechange(long i, double tentnew, double tleanre, double tbirthnew, double dprob);  
  void addreminit();                               
  void addrem_prop();                              
  void simsumcalc();                              
  double Lout();                                     
  void plotparamlike(long p, double min, double max);
  long numwrong(long i, vector<EV> &vec);            
  void sing_prop(long i);
  void twothree_prop(long i);                       
  void pair_prop(long i);
  void gap_prop(long i);                          
  void initparamsamp();                             
  void initproc();                                  
  double fit();                                     
  long localresamp(long i, long cl, long estart, double t, double tend, long fixend, long evfac);
  void indsim_prop(long i);                        
  double priorsamp(long pr);
  double fixnmtrans(vector<EV> &ev);       
  void addfixnmfute(long tr, double t, vector<EV> &ev, long e);
  long addevnew(long tr, double t);
  double startseq(long i);
  void nmcha(vector <EV> &ev, long e, long dir);
  void deriveplotinit();                          
  void deriveplotcalc();                            
  void derivepl();                                 
  double probuocap(long i, vector<EV> &vec);          
	
	// CORONA SPECIFIC
	void checkinflist(short num);
	void initinflist();                      
	void addreminfinit();
	void getprob(long r);
	void addreminf();
	void multisec();
	void multisecend();
	void multimoveinit();
	void multimove();
	double movecalcgrad(long tr, double t);
	double transrate(long eq, double t);
	void checklikedisc(long num);
	void addparam();
	void test(short num);
};

const long chainmax = 1;                       // Stores a single MCMC chain
Chain* ch[chainmax];
long nchain;

// used in equation.h
vector<string> getcommasep(string a);          // Takes a string and splits it up by the commas
long getclassval(long cl, string st);          // Returns the value of a classification from a string
vector<long> getallcomp(string st);            // Gets all the comprtments from a string

// used in readinput.h
string get(XMLNode* node, string attr);        // This gets an attribute
double getnum(XMLNode* node, string attr);     // Gets a number fron an XML attribute
void getalltrans(XMLNode* node, XMLNode* node2, long cl, long in, string rate);   // works out all the compartments consistent with a transition
long exist(XMLNode* node, string attr);        // Determines if an attribute exists
long getint(XMLNode* node, string attr);       // Gets an integer XML attribute
vector<string> getlinesep(string a);           // Seperates using the | mark
string repla(string st, string sub1, string sub2); // Replaces one substring with another
long getparam(string na);                      // Gets the parametere number fron the name
long addind(string id, long c);                // Adds an individual

// using in output.h
void outputmodel();                            // Outputs the model (used in debugging)
void outcomp(long c);                          // Outputs the state of a compartment
void traceinit();                              // Initialises the trace plot
void oe(string st, vector<EV> ev);             // Outputs an event sequence 
void tracefileinit();                          // Initialises tracefile output

// udes in prior
void priorinit();                              // Initialises the prior

// used in likelihoodinit.h
void likelihoodinit();                         // Initialises quantities for fast likelihood calculations

// used in ratecalc.h
double ratecalc(long eq, double *popnum, double *param);  // Calculates the rate of an event
double ratecalcdep(long d, double *popnum, double *param); // uses ncompind and popnum from the reduced set for a link
double ratecalcnotdep(long eq, double *paramval); // uses ncompind and popnum from the reduced set for a link

// used in observation.h
void observationinit();                        // Initialises the observation model
double popcalc(long p, double val);            // Calculates the observation probability for a population
double dercalc(short j, double *popnum, double *param);    // Calculates the observation probability for a derived quantity

// param_prop.h
void paraminit();                              // Initialisation for parameter proposals
long notinsideprior(long p, double val);       // Determines if a value is within the prior

// part_prop.h
void partinit();                               // Initialises quantities used when sampling from partivles

// ratecalc.h
double ratecalcdeptakeoff(long d, double *popnum, double *param, long c);    // Removes the individual at c when the quantity is calculated

// sim,h
void siminit();                                // Initialises variables for use in simuations

// dist.h
double gammaprob(double x, double a, double b);// The log of the probability from the gamma distribution
double dgammaprob(double x, double x2, double a, double b);  // The difference in the log of the probability from the gamma distribution

double gammaup(double x, double a, double b);  // The log of the integral of the gamma function from x up to infinity
double gammasamp(double a, double b);          // Draws a sample from x^(a-1)*exp(-b*x)
void gammainit();                              // Initialises the gamma distribution

double weibullprob(double x, double lam, double kk);  // The log of the probability of sampling from the Weibull distribution
double weibullup(double x, double lam, double kk);    // The log of the integral of the Weibull function from x up to infinity
double weibullsamp(double lam, double kk);            // Samples from the Weibull distribution
double dweibullprob(double x, double x2, double lam, double kk);   // The differnce log of the probability of sampling from the Weibull distribution
double dweibullup(double x, double x2, double lam, double kk);     // The difference in the log of the integral of the Weibull function from x up to infinity

// startseq.h
void startseqinit();                           // Initialises quantitites used in setingup the initial sequence
long getctop(double t, double tbirth);         // Returns the compartment part related to the top three classifications
long checkforfixedinpath(long c, long cc, long ctopp);// Considers if fixed events stop a path from occuring

// check.h
void checkevseq(vector<EV> &vec);              // CHecks that an event sequnece is consistent and correct
void checkevseqsimp(vector<EV> &vec);

// multimove.h
void addsettime();                             // Adds settime transitions into an event sequence
void remsettime();                             // Removes settime transitions into an event sequence

// Common function
double ran(){ 
	double v = (double(rand())*29999+rand())/(30000.0*RAND_MAX); if(v == 0 || v == 1) return 0.1;
	else return v;
}

double normal(float mu, float sd){ return mu + sd*sqrt(-2*log(ran()))*cos(2*3.141592654*ran());}
long len(string s){ return s.length();}

void storeemsg(string msg)
{
	stringstream sdi; sdi << root << "/error" << simnum << ".txt";
	ofstream err(sdi.str().c_str());
  err << msg << "    Sample number:" << samp << "\n";
	err.flush();
}

void emsg(string msg)
{
  cout << "e|ERROR MSG: " << msg << "    Sample number:" << samp << "\n";
 
	storeemsg(msg);
	
	#ifdef MP
	MPI_Finalize();
	#endif
	
  exit (EXIT_FAILURE); 
}

void emsg(string msg, double v1){ stringstream ss; ss << msg << ":" << v1; emsg(ss.str());}
void emsg(string msg, double v1, double v2)
{
	stringstream ss; ss << msg << ":" << v1 << "," << v2; emsg(ss.str());
}

void emsg(string msg, double v1, double v2, double v3)
{ 
	stringstream ss; ss << msg << ":" << v1 << "," << v2 << "," << v3; emsg(ss.str());
}

void emsg2(string msg)
{
  warning.push_back(msg);
}
