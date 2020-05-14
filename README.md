
# CoronaBICI

![](https://github.com/ScottishCovidResponse/CoronaBICI/workflows/CI/badge.svg?branch=master) (See the [build and test results](https://github.com/ScottishCovidResponse/CoronaBICI/actions?query=workflow%3ACI).)

C. M. Pooley† [1], Andrea Doeschl-Wilson [2], and Glenn Marion [1]

[1] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK 

[2] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

† Corresponding author

Email: [chris.pooley@bioss.ac.uk](mailto:chris.pooley@bioss.ac.uk)

BICI, which stands for "Bayesian Individual-based Compartmental Inference", is a general purpose software for analysing individual and population-wide data using compartmental models. CoronaBICI is a version of BICI specifically designed to analyse Corona virus data (with bespoke efficient MCMC proposals suitable for large population sizes).

For contributing, please see the [Developer Guide](DeveloperGuide.md).

Analysis proceeds in two stages, summarised here:

```
Data/censustable.txt
    |
    |  # Generate a BICI model
    |  ./gen Scotland model
    ↓
Scotland_model.txt
    |
    |  # Convert to XML
    |  (BICI GUI)
    ↓
Scotland_bici_model.xml      Data/covid-19-cases-uk.txt
    |                                 |
    |_________________________________|
    |  # Add COVID case data
    |  ./gen Scotland input
    ↓
Scotland_bici_input.xml
```

In more detail:

### Generating an input file for BICI

This is performed in the "GenerateInputFile" directory. 
* The file "generateinputfile.cc" is compiled and run using the command:
  ```
  ./gen Scotland model
  ```
This programatically generates the file "Scotland_model.txt", which is a regional model of Scotland (see Appendix A in "BICI Manual_v1.0.pdf" for a description of how this file specifies the model). Note, using the argument "UK" instead of "Scotland" generates an equivalent file for the UK. Census information (taken from government websites for England, Scotland, Wales and Northern Ireland) was used to get the populations in the various regions.

* The file "Scotland_model.txt" is then read into the BICI GUI which outputs the same model but using XML format ("Scotland_bici_model.xml"). The BICI GUI itself is not currently included in Github because it was found to be too large to upload.

* The file "Scotland_bici_model.xml" then gets read by "generateinputfile.cc" a second time using the command:
./gen Scotland input

This time data is added to finally create the file "Scotland_bici_input.xml", which sumarises the model and data in a way that BICI can understand. The data itself consists of daily regional cases from the website https://smazeri.shinyapps.io/Covid19_Scotland/.   


### Running the analysis

This is done in the "Analysis" directory. There are two ways to run: 

1. If running as a single process the "bici.cc" file is compiled using
   ```
   make bici
   ```
and run using:
   ```
   ./bici Scotland_bici_input.xml 10000
   ```
   Here 10000 represents the number of samples on each chain. 

2. To run multiple chains under MPI the line #define MP 1 is uncommented in bici.cc and code is compiled using:
   ```
   make bici_mpi
   ```
   and run using:        
   ```
   mpirun -n 4 ./bici_mpi Scotland_bici_input.xml 10000
   ```
Here -n 4 represents the number of parallel MCMC chains and 10000 represents the number of samples on each chain. 

### Output

BICI creates the "Scotland_bici" subdirectory in "Outputs" and  writes to a number of different file during program execution:

1. "Trace[X].txt" (where X is the number of the chain) gives trace plots of the model parameters. These can be loaded into other software to visually check that MCMC is mixing well (e.g. Tracer).

2. "Diagnostic[X].txt" gives information about the MCMC proposals to help identify inefficiencies in the code

3. "Bici[X].txt" is an output file which can be read into the BICI GUI for visualisation.

4. "Stats.txt" calculates the posterior means, 90% credible intervals, effective sample sizes and Gelman-Rubin diagnostic statics (these last two are used to confirm that MCMC is well mixed).

### Methodology

This is yet to be written up for BICI, but the basic framework is presented here (for a simple SIR model):

Estimating individuals' genetic and non-genetic effects underlying infectious disease transmission from temporal epidemic data, C. M. Pooley, G. Marion, S. C. Bishop, R. I. Bailety, and A. B. Doeschl-Wilson, [bioRxiv 618363](https://www.biorxiv.org/content/10.1101/618363v3.full) (2019).

