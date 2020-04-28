
# CoronaBICI

C. M. Pooley* [1], Andrea Doeschl-Wilson [2], and Glenn Marion [1]

[1] Biomathematics and Statistics Scotland, James Clerk Maxwell Building, The King's Buildings, Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK 

[2] The Roslin Institute, The University of Edinburgh, Midlothian, EH25 9RG, UK. 

* Corresponding author
Email: chris.pooley@bioss.ac.uk


BICI, which stands for "Bayesian Individual-based Compartmental Inference", is a general purpose software for analysing individual and population-wide data using compartmental models. CoronaBICI is a version of BICI specifically designed to analyse Corona virus data (with bespoke efficient MCMC proposals suitable for large population sizes).

Analysis proceeds in two stages:

## 1) GENERATING AN INPUT FILE FOR BICI

This is performed in the "GenerateInputFile" directory. 

a) The file "generateinputfile.cc" is compiled and run using the command:

./gen Scotland model

This programatically generates the file "Scotland_model.txt", which is a regional model of Scotland (see Appendix A in "BICI Manual_v1.0.pdf" for a description of how this file specifies the model). Note, using the argument "UK" instead of "Scotland" generates an equivalent file for the UK. Census information (taken from government websites for England, Scotland, Wales and Northern Ireland) was used to get the populations in the various regions.

b) The file "Scotland_model.txt" is then read into the BICI GUI which outputs the same model but using XML format ("Scotland_bici_model.xml"). The BICI GUI itself is not currently included in Github because it was found to be too large to upload.

c) The file "Scotland_bici_model.xml" then gets read by "generateinputfile.cc" a second time using the command:

./gen Scotland input

This time data is added to finally create the file "Scotland_bici_input.xml", which sumarises the model and data in a way that BICI can understand. The data itself consists of daily regional cases from the website https://smazeri.shinyapps.io/Covid19_Scotland/.   


## 2) RUNNING THE ANALYSIS

This is done in the "Analysis" directory. The "bici.cc" file is compiled using

mpic++ bici.cc header/tinyxml2.cc -O3 -o bici

and run using:        

mpirun -n 4 ./bici Scotland_bici_input.xml 10000

Here -n 4 represents the number of parallel MCMC chains and 10000 represents the number of samples on each chain. 
BICI creates the "Scotland_bici" subdirectory in "Outputs" and  writes to a number of different file during program execution:

a) "Trace[X].txt" (where X is the number of the chain) gives trace plots of the model parameters. These can be loaded into other software to visually check that MCMC is mixing well (e.g. Tracer).

b) "Diagnostic[X].txt" gives information about the MCMC proposals to help identify inefficiencies in the code

c) "Bici[X].txt" is an output file which can be read into the BICI GUI for visualisation.

d) "Stats.txt" calculates the posterior means, 90% credible intervals, effective sample sizes and Gelman-Rubin diagnostic statics (these last two are used to confirm that MCMC is well mixed).

## METHODOLOGY

This is yet to be written up for BICI, but the basic framework is presented here (for a simple SIR model):

Estimating individuals' genetic and non-genetic effects underlying infectious disease transmission from temporal epidemic data,
C. M. Pooley, G. Marion, S. C. Bishop, R. I. Bailety, and A. B. Doeschl-Wilson
bioRxiv 618363 (2019).

