
# How to run BICI on CSD3

- Get an account on CSD3 (see https://github.com/ScottishCovidResponse/SCRCIssueTracking/wiki/DiRAC-Access)
- Create a directory for your user in the project space ~/rds/rds-dirac-dc003
  ```
  mkdir ~/rds/rds-dirac-dc003/$USER
  ```
- Copy the source code directory to CSD3 (not in the rds directory
  above), for example, from your local machine,
  ```
  rsync -avz BICI/ login-cpu.hpc.cam.ac.uk:BICI/
  ```
- Change into the source code directory on CSD3 and compile as usual:
  ```
  cd BICI/Analysis
  mpicxx bici.cc header/tinyxml2.cc -O3 -o bici
  ```
- Run the "submit-bici" script to submit your job:
  ```
  ./submit-bici --input Scotland_bici_input.xml --samples 1000 --nprocs 4 --walltime 1:00:00 --dir ~/rds/rds-dirac-dc003/$USER/bici_001
  ```
