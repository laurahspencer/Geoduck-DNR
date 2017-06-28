To "Hack" a blib file:

Go to Pecan4 output director.
Under "Percolatorr" director there's a series of job files.
Take the first geoduck job file, and say: `chmod +x [jobfilename]` execute
  This tells Linux that it's an executable file, and run that.
  Then execute `./[jobfilename]`  
  Then back out of the percolator director, and go to the pecan2blib directory. 
  The do the same `chmod +x [jobfilename]` in the pecan2blib directory
  
  Overall, run the: 
  /geoduck.job file to aggregate all the separate isolation window results, then...
  /percolator.job file to aggregate all the samples, then..
  /pecan2blig.job file to condense everything into one .blib file. 
  
ALSO: in all the percolator sample.job and joint.job files have the -Q input. If the error log shows that it doesn't know what -Q is, remove it from all the job 
  
  