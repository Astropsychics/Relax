
Information on the Relax Monte Carlo transport simulation software

Note: All the table data found in the Tables directory is compressed via gzip before
being placed in the github repo. This is due to large file sizes for some of the scattering
data tables used. They must be uncompressed locally on your machine by

> ls
Data Execute Inputs Modules README.txt Routines Tables
> cd Tables
> gzip -d -k *

which keeps the original zipped files as well as creating unziped .dat files.
Keeping the original zipped files ensures that the git repo is happy that it is
not missing files which are included within the repo (*.gz files). 

Similarly, before being pushed to the repository, if any changes are commit or new tables
added they must be compressed first

> cd Tables
> gzip *



