The purpose of these scripts is to expedite the retreival of CRAB3 output from the T2 at which they are stored.  The current version of these scripts is only capable of retreiving files from Vanderbilt.  Future versions of this code will extend to other T2s.  The workflow is as such:

1. Use an lcg-ls command to generate a master list of all files to be transferred.  The file must be called MASTER.txt.  IF logs are transferred or jobs have failed, make sure the "log" and "failed" directories are not within this list!  Please open the MASTER.txt file in your favorite editor and remove these lines manually.  Future versions of this code will have safety catches for this.
   Example: $ lcg-ls -b -D srmv2 'srm://se1.accre.vanderbilt.edu:6288/srm/v2/server?SFN=/lio/lfs/cms/store/user/<path-to-job-results>' >> MASTER.txt

2. Run splitPD_Vandy.py to divide the master list into five sub-lists named files<1-5>.txt.  This script will not run if MASTER.txt does not exist.  
   $ python splitPD_Vandy.py

3. In up to five separate shells, run the transfer*.py scripts to copy the files from the T2 to the directory in which these scripts are run.  These scripts will not run if files<1-5>.txt do not exist.
   $ python transfer1.py

4. When all files are transferred one can utilize the CRABgenocide.py script to remove the files from the T2.  This will prevent future angry e-mails complaining about disk space issues.  WARNING: This script uses lcg-del to remove the files from the remote storage element.  Some T2s do not want users manually deleting files from their storage elements (MIT for example).  Please check the rules for deleting files at the SE you are usuing before running this script!  This script will not run if MASTER.txt does not exist.
   $ python CRABgenocide.py
