# diff_tagg_ana 

## Instructions to get started

Contact: Kong Tu (BNL) and Bill Li (SBU)

Disclaimer: this only applies at BNL RCF. If one wants to try virtualbox, please go to this link,https://github.com/ECCE-EIC/Singularity/blob/master/VirtualBox.md.

- Go to your working directory; in the following example, I do it in my home directory for simplicity.
"/direct/eic+u/zhoudunming"

- copy two files from me, 
  - ``cp /direct/eic+u/zhoudunming/singularity_shell.sh ./``
  - ``cp /direct/eic+u/zhoudunming/setup.sh ./``

- source both by, 
  - ``source singularity_shell.sh``
  - ``source setup.sh``

- Create a directory for this work, I called it "ECCE"
  - ``mkdir ECCE``
  - ``cd ECCE``

- Get the repository
  - ``git clone https://github.com/KongTu/diff_tagg_ana.git``

- Build the project,
  - ``cd diff_tagg_ana``
  - ``chmod +x autogen.sh``
  - ``mkdir build``
  - ``cd build``
  - ``../autogen.sh --prefix=$MYINSTALL``
  - ``make install``

- Running DST Tree reader
  - ``cd ../macros``
  - ``bash copy.sh`` (this is to copy a test file for you to continue; it will store the DST tree under "diff_tagg_ana/DSTrees"
  - Now the place for listing the input DST Trees is in "myFileList.txt". The default has the "test.root" that you just created. Therefore, you are good for running it directly
  - ``root -l Fun4all_reana.C``, now the number of events = 100 in this macro; set to -1 for all events.
  - The output of this will go to "../flatTree" with the name "phi.root" (because that's the test file)

- Running Flat Tree reader (last step!)
  - ``root -l readFlatTree.C``
  - Check output under "results" folder.
  - you are done.

Happy analyzing!
  
