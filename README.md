# GOrGONA
A Global Optimization (GO) code for NanoAlloys (NA) exploiting Grouping (Gr) Algorhytms

This file illustrates the basic usage of the code. Please follow the steps as described. \
For questions or problems, please contact me at the following mail: giovanni.barcaro@cnr.it

Step 0: Go into the folder EAM-INTERNAL and type 'make' to compile the program (Makefile is editable): this step will create the 'bh_v07' executable which will be used as a Basin Hopping tool exploiting Grouping algorithm for the optimization of the chemical ordering. \
Once created, copy the 'bh_v07' executable to the folder "parallel_unseeded_XxYy" and "parallel_seeded_XxYy".

Step 1: rename the folders "parallel_unseeded_XxYy" and "parallel_seeded_XxYy" by choosing the couple of metals that you want to investigate. For example, Xx=Ag and Yy=Cu.

Step 2: open the bahs script 01_search_unbiased. This script will drive the first generation of structures for a chosen size of the cluster and on the desired range of compositions. \
Most of the parameters can be changed, but only by expert users (see Section "Expert users"). For beginners, we suggest to leave these parameters unchanged and follow these steps:

2.1: (line 9) change 'parallel_unseeded_XxYy' by specifying the two metals you have chosen;\
2.2: (line 16) change "SizeEnd1" by choosing the size of the cluster;\
2.3: (line 20) change "SizeStart1" in 0, specify the same "SizeEnd1" and choose the mesh to sweep composition range by specifying "Pruning1": for example, let's suppose that "SizeEnd1"=100; if you choose "Pruning1"=20, you will run parallel BH simulations at the following compositions: (20,80), (40,60), (60,40) and (80,20) for mixed clusters and two more runs for the pure compositions (100,0) and (0,100).

Step 3: go into the folder parallel_unseeded_XxYy and:

3.1: put in the folder the correct force field file; a collection of force fields is provided in the database_ff folder;\
3.2: open the file "analyze", where you have to:\
3.2.1: specify the correct couple of metals by replacing Xx and Yy; \
3.2.2: specify the correct value of "SizeEnd1".

3.3: open the files "clean" and "crea", where you have to specify the correct value of "SizeEnd1".

3.4: open the files "input_bh..." and choose the desired number of BH steps ("nbh1"=100, line 8) and specify the correct name of the force field file (line 6).

3.5: open the program "parallel_unseeded.f90" and:\
3.5.1: specify the correct couple of metals by replacing Xx and Yy;\
3.5.2: pay attention to tuning the dimension of the box (lines 17-19) where the random starting structure is created. The following values are suggested: for a cluster of 200 atoms choose (-6 to +6); for a cluster of 800 atoms choose (-10 to +10);\
3.5.3: compile the program using the following command: "gfortran -o parallel_unseeded.x parallel_unseeded.f90"

Step 4: compile the program read.f90 using the following command: "gfortran -o read.x read.f90".

Step 5: execute the script 01_search_unbiased in batch mode with the command: " nohup ./01_search_unbiased < /dev/null >& LOG_U &.\

 Step04: script 021_identify_pures 
 - modify the name of the metals and the dimension of the cluster
 - run the script
 - check the file res-pures.dat \
 Step05: execute ./022_take.x (nothing to change here) and check the script 023_copy.sh \
 Step06: execute ./023_copy.sh \
 Step07: open the script 024_search_unbiased.REVISE 
 - change the size of the system and the range
 - run the script \
 Step08: execute ./02_best.x (nothing to change here) and check the file result.out \
 Step09: execute ./03_take.x (nothing to change here) and check the script 04_extract.sh \
 Step10: execute the script ./04_extract.sh \
 Step11: open the script 05_search_biased
 - change size of the cluster
 - change range as a function of particle size
 - be careful to the name of the folder "parallel_seeded_XXXX": is it the correct mix of metals?
 - pay attention to the number of cycles that the control loop has to perform: I suggest 100000 \
 Step12: RENAME and go into the folder parallel_seeded_XXXX
 - put in the folder the correct potential file; (be careful that the format is consistent, in particular three lines for the cut-off are necessary)
 - in the file "analize" specify the correct couple of metals and the range: the count MUST arrive to full size!
 - in the files "clean" and "crea" choose the correct range
 - in the files "input_bh..." choose the correct number of steps and specify the correct file for the potential
 - in the file "parallel_seeded.f90" specify the correct couple of metals; then recompile the program
 - the files "range.dat" and "param.dat" are automatically tuned by the script \
 Step13: send the script 05_search_biased \
(Step14: open the script 051_search_biased.REVISE) \
(- change the size of the system and the range) \
(- run the script) \
 Step15: execute ./06_best.x (nothing to change here) and check the file result.out \
 Step16: execute ./07_take.x (nothing to change here) and check the script 08_extract.sh \
 Step17: execute the script ./08_extract.sh
