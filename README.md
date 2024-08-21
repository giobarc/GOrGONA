# GOrGONA
A Global Optimization (GO) code for NanoAlloys (NA) exploiting Grouping (Gr) Algorhytms

This file illustrates the basic usage of the code. Please follow the steps as described. \
For questions or problems, please contact me at the following mail: giovanni.barcaro@cnr.it

Step 0: Go into the folder EAM-INTERNAL and type 'make' to compile the program (Makefile is editable): this step will create the 'bh_v07' executable which will be used as a Basin Hopping tool exploiting Grouping algorithm for the optimization of the chemical ordering. \
Once created, copy the 'bh_v07' executable to the folder "parallel_unseeded_XxYy" and "parallel_seeded_XxYy".

Step 1: rename the folders "parallel_unseeded_XxYy" and "parallel_seeded_XxYy" by choosing the couple of metals that you want to investigate. For example, Xx=Ag and Yy=Cu.

Step 2: open the bash script 01_search_unbiased. This script will drive the first generation of structures for a chosen size of the cluster and on the desired range of compositions. \
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

3.4: open the files "input_bh..." and choose the desired number of BH steps ("nbh1", line 8) and specify the correct name of the force field file (line 6). For beginners, we suggest to use nbh values between 10000 (for clusters composed by about 200 atoms) and 30000/50000 for clusters composed by 600/800 atoms.

3.5: open the program "parallel_unseeded.f90" and:\
3.5.1: specify the correct couple of metals by replacing Xx and Yy;\
3.5.2: pay attention to tuning the dimension of the box (lines 17-19) where the random starting structure is created. The following values are suggested: for a cluster of 200 atoms choose (-6 to +6); for a cluster of 800 atoms choose (-10 to +10);\
3.5.3: compile the program using the following command: "gfortran -o parallel_unseeded.x parallel_unseeded.f90"

Step 4: compile the program read.f90 using the following command: "gfortran -o read.x read.f90".

Step 5: execute the script 01_search_unbiased in batch mode with the command: " nohup ./01_search_unbiased < /dev/null >& LOG_U &. This command will send (in parallel) the BH simulations. Time requested strongly depends on the size of the clusters and on the required number of steps.

Step 6: open the script "021_identify_pures" and: \
6.1: specify the correct 'SizeStart1' (lines 3, 7, 16 and 19) by using (I4) format: for example, if 'SizeStart1'=0, you have to specify 0000;\
6.2: consistently, specify the correct 'SizeEnd1' (lines 3, 12, 17 and 19) by using (I4) format: for example, if 'SizeEne1'=100, you have to specify 0100;\
6.3: specify the correct couple of metals by replacing Xx and Yy.

Step 7: run the script "021_identify_pures" using the following command: "./021_identify_pures". The file "res-pures.dat" will be created.

Step 8: compile the program 022_take.f90 using the following command: "gfortran -o 022_take.x 022_take.f90" and execute it using the following command: "./022_take.x". The script "./023_copy.sh" will be created.

Step 9: run the script "023_copy.sh" using the following command: "./023_copy.sh". In the main folder, two new folders will appear: "run_'SizeStart1'" and "run_'SizeEnd1'". In the case of SizeStart1=0 and SizeEnd1=100, they will be run_0000 and run_0100. The program will use the information in these folders for an energy reference of pure clusters.

Step 10: open the script 024_search_unbiased.REVISE and: \
10.1: specify the correct 'SizeStart1' (lines 7 and 9) and 'SizeEnd1' (lines 8 and 10) by using (I4) format.

Step 11: run the script "024_search_unbiased.REVISE" using the following command: "./024_search_unbiased.REVISE"; "curve-..." files will be updated.

Step 12: compile the program 02_best.f90 using the following command: "gfortran -o 02_best.x 02_best.f90" and execute it using the following command: "./02_best.x". The file "result.out" will be created. It corresponds to the convex hull obtained by comparing the results in the "curve-..." files.

Step 13: compile the program 03_take.f90 using the following command: "gfortran -o 03_take.x 03_take.f90" and execute it using the following command: "./03_take.x". The "DATA" folder will appear and the script "04_extract.sh" will be generated. By executing it ("./04_extract.sh") the DATA folder will be populated with xyz files corresponding to the lowest-energy structure at each investigated composition found in the BH runs.

First part of the search (unbiased one) is complete and we can proceed (if desired) to a refined search using the collected motifs in the DATA folder. This second part is defined "biased" as it exploits a database of known structures. To proceed follow these steps:

Step 14: open the bash script 05_search_biased and follow these steps:

14.1: change 'parallel_seeded_XxYy' by specifying the two metals you have chosen (lines 12, 40 and 70);\
14.2: change "SizeEnd2", "SizeEnd2" and "Pruning2": the first two specify the range of compositions to investigate; it can be the same range investigated in the unbiased search, but it can be also a different range; same considerations hold for "Pruning2". If "SizeStart1"=0 and "SizeEnd1"=100 and "Pruning1"=20, if we want to explore the same range, we can set "SizeStart2"=20, "SizeEnd2"=80 and "Pruning2"=20;\
14.3: change also "SizeStart1" and "SizeEnd1" by using (I4) format.

Step 15: go into the folder parallel_seeded_XxYy and:

Step 14: open the bash script 05_search_biased
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
