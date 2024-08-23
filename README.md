# GOrGONA: A Global Optimization (GO) code for NanoAlloys (NA) exploiting Grouping (Gr) Algorhthms
This version of the code exploits pre-existing codes by integrating them: BHGO (Basin Hopping Global Optimization) for Nano Alloys developed by the group of Prof. Ferrando (see Ref. [1]); AugGGO (Augmented Grouping GO ) developed by Dr. Barcaro and Dr. Fortunelli (see Refs. [2] and [3]).
This file illustrates the basic usage of the code. Please follow the steps as described. \
For questions or problems, please contact us at the following mail: giovanni.barcaro@cnr.it

Preliminary step: download all the files and the folders of the archive.

Step 0: Go into the folder EAM-INTERNAL and type "make" to compile the program ("Makefile" is editable): this step will create the "bh_v07" executable which will be used as the Basin Hopping tool exploiting Grouping algorithm for the optimization of the chemical ordering. The code has been tested by using "gfortran" compiler. We suggest to use such compiler, as suggested in the following description. \
Once created, copy the "bh_v07" executable to the folder "parallel_unseeded_XxYy" and "parallel_seeded_XxYy".

Step 1: rename the folders "parallel_unseeded_XxYy" and "parallel_seeded_XxYy" by choosing the couple of metals that you want to investigate. For example, Xx=Ag and Yy=Cu.

Step 2: open the bash script "01_search_unbiased". This script will drive the first generation of structures for a chosen size of the cluster and on the desired range of compositions. \
Most of the parameters can be changed, but only by expert users (please, contact us for more detailed information in this sense). For beginners, we suggest to leave these parameters unchanged and follow these steps:\
2.1: (line 9) change "parallel_unseeded_XxYy" by specifying the two metals you have chosen;\
2.2: (line 16) change "SizeEnd1" by choosing the size of the cluster;\
2.3: (line 20) change "SizeStart1" in 0, specify the same "SizeEnd1" chosen at the previous step and choose the mesh to sweep composition range by specifying "Pruning1": for example, if "SizeEnd1"=100 and if you choose "Pruning1"=20, you will run parallel BH simulations at the following compositions: (20,80), (40,60), (60,40) and (80,20) for mixed clusters and two more runs for the pure compositions (100,0) and (0,100).

Step 3: go into the folder "parallel_unseeded_XxYy" and:\
3.1: put in the folder the correct force field file; a collection of force fields is provided in the "database_ff" folder;\
3.2: open the file "analize", where you have to:\
3.2.1: specify the correct couple of metals by replacing "Xx" and "Yy"; \
3.2.2: specify the correct value of "SizeEnd1".\
3.3: open the files "clean" and "crea", where you have to specify the correct value of "SizeEnd1".\
3.4: open the files "input_bh..." and choose the desired number of BH steps ("nbh1", line 8) and specify the correct name of the force field file (line 6). For beginners, we suggest to use "nbh" values between 10000 (for clusters composed by about 200 atoms) and 30000/50000 for clusters composed by 600/800 atoms.\
3.5: open the program "parallel_unseeded.f90" and:\
3.5.1: specify the correct couple of metals by replacing "Xx" and "Yy";\
3.5.2: pay attention to tuning the dimension of the box (lines 17-19) where the random starting structure is created. The following values are suggested: for a cluster of 200 atoms choose (-6 to +6); for a cluster of 800 atoms choose (-10 to +10);\
3.5.3: compile the program using the following command: "gfortran -o parallel_unseeded.x parallel_unseeded.f90";\
3.5.4: compile the program "read.f90" using the following command: "gfortran -o read.x read.f90".

Step 4: execute the bash script "01_search_unbiased" in batch mode with the command: "nohup ./01_search_unbiased < /dev/null >& LOG_U &". This command will send the BH simulations (simulations will run in parallel at different compositions by spanning the chosen range). Time requested strongly depends on the size of the clusters ("SizeEnd1") and on the required number of steps ("nbh1").

Step 5: open the bash script "021_identify_pures" and: \
5.1: specify the correct "SizeStart1" (lines 3, 7, 16 and 19) by using (I4) format: for example, if "SizeStart1"=0, you have to specify "0000";\
5.2: consistently, specify the correct "SizeEnd1" (lines 3, 12, 17 and 19) by using (I4) format: for example, if "SizeEne1"=100, you have to specify "0100";\
5.3: specify the correct couple of metals by replacing "Xx" and "Yy".

Step 6: run the bash script "021_identify_pures" using the following command: "./021_identify_pures". The file "res-pures.dat" will be created.

Step 7: compile the program "022_take.f90" using the following command: "gfortran -o 022_take.x 022_take.f90" and execute it using the following command: "./022_take.x". The script "./023_copy.sh" will be created.

Step 8: run the bash script "023_copy.sh" using the following command: "./023_copy.sh". In the main folder, two new folders will appear: "run_'SizeStart1'" and "run_'SizeEnd1'". In the case of "SizeStart1"=0 and "SizeEnd1"=100, they will be "run_0000" and "run_0100". The program will use the information in these folders for a common energy reference of pure clusters.

Step 9: open the bash script "024_search_unbiased.REVISE" and: \
9.1: specify the correct "SizeStart1" (lines 7 and 9) and "SizeEnd1" (lines 8 and 10) by using (I4) format.

Step 10: run the script "024_search_unbiased.REVISE" using the following command: "./024_search_unbiased.REVISE"; "curve-..." files will be updated.

Step 11: compile the program "02_best.f90" using the following command: "gfortran -o 02_best.x 02_best.f90" and execute it using the following command: "./02_best.x". The file "result.out" will be created. It corresponds to the convex hull obtained by comparing the results in the "curve-..." files.

Step 12: compile the program "03_take.f90" using the following command: "gfortran -o 03_take.x 03_take.f90" and execute it using the following command: "./03_take.x". The "DATA" folder will appear and the script "04_extract.sh" will be generated. By executing it ("./04_extract.sh") the DATA folder will be populated with xyz files corresponding to the lowest-energy structure at each investigated composition found in the BH runs.

First part of the search (unbiased one) is complete and we can proceed (if desired) to a refined search using the collected motifs in the DATA folder. This second part is defined 'biased' as it exploits a database of known structures. To proceed follow these steps:

Step 13: open the bash script "05_search_biased" and follow these steps:\
13.1: change "parallel_seeded_XxYy" by specifying the two metals you have chosen (lines 12, 40 and 70);\
13.2: change "SizeEnd2", "SizeEnd2" and "Pruning2": the first two specify the range of compositions to investigate; it can be the same range investigated in the unbiased search, but it can be also a different range; same considerations hold for "Pruning2". If "SizeStart1"=0 and "SizeEnd1"=100 and "Pruning1"=20, if we want to explore the same range, we can set "SizeStart2"=20, "SizeEnd2"=80 and "Pruning2"=20;\
13.3: change also "SizeStart1" and "SizeEnd1" by using (I4) format (see Step 5.1).

Step 14: go into the folder "parallel_seeded_XxYy" and:\
14.1: put in the folder the correct force field file; \
14.2: open the file "analize", where you have to:\
14.2.1: specify the correct couple of metals by replacing "Xx" and "Yy"; \
14.2.2: specify the correct value of "SizeEnd2".\
14.3: open the files "clean" and "crea", where you have to specify the correct value of "SizeEnd2".\
14.4: open the file "input_bh.in" and choose the desired number of BH steps ("nbh2", line 8) and specify the correct name of the force field file (line 6). For beginners, we suggest to use "nbh2" values between 5000 (for clusters composed by about 200 atoms) and 15000/20000 for clusters composed by 600/800 atoms.\
14.5: open the program "parallel_seeded.f90" and:\
14.5.1: specify the correct couple of metals by replacing "Xx" and "Yy";\
14.5.2: compile the program using the following command: "gfortran -o parallel_seeded.x parallel_seeded.f90";
14.6: compile the program "read.f90" using the following command: "gfortran -o read.x read.f90".

Step 16: execute the bash script "05_search_biased" in batch mode with the command: "nohup ./05_search_biased < /dev/null >& LOG_B &". This command will send (in parallel) the seeded BH simulations optimizing chemical order only. 

Step 17: compile the program "06_best.f90" using the following command: "gfortran -o 06_best.x 06_best.f90" and execute it using the following command: "./06_best.x"; a new "result.out" will be created by putting together all the results coming from the "curve-.." files.

Step 18: compile the program "07_take.f90" using the following command: "gfortran -o 07_take.x 07_take.f90" and execute it using the following command: "./07_take.x"; a new folder, "DBFN", will appear and the script "08_extract.sh" will be generated. By executing it ("./08_extract.sh") the DBFN folder will be populated with xyz files corresponding to the lowest-energy structure at each investigated composition found in the refined BH runs.

References\
[1] D. Rapetti, C. Roncaglia, R. Ferrando, "Optimizing the Shape and Chemical Ordering of Nanoalloys with Specialized Walkers", ADVANCED THEORY AND SIMULATIONS, 2300268 (2023)\
[2] G. Barcaro, L. Sementa, A. Fortunelli, "A grouping approach to homotop global optimization in alloy nanoparticles", PHYSICAL CHEMISTRY CHEMICAL PHYSICS, 16, 24256-24265 (2014)\
[3] D, Fioravanti, G. Barcaro, A. Fortunelli, "An augmented (multi-descriptor) grouping algorithm to optimize chemical ordering in nanoalloys", PHYSICAL CHEMISTRY CHEMICAL PHYSICS, 23, 23075-23089 (2021)

Acknowledgments\
We acknowledge financial support within Task WP5.2.1 of Work Package 5 on "Materials Foundry" of Spoke 7 of the ICSC, the Italian "Centro Nazionale di Ricerca in High Performance Computing, Big Data and Quantum Computing", funded by the the European Union (Next Generation EU, grant number CN00000013) - PNRR, Mission 4 Component 2 Investment 1.4.
