# fully-dynamic-k-center
This is the source code for an Algorithms submission.    
This repository contains the source code of a program related to an article submitted to PKDD22. 

Three algorithms are implemented in the package, including fully dynamic algorithm(FD), offline greedy algorithm(OG), and offline linear program algorithm(LP). All algorithms approximate the robust k-center problem.

-- Compilation --

To compile it, use 'make' in the root directory.
It was only tested on ubuntu/debian/fedora. 

-- Command Line arguments -- 

The different algorithms implemented are:
    - FD, the corresponding option is -m.
    - OG, the corresponding option is -p.
    - LP, the corresponding option is -s.

One must specify one of those mandatory options to run the program.

The program takes one of the mendatory options and 7 others mandatory arguments:
    - for FD, it requires the arguments in this order:
      	  ./k-center -m k t d_min d_max data_file query_file tau
    - for both OG and LP, it requires arguments in this order:
      	  ./k-center -p/-s k tau d_min d_max data_file query_file outliers

Because FD can maintain the layers dynamically, people may change the number of outliers from time to time. The value of outliers for FD can be set in function fully_adv_write_log() in algo_fully_adv.c when you need to output the solution.

The option -l file_name create a log file with one line per operation having the following format:
level_of_solution radius_of_level true_radius_of_level

-- Data file format -- 

For all options -m, -p and -s:
The data_file should have one entry per line and every line should have the following format :
timestamp tab longitude space latitude

The timestamp should be integer that fits in an unsigned int.
The index of point is fixed by their order in the file, ranging from 0 to n-1 where n is the numbe of lines of the data file.

As for now, the euclidean distance is used to compute the distance between 2 points.

-- Queryfile format -- 

The queryfile should be a binary file containing a sequence of unsigned int. The endianness and size of elements in the query file should correspond to the endianness and size of the computer it is run on. 

With all options of -m, -p and -s, the parity of an occurence determines if the corresponding point should be added or removed. If odd, the point is inserted, if even, it is removed.

-- Example --

In the 'example ' file, a small dataset (example_dataset.txt) which includes 4000 points and a query (example_query) which inserts all of the points are provided for test. The distribution of the points can be seen in example_dataset.png.

To test the code, you may first compile the code. Please refer to  'Compilation ' section to learn how to compile. Then you can follow 'Command Line arguments' section to learn how to use the commands.

For example, if you input:

./k-center -l 'log' -p 10 0.1 1 1000 'example/example_dataset.txt' 'example/example_query' 100

You will run the OG algorithm with k=10, dmin=1, dmax=1000, tau=0.1 (distance precision used for binary search) and outliers=100. The data file and query file are 'example/example_dataset.txt' and 'example/example_query' respectively. And you may see a log file named 'log' is generated.
It will show 'true_radius: 347.852864' if you make every steps correctly.

-- Other -- 

Below 'lib' path, there are some Gurobi libraries, which are used to solve linear program sub-routine. A Gurobi license must be exported first for the running of LP algorithm. Please refer to https://www.gurobi.com/ to learn more about the license.
