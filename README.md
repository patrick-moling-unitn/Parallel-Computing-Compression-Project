# Parallel Computing project - Compression

### Goal of this project
Design and develop a parallel compression tool that efficiently compresses large datasets using multiple processing units. File compression algorithms (e.g., Huffman, LZ77, LZMA) are computationally intensive and can benefit significantly from concurrent processing by dividing the input data into chunks that can be compressed independently. 

### How to reproduce the results
Download the files from git and copy the folder into the cluster
> scp -r ./YourGitFolderName/ yourname.yoursurname@hpc.unitn.it:~/YourDestinationFolder

Inside your cluster, in YourDestinationFolder, compile the **main.cpp**
> g++ -std=c++11 -fopenmp -o main main.cpp mmio.c

In case of errors due to Windows' *^M* execute the following commands and repeat the compilation attempt
> dos2unix schedulable_job.pbs
> 
> dos2unix main.cpp

Run the job you'd like to
> qsub schedulable_job.pbs

Check if the job has finished
> qstat -u yourname.yoursurname

Verify the output once done
> cat execution.out