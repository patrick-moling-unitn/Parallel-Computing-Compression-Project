#include <iostream>
#include <fstream>
#include <string>

//Citation: https://www.researchgate.net/publication/371729160_Parallel_Data_Compression_Techniques

// Huffman algorithm available on "https://github.com/bhumijgupta/huffman-compression-library"
#include "huffmantool.h"

#include <algorithm>  //std::sort

//For OpenMP
//>
#ifdef _OPENMP
#include <omp.h>
#endif
//<

//For MPI
//>
#include <stdio.h>
#include <mpi.h>
//<

using namespace std;

#define EXECUTION_TIMES 10

void AnalizeCompressionRatio(string original, string compressed, string decompressed){
	float originalSize = sizeof(original) * original.size(), compressedSize = sizeof(compressed) * compressed.size();
	float compressionRatio = (1 - compressedSize / originalSize) * 100;
	cout << "--------------------------------" << endl;
	string dataStatus = original == decompressed ? "working as expected! [Success]" : "corrupting the original data! [Error]";
	cout << "Data decompression is " << dataStatus << endl;
	cout << "Original size is: " << originalSize << " bytes" << endl;
	cout << "Compressed size is: " << compressedSize << " bytes" << endl;
	cout << "Compression ratio is: " << compressionRatio << "%";
}

// 1. Implement a sequential compression algorithm as a baseline.
// 2. Parallelize the compression process using OpenMP, distributing chunks of the across available threads.
// 3. Implement a distributed-memory version using MPI, where each process handles a portion of the input data and communicates compressed results to a master node.

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int my_rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
	bool testing = false, openMPlog = true, writeResultFiles = true;
	int numberOfThreads;
	string filename = "";
	
	if (my_rank == 0) 
	{
		if (argc < 3 && !testing){
			fprintf(stderr, "Usage: %s [to-compress-filename][number-of-threads]\n", argv[0]);
			numberOfThreads = -1;
		}else
		{
			if (testing)
			{
				string input = "";
				cout << "Please enter the amount of OpenMP threads to use: ";
				cin >> input;
				numberOfThreads = std::stoi(input);
				
				cout << "Please write the filename of what you'd like to compress (including type: '.txt/.raw/.bmp/etc.'): ";
				cin >> filename;
			}else
			{
				numberOfThreads = std::stoi(argv[2]);
				filename = argv[1];
			}
		}
	}
	
	//--We broadcast the number of threads we want to use for each processor
	MPI_Bcast(&numberOfThreads, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (numberOfThreads == -1) MPI_Abort(MPI_COMM_WORLD, 1);
			
	omp_set_num_threads(numberOfThreads);
	
	if (openMPlog){
		#pragma omp parallel
		{
		    int thread_id = omp_get_thread_num();
		    int num_threads = omp_get_num_threads();
		    if (thread_id == 0) {
		        cout << "Using " << num_threads << " threads for the algorithm [rank:" << my_rank << "]" << endl;
		    }
		}
	}
	
	// Website containing free raw pictures https://www.lapseoftheshutter.com/free-raw-landscape-images-for-retouching/ -> cr2 to bmp conversion needed!
	// Why doesn't this algorithm work with normal jpg/png/cr2? Short answer: those formats are already compressed (jpg/png) and/or ordered (cr2)
	
	huffmantool ht;
	string fileContent = "";
	bool sourceIsImage = false, fatalError = false;
	
	if (my_rank == 0) 
	{
		//--First we read the file
	    std::ifstream reader;
		reader.open(filename, std::ios::in | std::ios::binary);
	    if (!reader.is_open())
	    {
	        std::cerr << "ERROR: No such file exists or cannot open file " + filename;
	        fatalError = true;
	    }else
		{
			string line;
		    while (getline(reader, line)) {
		    	fileContent += line + "\n";
		    }
		    sourceIsImage = ht.isImageFile(filename);
		}
	}
	
	MPI_Bcast(&fatalError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (fatalError) MPI_Abort(MPI_COMM_WORLD, 1);
	
	//--We broadcast wether the source file is an image or not 
	MPI_Bcast(&sourceIsImage, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//--We allocate the variables for calculating the latency
	std::chrono::time_point<std::chrono::high_resolution_clock> startTime, endTime;
	double* times = (my_rank == 0) ? new double[EXECUTION_TIMES] : 0;
	
	int totalSize = fileContent.size();
	int chunkSize = totalSize / comm_size;
	int disparity = totalSize % comm_size;
	
	std::vector<int> rank_data_counts(comm_size), rank_offsets(comm_size), compressed_sizes(comm_size);
	int offset = 0;
	
	for (int i = 0; i < comm_size; i++) {
	    rank_data_counts[i] = chunkSize + (i < disparity ? 1 : 0);  
	    rank_offsets[i] = offset;
	    offset += rank_data_counts[i];
	}
	
	int localSize = rank_data_counts[my_rank]; 
	std::vector<char> localChunk(localSize); 
	
	//--We run for [EXECUTION_TIMES] times the compression algorithm
	std::vector<char> compressed_buffer(comm_size);
	string compressedData;
	for (int i = 0; i < EXECUTION_TIMES; i++){
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
		
		if (my_rank == 0) startTime = std::chrono::high_resolution_clock::now();
		
		MPI_Scatterv(fileContent.data(), rank_data_counts.data(), rank_offsets.data(), MPI_CHAR, localChunk.data(), localSize, MPI_CHAR, 0, MPI_COMM_WORLD);
		
		std::string localCompressed = ht.compressString(std::string(localChunk.begin(), localChunk.end()), sourceIsImage);
		
		int local_size = localCompressed.size();
		MPI_Gather(&local_size, 1, MPI_INT, compressed_sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		cout << "Local size of rank " << my_rank << ": " << local_size << endl;

	    if (my_rank == 0) {
	        int offset = 0;
		    for (int i = 0; i < comm_size; i++) {
		        rank_offsets[i] = offset;
		        offset += compressed_sizes[i];
		    }
		    compressedData.resize(offset);
	    
			cout << "Final calculated size: " << offset << endl;
	    }
	    
	    MPI_Abort(MPI_COMM_WORLD, 1);
	    
		MPI_Gatherv(localCompressed.data(), local_size, MPI_CHAR, compressed_buffer.data(), 
			compressed_sizes.data(), rank_offsets.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
		
		if (my_rank == 0)
		    compressedData = std::string(compressed_buffer.begin(), compressed_buffer.end());
		
		if (my_rank == 0) 
		{
			endTime = chrono::high_resolution_clock::now();
			times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
		}
	}
	
	double compression_percentile_value;
	
	if (my_rank == 0) 
	{
		//--We get the 90* percentile of the runs
	    std::sort(times, times + EXECUTION_TIMES);
	    int compression_percentile_index = 0.9 * (EXECUTION_TIMES - 1); 
	    compression_percentile_value = times[compression_percentile_index];
	}
	
	std::ofstream writer;
	
	if (my_rank == 0 && writeResultFiles) 
	{
		//--We write the compressed data
	    writer.open("compressed_"+filename.substr(0, filename.find_last_of("."))+".data", std::ios::out | std::ios::binary | std::ios::trunc);
	    writer.write(compressedData.data(), compressedData.size());
	    writer.close();
	}
    
	//--We run for [EXECUTION_TIMES] times the decompression algorithm
	string decompressed = "";
	if (my_rank == 0) // < AT THE MOMENT WE ARE JUST TESTING COMPRESSION, DECOMPRESSION WILL BE DONE ONLY BY RANK 0!!!
	{
	    for (int i = 0; i < EXECUTION_TIMES; i++){
			MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
			
			if (my_rank == 0) startTime = std::chrono::high_resolution_clock::now();
			
			decompressed = ht.decompressString(compressedData, sourceIsImage); // > We need to scatter this!!! <
			
			if (my_rank == 0) 
			{
				endTime = chrono::high_resolution_clock::now();
				times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
			}
		}
	}
	
	if (my_rank == 0)
	{
		//--We get the 90* percentile of the runs
	    std::sort(times, times + EXECUTION_TIMES);
	    int decompression_percentile_index = 0.9 * (EXECUTION_TIMES - 1); 
	    double decompression_percentile_value = times[decompression_percentile_index];
	
		delete[] times;
	
		cout << "Executed compression for [90* percentile]: " << compression_percentile_value / 1000000 << "ms" << endl;
		
		cout << "Executed decompression for [90* percentile]: " << decompression_percentile_value / 1000000 << "ms" << endl;
	}
	
	if (my_rank == 0 && writeResultFiles) 
	{
		//--We write the decompressed data
	    writer.open("decompressed_"+filename, std::ios::out | std::ios::binary | std::ios::trunc);
	    writer.write(decompressed.data(), decompressed.size());
	    writer.close();
	}
	
	if (my_rank == 0) AnalizeCompressionRatio(fileContent, compressedData, decompressed);
	
	MPI_Finalize();
}
