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
	float originalSize = (float)original.size(), compressedSize = (float)compressed.size();
	float compressionRatio = (1 - ((float)compressedSize / originalSize)) * 100.0f;
	cout << "--------------------------------" << endl;
	string dataStatus = original == decompressed ? "working as expected! [Success]" : "corrupting the original data! [Error]";
	cout << "Data decompression is " << dataStatus << endl;
	cout << "Original size is: " << originalSize << " bytes" << endl;
	cout << "Compressed size is: " << compressedSize << " bytes" << endl;
	cout << "Compression ratio is: " << compressionRatio << "%" << endl;
	cout << "--------------------------------" << endl;
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
	
	bool testing = false, openMPlog = true, writeResultFiles = true, debuggingLogs = true;
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
	int sourceIsImage = 0, fatalError = 0;
	
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
		    	if (fileContent.size() > 0) fileContent += "\n";
		    	fileContent += line;
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
	
	vector<int> rank_chunk_sizes(comm_size), rank_chunk_offsets(comm_size), compressed_sizes(comm_size), compressed_offsets(comm_size);
	
	//--We calculate the chunks to scatter
	if (my_rank == 0)
	{
		int totalSize = fileContent.size();
		int chunkSize = totalSize / comm_size;
		int disparity = totalSize % comm_size;
		
		int offset = 0;
		
		for (int i = 0; i < comm_size; i++) {
		    rank_chunk_sizes[i] = chunkSize + (i < disparity ? 1 : 0);  
		    rank_chunk_offsets[i] = offset;
		    offset += rank_chunk_sizes[i];
		}
	}
	
	//--We share just the chunk size to 
	MPI_Bcast(rank_chunk_sizes.data(), comm_size, MPI_INT, 0, MPI_COMM_WORLD);
	
	int localSize = rank_chunk_sizes[my_rank]; 
	vector<char> localChunk(localSize, 0);
	
	//--We run for [EXECUTION_TIMES] times the compression algorithm
	vector<char> compressed_buffer;
	string compressedData;
	for (int i = 0; i < EXECUTION_TIMES; i++){
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
		
		if (my_rank == 0) startTime = std::chrono::high_resolution_clock::now();
		
		MPI_Scatterv(fileContent.data(), rank_chunk_sizes.data(), rank_chunk_offsets.data(), MPI_CHAR, 
					 localChunk.data(), localSize, MPI_CHAR, 0, MPI_COMM_WORLD);
		
		string localCompressed = ht.compressString(std::string(localChunk.begin(), localChunk.end()), sourceIsImage) + '\x1E';
		
		int local_compressed_size = localCompressed.size();
		MPI_Gather(&local_compressed_size, 1, MPI_INT, compressed_sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	    if (my_rank == 0) {
	        int offset = 0;
		    for (int i = 0; i < comm_size; i++) {
		        compressed_offsets[i] = offset;
		        offset += compressed_sizes[i];
		    }
		    compressed_buffer.resize(offset);
	    }
	    
		MPI_Gatherv(localCompressed.data(), local_compressed_size, MPI_CHAR, 
			compressed_buffer.data(), compressed_sizes.data(), compressed_offsets.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
		
		if (my_rank == 0) 
		{
			endTime = chrono::high_resolution_clock::now();
			times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
		    compressedData = std::string(compressed_buffer.begin(), compressed_buffer.end());
		}
		
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
		if (my_rank == 0 && debuggingLogs) 
		{
			cout << ">>> Iteration [" << i << "], compressed data: ";
			for (int i=0; i < compressedData.size(); i++){
				cout << ((int)compressedData[i]) << " ";
			}
			cout << endl;
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
	
	rank_chunk_sizes.resize(comm_size);
	rank_chunk_offsets.resize(comm_size);
	vector<int> decompressed_sizes(comm_size), decompressed_offsets(comm_size);
	
	if (my_rank == 0)
	{
		size_t index_of;
		
		int offset = 0;
		
		for (int i = 0; i < comm_size; i++) {
    		index_of = compressedData.find('\x1E', offset);
    		
			if (index_of != std::string::npos)
			{
	    		compressedData.erase(index_of, 1);
	    		
			    rank_chunk_sizes[i] = index_of - offset;  
			    rank_chunk_offsets[i] = offset;
				//cout << "compressedData: " << compressedData << "\tindexOf: " << index_of << "\toffset: " << offset << "\tsize: "<< rank_chunk_sizes[i] << endl;
			}else if (i == comm_size-1)
			{
			    rank_chunk_sizes[i] = compressedData.size() - offset;
			}else 
			{
	        	std::cerr << "ERROR: The chunk calculation on the decompression produced an exception!";
	        	MPI_Abort(MPI_COMM_WORLD, 1);
			}
			offset += rank_chunk_sizes[i];
		}
	}
	
	MPI_Bcast(rank_chunk_sizes.data(), comm_size, MPI_INT, 0, MPI_COMM_WORLD);
	
	localSize = rank_chunk_sizes[my_rank]; 
	localChunk.resize(localSize, 0);
    
	//--We run for [EXECUTION_TIMES] times the decompression algorithm
	vector<char> decompressed_buffer;
	string decompressedData = "";
    for (int i = 0; i < EXECUTION_TIMES; i++){
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
		
		if (my_rank == 0) startTime = std::chrono::high_resolution_clock::now();
		
		MPI_Scatterv(compressedData.data(), rank_chunk_sizes.data(), rank_chunk_offsets.data(), MPI_CHAR, 
					 localChunk.data(), localSize, MPI_CHAR, 0, MPI_COMM_WORLD);
		
		string localDecompressed = ht.decompressString(std::string(localChunk.begin(), localChunk.end()), sourceIsImage); 
	    
	    int local_decompressed_size = localDecompressed.size();
		MPI_Gather(&local_decompressed_size, 1, MPI_INT, decompressed_sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	    if (my_rank == 0) {
	        int offset = 0;
		    for (int i = 0; i < comm_size; i++) {
		        decompressed_offsets[i] = offset;
		        offset += decompressed_sizes[i];
		    }
		    decompressed_buffer.resize(offset, 0);
	    }
	    
		MPI_Gatherv(localDecompressed.data(), local_decompressed_size, MPI_CHAR, 
			decompressed_buffer.data(), decompressed_sizes.data(), decompressed_offsets.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
		
		if (my_rank == 0) 
		{
			endTime = chrono::high_resolution_clock::now();
			times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
		    decompressedData = std::string(decompressed_buffer.begin(), decompressed_buffer.end());
		}
		
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
		if (my_rank == 0 && debuggingLogs) 
			cout << ">>> Iteration [" << i << "], decompressed data:" << decompressedData << endl;
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
	    writer.write(decompressedData.data(), decompressedData.size());
	    writer.close();
	}
	
	if (my_rank == 0) AnalizeCompressionRatio(fileContent, compressedData, decompressedData);
	
	MPI_Barrier(MPI_COMM_WORLD); // Wait for all the processes
	MPI_Finalize();
}
