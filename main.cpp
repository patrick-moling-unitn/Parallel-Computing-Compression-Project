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

int main(int argc, char** argv)
{
	bool testing = true;
	
	if (argc < 3 && !testing){
		fprintf(stderr, "Usage: %s [to-compress-filename][number-of-threads]\n", argv[0]);
		exit(1);
	}else{
		int targetThreads;
		if (testing)
		{
			string input = "";
			cout << "Please enter the amount of OpenMP threads to use: ";
			cin >> input;
			targetThreads = std::stoi(input);
		}
		const int numberOfThreads = testing ? targetThreads : std::stoi(argv[2]);
		omp_set_num_threads(numberOfThreads);
	}
	
	#pragma omp parallel
	{
	    int thread_id = omp_get_thread_num();
	    int num_threads = omp_get_num_threads();
	    if (thread_id == 0) {
	        cout << "Using " << num_threads << " threads for the algorithm" << endl;
	    }
	}
	
	string filename = "";
	cout << "Please write the filename of what you'd like to compress (including type: '.txt/.raw/.bmp/etc.'): ";
	cin >> filename;
	
	// Website containing free raw pictures https://www.lapseoftheshutter.com/free-raw-landscape-images-for-retouching/ -> cr2 to bmp conversion needed!
	// Why doesn't this algorithm work with normal jpg/png/cr2? Short answer: those formats are already compressed (jpg/png) and/or ordered (cr2)
	
	huffmantool ht;
	
	//--First we read the file
    std::ifstream reader;
	reader.open(filename, std::ios::in | std::ios::binary);
    if (!reader.is_open())
    {
        std::cerr << "ERROR: No such file exists or cannot open file " + filename;
        return 0;
    } 
	string line, fileContent = "";
    while (getline(reader, line)) {
    	fileContent += line + "\n";
    }
    bool sourceIsImage = ht.isImageFile(filename);
	
	//--We allocate the variables for calculating the latency
	std::chrono::time_point<std::chrono::high_resolution_clock> startTime, endTime;
	double* times = new double[EXECUTION_TIMES];
	
	//--We run for [EXECUTION_TIMES] times the compression algorithm
	string compressedData;
	for (int i = 0; i < EXECUTION_TIMES; i++){
		startTime = std::chrono::high_resolution_clock::now();
		compressedData = ht.compressString(fileContent, sourceIsImage);
		endTime = chrono::high_resolution_clock::now();
		times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
	}
	
	//--We get the 90* percentile of the runs
    std::sort(times, times + EXECUTION_TIMES);
    int compression_percentile_index = 0.9 * (EXECUTION_TIMES - 1); 
    double compression_percentile_value = times[compression_percentile_index];
	
	//--We write the compressed data
    std::ofstream writer;
    writer.open("compressed.data", std::ios::out | std::ios::binary | std::ios::trunc);
    writer.write(compressedData.data(), compressedData.size());
    writer.close();
    
	//--We run for [EXECUTION_TIMES] times the decompression algorithm
	string decompressed;
    for (int i = 0; i < EXECUTION_TIMES; i++){
		startTime = std::chrono::high_resolution_clock::now();
		decompressed = ht.decompressString(compressedData, sourceIsImage);
		endTime = chrono::high_resolution_clock::now();
		times[i] = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
	}
	
	//--We get the 90* percentile of the runs
    std::sort(times, times + EXECUTION_TIMES);
    int decompression_percentile_index = 0.9 * (EXECUTION_TIMES - 1); 
    double decompression_percentile_value = times[decompression_percentile_index];
	
	delete[] times;
	
	cout << "Executed compression for [90* percentile]: " << compression_percentile_value / 1000000 << "ms" << endl;
	
	cout << "Executed decompression for [90* percentile]: " << decompression_percentile_value / 1000000 << "ms" << endl;
	
	//--We write the decompressed data
    writer.open("decompressed_"+filename, std::ios::out | std::ios::binary | std::ios::trunc);
    writer.write(decompressed.data(), decompressed.size());
    writer.close();
	
	AnalizeCompressionRatio(fileContent, compressedData, decompressed);
}
