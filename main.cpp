#include <iostream>
#include <fstream>
#include <string>

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

void AnalizeCompressionRatio(string original, string compressed){
	float originalSize = sizeof(original) * original.size(), compressedSize = sizeof(compressed) * compressed.size();
	float compressionRatio = (1 - compressedSize / originalSize) * 100;
	cout << "--------------------------------" << endl;
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
	
	/*
	ht.compressFile(filename, "compressed.data");
	ht.decompressFile("compressed.data", "decompressed_"+filename, filename);*/
	
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
	
	//*This is gonna be very useful later on for MPI! Just split the string into N chunks.
	int n = 10;
	double* times = new double[n];
	string compressedData;
	for (int i = 0; i < n; i++){
		auto start = std::chrono::high_resolution_clock::now();
		compressedData = ht.compressString(fileContent, sourceIsImage);
		auto end = chrono::high_resolution_clock::now();
		double time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
		times[i] = time;
	}
	
    std::sort(times, times + 10);
    int percentile_index = 0.9 * (n - 1);  // 90* percentile
    double percentile_value = times[percentile_index];
	
	cout << "Executed for [90* percentile]: " << percentile_value / 1000000 << "ms" << endl;
	
	cout << "Data size: " << compressedData.size() << endl;
	
	//--Then we write the compressed data
    std::ofstream writer;
    writer.open("compressed.data", std::ios::out | std::ios::binary | std::ios::trunc);
    writer.write(compressedData.data(), compressedData.size());
    writer.close();
    
    //*Also useful for MPI splitting
	string decompressed = ht.decompressString(compressedData, sourceIsImage);
	
	//--Then we write the decompressed data
    writer.open("decompressed_"+filename, std::ios::out | std::ios::binary | std::ios::trunc);
    writer.write(decompressed.data(), decompressed.size());
    writer.close();
	
	//cout << ht.decompressString(compressedData) << endl;
	AnalizeCompressionRatio(fileContent, compressedData);

	//ht.benchmark("test.txt");
}
