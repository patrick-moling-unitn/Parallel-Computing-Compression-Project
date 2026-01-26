#include <iostream>
#include <fstream>
#include <string>

// Huffman algorithm available on "https://github.com/bhumijgupta/huffman-compression-library"
#include "huffmantool.h"

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

int main(int argc, char** argv)
{
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
	string compressedData = ht.compressString(fileContent, sourceIsImage);
	
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
