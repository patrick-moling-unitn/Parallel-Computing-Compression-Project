#include <iostream>
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
	ht.compressFile(filename, "compressed.data");
	ht.decompressFile("compressed.data", "decompressed_"+filename, filename);
	
	string originalData = "Hello World!";
	string compressedData = ht.compressString(originalData);
	cout << ht.decompressString(compressedData) << endl;
	AnalizeCompressionRatio(originalData, compressedData);

	//ht.benchmark("test.txt");
}
