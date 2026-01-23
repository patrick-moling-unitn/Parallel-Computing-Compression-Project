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
	huffmantool ht;
	ht.compressFile("test.txt", "compressed.txt");
	ht.decompressFile("compressed.txt", "decompressed.txt");
	
	string originalData = "Hello World!";
	string compressedData = ht.compressString(originalData);
	cout << ht.decompressString(compressedData) << endl;
	AnalizeCompressionRatio(originalData, compressedData);

	//ht.benchmark("test.txt");
}
