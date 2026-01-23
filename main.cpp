#include <iostream>
#include <fstream>
#include <string>

// Huffman algorithm available on "https://github.com/bhumijgupta/huffman-compression-library"
#include "huffmantool.h"

using namespace std;

// 1. Implement a sequential compression algorithm as a baseline.


int main(int argc, char** argv)
{
	huffmantool ht;
	ht.compressFile("test.txt", "compressed.txt");
	ht.decompressFile("compressed.txt", "decompressed.txt");
	
	string compressedData = ht.compressString("Hello World!");
	cout << ht.decompressString(compressedData);

	//ht.benchmark("test.txt");
}
