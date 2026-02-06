// The MIT License (MIT)

// Copyright (c) 2020 BHUMIJ GUPTA

//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.

/* code */

// [1] See Wiki for more info
#ifndef HUFFMAN_TOOL_H
#define HUFFMAN_TOOL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <queue>
#include <string>
#include <bitset>
#include <chrono>
#include <sstream>

// ******************************
// 	CUSTOM INCLUDES & MACROS
// ******************************

#include "src/scanner.h"
#include "src/cfp.h"

#define cfp charFreqPair

// ******************************
// 	CLASS DECLARATION
// ******************************

/**
 * @brief This class contains methods for file compression and decompression
 * 
 */
class huffmantool
{
    void traverse(cfp const *, std::unordered_map<char, std::string> &, std::string const);
    void traverse(cfp const *, std::unordered_map<std::string, char> &, std::string const);
    cfp *readTree(std::istream &);
    void writeTree(std::ostream &, cfp const *);
    void prettyPrint(std::string const);
    void prettyPrint(int const);
    void printSeparator();
    int lposSlash(std::string const);

public:
    std::string compressFile(std::string, std::string);
    std::string decompressFile(std::string, std::string, std::string);
    //>> Added for future MPI implementations (Division for chuncks)
    std::string compressString(const std::string &input, bool sourceIsImage);
    std::string decompressString(const std::string &input, bool sourceIsImage);
    //<<
    void benchmark(std::string);
    
    bool isImageFile(const std::string& filename);
};

// ******************************
// 	METHODS DEFINATION
// ******************************

void huffmantool::traverse(cfp const *head, std::unordered_map<char, std::string> &charKeyMap, std::string const s)
{
    if (head->left == NULL && head->right == NULL)
    {
        charKeyMap[head->getChar()] = s;
        return;
    }
    traverse(head->left, charKeyMap, s + "0");
    traverse(head->right, charKeyMap, s + "1");
}

void huffmantool::traverse(cfp const *head, std::unordered_map<std::string, char> &keyCharMap, std::string const s)
{
    if (head->left == NULL && head->right == NULL)
    {
        keyCharMap[s] = head->getChar();
        return;
    }
    traverse(head->left, keyCharMap, s + "0");
    traverse(head->right, keyCharMap, s + "1");
}

cfp *huffmantool::readTree(std::istream &reader)
{
    uint8_t marker;
    reader.read(reinterpret_cast<char*>(&marker), 1);

    if (marker == 1)
    {
        uint8_t c;
        reader.read(reinterpret_cast<char*>(&c), 1);
        return new cfp(static_cast<char>(c), 0);
    }

    cfp *head = new cfp('~', 0);
    head->left = readTree(reader);
    head->right = readTree(reader);
    return head;
}

void huffmantool::writeTree(std::ostream& writer, const cfp* head)
{
    if (!head->left && !head->right)
    {
        uint8_t marker = 1;
        uint8_t c = static_cast<uint8_t>(head->getChar());
        writer.write(reinterpret_cast<char*>(&marker), 1);
        writer.write(reinterpret_cast<char*>(&c), 1);
        return;
    }
    uint8_t marker = 0;
    writer.write(reinterpret_cast<char*>(&marker), 1);
    writeTree(writer, head->left);
    writeTree(writer, head->right);
}

void huffmantool::prettyPrint(std::string const out)
{
    // [5]
    std::cout << std::left << std::setw(30) << out;
}
void huffmantool::prettyPrint(int const out)
{
    std::cout << std::left << std::setw(30) << out;
}
void huffmantool::printSeparator()
{
    for (int i = 0; i < 80; i++)
        std::cout << "-";
    std::cout << "\n";
}

int huffmantool::lposSlash(std::string const filename)
{
    int pos = -1;
    for (unsigned int i = 0; i < filename.length(); i++)
    {
        if (filename[i] == '/')
            pos = i;
    }
    return pos;
}

bool huffmantool::isImageFile(const std::string& filename) {
    std::string extension = filename.substr(filename.find_last_of(".") + 1);
    for (auto &character : extension) character = std::tolower(character);
    return (extension == "bmp" || extension == "raw" || extension == "pgm" || extension == "ppm");
}

/**
 * @brief Compress the given file
 * 
 * @param sourceFile Path of the file to be compressed (relative or absolute)
 * @param compressedFile Path of the compressed file. Default: path/"compressed_"+sourceFile.extension
 * @return std::string Returns the path of the compressed file. If some error occurs, returns empty string
 */
std::string huffmantool::compressFile(std::string sourceFile, std::string compressedFile = "")
{
    if (compressedFile == "")
    {
        int pos = lposSlash(sourceFile);
        compressedFile = sourceFile.substr(0, pos + 1) + "compressed_" + sourceFile.substr(pos + 1);
    }
    std::ifstream reader;
    reader.open(sourceFile, std::ios::in | std::ios::binary);
    if (!reader.is_open())
    {
        std::cerr << "ERROR: No such file exists or cannot open file " + sourceFile;
        return "";
    }
    bool sourceIsImage = isImageFile(sourceFile);
    // map for index in vector
    std::unordered_map<char, int> *index = new std::unordered_map<char, int>;
    std::vector<cfp *> freq_store;
    char ch;
    int numChars = 0;
    // [2]
    if (sourceIsImage){
    	char prev = 0;
	    while (reader.read(&ch, 1)) {
	        char value = ch - prev;
	        prev = ch;
	
	        numChars++;
	        if (index->count(value) > 0) {
	            int ind = index->at(value);
	            freq_store[ind]->setFreq(freq_store[ind]->getFreq() + 1);
	        } else {
	            index->insert({value, freq_store.size()});
	            freq_store.push_back(new cfp(value, 1));
	        }
	    }
	}else
	{
	    while (reader.read(&ch, 1))
	    {
	        numChars++;
	        if (index->count(ch) > 0)
	        {
	            int ind = index->at(ch);
	            freq_store[ind]->setFreq(freq_store[ind]->getFreq() + 1);
	        }
	        else
	        {
	            index->insert({ch, freq_store.size()});
	            freq_store.push_back(new cfp(ch, 1));
	        }
	    }
	}
    reader.close();
    delete index;
    if (freq_store.size() <= 1)
    {
        std::cout << "INFO: No need for encryption\n";
        return "";
    }
    // Create min priority queue for cfp pair
    std::priority_queue<cfp *, std::vector<cfp *>, pairComparator> freq_sorted;
    for (auto i : freq_store)
    {
        freq_sorted.push(i);
    }
    cfp *head;
    while (!freq_sorted.empty())
    {
        cfp *first = freq_sorted.top();
        freq_sorted.pop();
        cfp *second = freq_sorted.top();
        freq_sorted.pop();
        cfp *newPair = new cfp('~', first->getFreq() + second->getFreq());
        newPair->left = first;
        newPair->right = second;
        freq_sorted.push(newPair);
        if (freq_sorted.size() == 1)
        {
            head = newPair;
            break;
        }
    }
    std::unordered_map<char, std::string> charKeyMap;
    traverse(head, charKeyMap, "");

    // Read from file and write compressed to new file
    std::ofstream writer;
    writer.open(compressedFile, std::ios::out | std::ios::binary | std::ios::trunc);
    reader.open(sourceFile, std::ios::in | std::ios::binary);

    // First write the tree in preorder form
    writeTree(writer, head);
    // delete head;
    // Write numChars to check file integrity
    writer.write(reinterpret_cast<char*>(&numChars), sizeof(numChars));
    char chr = 0;
    int bufferSize = 8;
    
    if (sourceIsImage){
    	char prev = 0;
	    while (reader.read(&ch, 1)) {
	        char value = ch - prev;
	        prev = ch;
	
	        std::string bin = charKeyMap[value];
	        for (unsigned int i = 0; i < bin.length(); i++)
	        {
	            chr = (chr << 1) ^ (bin[i] - '0');
	            bufferSize--;
	            if (bufferSize == 0)
	            {
	                writer.write(&chr, 1);
	                chr = 0;
	                bufferSize = 8;
	            }
	        }
	    }
	}else
	{
	    while (reader.read(&ch, 1))
	    {
	        std::string bin = charKeyMap[ch];
	        for (unsigned int i = 0; i < bin.length(); i++)
	        {
	            chr = (chr << 1) ^ (bin[i] - '0');
	            bufferSize--;
	            if (bufferSize == 0)
	            {
	                writer.write(&chr, 1);
	                chr = 0;
	                bufferSize = 8;
	            }
	        }
	    }
	}
    if (bufferSize)
    {
        chr = chr << bufferSize;
        writer.write(&chr, 1);
    }
    writer.close();
    reader.close();
    return compressedFile;
}

/**
 * @brief Decompress the given file
 * 
 * @param compressedFile Path of the file to be decompressed (absolute or relative)
 * @param retrievedFile Path of the decompressed file. Default: path/"decompressed_"+compressedFile.extension
 * @return std::string Returns the path of the retrieved file. If some error occurs, returns empty string
 */
std::string huffmantool::decompressFile(std::string compressedFile, std::string retrievedFile = "", std::string originalFileName = "")
{
    bool sourceIsImage = isImageFile(originalFileName);
    if (retrievedFile == "")
    {
        int pos = lposSlash(compressedFile);
        retrievedFile = compressedFile.substr(0, pos + 1) + "decompressed_";
        if (compressedFile.length() - pos + 1 >= 11 && compressedFile.substr(pos + 1, 11) == "compressed_")
            retrievedFile += compressedFile.substr(pos + 11 + 1);
        else
            retrievedFile += compressedFile.substr(pos + 1);
    }
    std::ifstream reader;
    std::ofstream writer;
    reader.open(compressedFile, std::ios::in | std::ios::binary);
    if (!reader.is_open())
    {
        std::cerr << "ERROR: No such file exists or cannot open file " + compressedFile;
        return "";
    }
    // create huffman tree from file
    cfp *head = readTree(reader);
    // Create key char map for decompression
    std::unordered_map<std::string, char> keyCharMap;
    traverse(head, keyCharMap, "");
    // delete head;
    // Read total number of characters
    int totalChars;
	reader.read(reinterpret_cast<char*>(&totalChars), sizeof(totalChars));
    writer.open(retrievedFile, std::ios::out | std::ios::binary | std::ios::trunc);
    std::string key = "";
    int readChars = 0;
    char ch;
    if (sourceIsImage){
    	char prev = 0;
	    while (reader.read(&ch, 1) && readChars != totalChars) {
	        std::string bin_read = std::bitset<8>(static_cast<unsigned char>(ch)).to_string();
	        for (unsigned int i = 0; i < bin_read.length(); i++) {
	            key += bin_read[i];
	            if (keyCharMap.count(key) > 0) {
	                char value = keyCharMap[key];
	                char out = value + prev;  
	                prev = out;
	                writer.write(&out, 1);
	                key = "";           
	                readChars++;
	                if (readChars == totalChars)
	                    break;
	            }
	        }
	    }
	}else 
	{
		while (reader.read(&ch, 1) && readChars != totalChars)
	    {
	        std::string bin_read = std::bitset<8>(static_cast<unsigned char>(ch)).to_string();
	        for (unsigned int i = 0; i < bin_read.length(); i++)
	        {
	            key += bin_read[i];
	            if (keyCharMap.count(key) > 0)
	            {
	            	char out = keyCharMap[key];
					writer.write(&out, 1);
	                key = "";
	                readChars++;
	                if (readChars == totalChars)
	                    break;
	            }
	        }
	    }
	}
    
    reader.close();
    writer.close();
    if (readChars != totalChars)
    {
        std::cerr << "ERROR: Compressed file is corrupted\n";
        return "";
    }
    return retrievedFile;
};

//ChatGPT Helped in rewriting the compression just for a string - very good understanding of streams was the key
std::string huffmantool::compressString(const std::string &input, bool sourceIsImage)
{
    // Mappa per l'indice nel vettore
    std::unordered_map<char, int> *index = new std::unordered_map<char, int>;
    std::vector<cfp *> freq_store;
    char ch;
    int numChars = 0;
    
    std::cout << "A";
    
    // Passo 1: Creiamo la mappa delle frequenze
    if (sourceIsImage){
		#pragma omp parallel
		{
		    std::unordered_map<char,int> local_index;
		    std::vector<cfp*> local_freq_store;  
	    	char prev = 0;
    
			#pragma omp for private(prev, local_freq_store)
	        for	(int i=0; i < input.length(); i++)
	        {
	            char value = input[i] - prev;
	            prev = input[i];
	
				#pragma omp atomic
	            numChars++;
	            if (index->count(value) > 0) {
	                int ind = index->at(value);
	                local_freq_store[ind]->setFreq(local_freq_store[ind]->getFreq() + 1);
	            } else {
	            	#pragma omp critical
					{
	                	index->insert({value, local_freq_store.size()});
					}
	                local_freq_store.push_back(new cfp(value, 1));
	            }
	        }
    	}
    }else
    { 
		#pragma omp parallel
		{
		    std::unordered_map<char,int> local_index;
		    std::vector<cfp*> local_freq_store;  
		    
    		#pragma omp single
	        {
    			std::cout << "B";
    		}	
    
			#pragma omp for
	        for	(int i=0; i < input.length(); i++)
	        {
	        	char value = input[i];
				#pragma omp atomic
	            numChars++;
	            
		        auto search = local_index.find(value);
		        if (search != local_index.end()) {
            		int ind = search->second;
	                local_freq_store[ind]->setFreq(local_freq_store[ind]->getFreq() + 1);
	            } else {
            		int pos = local_freq_store.size();
		            local_freq_store.push_back(new cfp(value, 1));
		            local_index.insert({value, pos});
	            }
	        }
	        
    		#pragma omp single
	        {
    			std::cout << "C";
    		}	
	        
	        #pragma omp critical
		    {
		    	//Joining back local results
		        for (auto& kv : local_index) {
		            char value = kv.first;
		            int pos = kv.second;
		            if (index->count(value) > 0) {
		                int global_pos = index->at(value);
		                freq_store[global_pos]->setFreq(freq_store[global_pos]->getFreq() + local_freq_store[pos]->getFreq());
		                delete local_freq_store[pos];
		            } else {
		                int global_pos = freq_store.size();
		                index->insert({value, global_pos});
		                freq_store.push_back(local_freq_store[pos]);
		            }
		        }
		    }
		}
    }
    
    std::cout << "D";
    printf("\n\n\n");
	    
	delete index;

    if (freq_store.size() <= 1) {
        std::cout << "INFO: No need for encryption\n";
        return "";  // Se c'e' solo un carattere, non serve compressione
    }

    // Creazione della coda di priorita' per i nodi di Huffman
    std::priority_queue<cfp *, std::vector<cfp *>, pairComparator> freq_sorted;
    
    for (auto i : freq_store) {
        freq_sorted.push(i);
    }

    cfp *head;
	#pragma omp parallel
	{
		while (!freq_sorted.empty()) {
	        cfp *first = freq_sorted.top();
			#pragma omp critical
			{
		        freq_sorted.pop();
			}
	        cfp *second = freq_sorted.top();
			#pragma omp critical
			{
	        	freq_sorted.pop();
			}
	        cfp *newPair = new cfp('~', first->getFreq() + second->getFreq());
	        newPair->left = first;
	        newPair->right = second;
	        
	        #pragma omp critical
			{
	        	freq_sorted.push(newPair);
			}
	        if (freq_sorted.size() == 1) {
	            head = newPair;
	            break;
	        }
	    }
	}

    std::unordered_map<char, std::string> charKeyMap;
    traverse(head, charKeyMap, "");

    // Passo 2: Comprimiamo la stringa e memorizziamo il risultato
    std::ostringstream compressedData(std::ios::binary);  // Usare std::ostringstream per binario
    
    // Scriviamo l'albero in formato binario (pre-order)
    writeTree(compressedData, head);
    
    // Scriviamo il numero di caratteri
    compressedData.write(reinterpret_cast<char*>(&numChars), sizeof(numChars));
    
    char chr = 0;
    int bufferSize = 8;

    // Scriviamo la stringa compressa nel stringstream
    if (sourceIsImage){
        char prev = 0;
		#pragma omp parallel for private(prev, bufferSize, chr)
        for (unsigned int i = 0; i < input.length(); i++)
        {
            char value = input[i] - prev;
            prev = input[i];

            const std::string &bin = charKeyMap[value];
            for (unsigned int j = 0; j < bin.length(); j++)
            {
                chr = (chr << 1) | (bin[j] - '0');
                bufferSize--;
                if (bufferSize == 0)
                {
					#pragma omp critical
					{
                    	compressedData.write(&chr, 1);
					}
                    chr = 0;
                    bufferSize = 8;
                }
            }
        }
    }else
	{
		#pragma omp parallel for private(bufferSize, chr)
        for (unsigned int i = 0; i < input.length(); i++) {
            std::string bin = charKeyMap[input[i]];
            for (unsigned int j = 0; j < bin.length(); j++) {
                chr = (chr << 1) ^ (bin[j] - '0');
                bufferSize--;
                if (bufferSize == 0) {
					#pragma omp critical
					{
                    	compressedData.write(&chr, 1);
					}
                    chr = 0;
                    bufferSize = 8;
                }
            }
        }
    }
    
    // Aggiungiamo l'eventuale byte residuo
    if (bufferSize) {
        chr = chr << bufferSize;
        compressedData.write(&chr, 1);
    }

    return compressedData.str();  // Restituiamo la stringa compressa
}

//ChatGPT Helped in rewriting the compression just for a string - very good understanding of streams was the key
std::string huffmantool::decompressString(const std::string &compressedInput, bool sourceIsImage)
{
    std::stringstream reader(compressedInput, std::ios::in | std::ios::binary);
    std::stringstream output;

    // Ricreiamo l'albero Huffman
    cfp *head = readTree(reader);

    // Creiamo la mappa chiave->char
    std::unordered_map<std::string, char> keyCharMap;
    traverse(head, keyCharMap, "");

    // Leggiamo il numero totale di caratteri
    int totalChars;
    reader.read(reinterpret_cast<char*>(&totalChars), sizeof(totalChars));

    std::string key = "";
    int readChars = 0;
    char ch;

	// We can't parallelize this since Huffamn decoding works serial bit after bit
    if (sourceIsImage){
        char prev = 0;
        
        while (reader.read(&ch, 1) && readChars < totalChars)
        {
            std::string bin_read =
                std::bitset<8>(static_cast<unsigned char>(ch)).to_string();

            for (char bit : bin_read)
            {
                key += bit;
                if (keyCharMap.count(key))
                {
                    char value = keyCharMap[key];
                    char out = value + prev;
                    prev = out;

                    output.write(&out, 1);
                    key.clear();
                    readChars++;

                    if (readChars == totalChars)
                        break;
                }
            }
        }
    }else 
    {
        while (reader.read(&ch, 1) && readChars < totalChars)
        {
            std::string bin_read = std::bitset<8>(static_cast<unsigned char>(ch)).to_string();
            
            for (unsigned int i = 0; i < bin_read.length(); i++)
            {
                key += bin_read[i];
                if (keyCharMap.count(key) > 0)
                {
                    char out = keyCharMap[key];
                    output.write(&out, 1);
                    key = "";
                    readChars++;
                    if (readChars == totalChars)
                        break;
                }
            }
        }
    }

    if (readChars != totalChars)
    {
        std::cerr << "ERROR: Compressed string is corrupted\n";
        return "";
    }

    return output.str();
}

/**
 * @brief Benchmarks the performance of the tool. Returns time taken for each step and compression ratio achieved.
 * 
 * @param sourceFile Path of the file to be used as source for compression and decompression
 */
void huffmantool::benchmark(std::string sourceFile)
{
    // [4]
    auto start_compression = std::chrono::high_resolution_clock::now();
    std::string compressedFile = compressFile(sourceFile);
    auto end_compression = std::chrono::high_resolution_clock::now();
    if (compressedFile == "")
        return;
    auto compression_time = std::chrono::duration_cast<std::chrono::microseconds>(end_compression - start_compression);

    auto start_decompression = std::chrono::high_resolution_clock::now();
    std::string retrievedFile = decompressFile(compressedFile, sourceFile);
    if (retrievedFile == "")
        return;
    auto end_decompression = std::chrono::high_resolution_clock::now();
    auto decompression_time = std::chrono::duration_cast<std::chrono::microseconds>(end_decompression - start_decompression);

    // Get file sizes
    scanner sc;
    int original_size = sc.getFileSize(sourceFile);
    if (original_size == -1)
        return;
    int compressed_size = sc.getFileSize(compressedFile);
    // If file not present or cannot open, it returns -1
    if (compressed_size == -1)
        return;
    int decompressed_size = sc.getFileSize(retrievedFile);
    // If file not present or cannot open, it returns -1
    if (decompressed_size == -1)
        return;

    printSeparator();
    std::cout << "                              B E N C H M A R K\n";
    printSeparator();
    std::cout << "\n";

    prettyPrint("Filetype");
    prettyPrint("Filename");
    prettyPrint("Filesize in bytes");
    std::cout << "\n\n";
    prettyPrint("Original");
    sourceFile = sourceFile.substr(lposSlash(sourceFile) + 1);
    prettyPrint(sourceFile);
    prettyPrint(original_size);
    std::cout << "\n";
    prettyPrint("Compressed");
    compressedFile = compressedFile.substr(lposSlash(compressedFile) + 1);
    prettyPrint(compressedFile);
    prettyPrint(compressed_size);
    std::cout << "\n";
    prettyPrint("Decompressed");
    retrievedFile = retrievedFile.substr(lposSlash(retrievedFile) + 1);
    prettyPrint(retrievedFile);
    prettyPrint(decompressed_size);

    float compression = 100.0 - ((float)compressed_size / original_size) * 100.0;
    std::cout << "\n\n";
    printSeparator();
    std::cout << "Time taken to compress file: " << compression_time.count() << " microseconds\n";
    std::cout << "Time taken to decompress file: " << decompression_time.count() << " microseconds\n";
    std::cout << "Compression: " << compression << "% \n\n";
};

#endif // HUFFMAN_TOOL_H
