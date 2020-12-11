#include <iostream>
#include <chrono> 
#include <string>
#include <vector>
#include <thread>
#include <future>
#include <math.h>
#include <fstream>
#include <armadillo>
#include "api/BamReader.h" 
#include "api/BamAux.h"
#include "lib/parser.h"
#include "lib/utils.h"
#include "lib/peaks.h"
#include "lib/graph.h"
//#include "lib/annotations.h"
#include "lib/alignments.h"

// Threads
void open_alignment(AlignmentFile *alignment) {
    alignment -> open();
}

// void open_annotation(AnnotationFile *annotation) {
//     annotation -> open();
// }

void count_thread(const ImpactArguments *args, int ref) {

    AlignmentFile alignment(args);
    alignment.open();
    alignment.get_counts(ref);
    alignment.close();
}


// Main 
int main(int argc, char const ** argv) {

	auto start = std::chrono::high_resolution_clock::now(); 

	ImpactArguments args;
	seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);


	// Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;
    }
    
    std::cerr << "[IMPACT]\n";

    std::cerr << "[Parsing Input Files...]\n";

    AlignmentFile alignment(&args);
    alignment.open();
    alignment.close(); 

    // Construct alignment object in thread
    // AlignmentFile alignment(&args);
    // std::thread align_thread(open_alignment, &alignment);

    // Construct annotation object in thread
    // AnnotationFile annotation(&args);
    // std::thread annotate_thread(open_annotation, &annotation);

    // join threads
    //align_thread.join();
    //annotate_thread.join();

    // int i = 0;
    // while (i < alignment.references.size()) {
    //     std::cerr << "[Counting from " << alignment.contig_cache[i] << "...]\n";
    //     alignment.get_counts(i);
    //     i ++;
    //     break;
    // }

    // Close alignment file
    //alignment.close();

    int i = 0;
    int n = 1; // alignment.references.size();


    while (i < n){

        std::cerr << "[Counting from " << alignment.contig_cache[i] << "...]\n";
        std::thread thread1(count_thread, &args, i);

        i++;

        // std::cerr << "[Counting from " << alignment.contig_cache[i] << "...]\n";
        // std::thread thread2(count_thread, &args, i);

        // i++; 

        // std::cerr << "[Counting from " << alignment.contig_cache[i] << "...]\n";
        // std::thread thread3(count_thread, &args, i);

        thread1.join();
        // thread2.join();
        // thread3.join();

        // i++;
    }

   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}