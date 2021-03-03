#include <iostream>
#include <chrono> 
#include <string>
#include <vector>
#include <thread>
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
// void open_alignment(AlignmentFile *alignment) {
//     alignment -> open();
// }

// void open_annotation(AnnotationFile *annotation) {
//     annotation -> open();
// }

void count_thread(AlignmentFile *alignment, int ref) {
    alignment -> open();
    alignment -> get_counts(ref);
    alignment -> close();
}


// Main 
int main(int argc, char const ** argv) {

    std::cerr << "[IMPACT]\n";

	auto start = std::chrono::high_resolution_clock::now(); 

    // Parse arguments
	ImpactArguments args;
	seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);


	// Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;
    }
    

    // Parse input files
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


    // Number of contigs for subdividing work
    int n = alignment.references.size();

    // Create vector of objects for multithreading
    std::vector<AlignmentFile*> alignments;
    for (int i = 0; i < n; i++) {
        alignments.push_back(new AlignmentFile(&args));
    }


    // initialize thread vector and process numbers
    int i = 0;
    int proc = args.threads;
    std::vector<std::thread> threads;

    // Multithreaded processing of alignmetns
    std::cerr << "[Processing Alignments...]\n";
    while (i < n) {

        proc = std::min(proc, n - i);

        // Add threads to thread vector
        for (int j = 0; j < proc; j++) {
            threads.push_back(std::thread(count_thread, alignments[i + j], i + j));
        }

        // Join threads
        for (auto &th : threads) {
            th.join();
        }

        threads.clear();
        i += proc;
    }

    
    // Report counts (this is single threaded for order reasons)
    std::cerr << "[Writing Results to STDOUT...]\n";
    for (int i = 0; i < n; i++) {
        alignments[i] -> print_counts();
    }


   	// The most complicated line of "get the time" I have ever seen. 
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}