#include <iostream>
#include <chrono> 
#include <string>
#include <vector>
#include <thread>
#include <future>
#include <math.h>
#include <fstream>
#include "lib/parser.h"
#include "lib/utils.h"
#include "lib/peaks.h"
#include "lib/alignments.h"
#include "lib/annotations.h"
#include "lib/count.h"

void open_alignment(AlignmentFile *alignment) {
    alignment -> open();
}

void open_annotation(AnnotationFile *annotation) {
    annotation -> open();
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

    // Construct alignment object in thread
    AlignmentFile alignment(&args);
    std::thread align_thread(open_alignment, &alignment);

    // Construct annotation object in thread
    AnnotationFile annotation(&args);
    std::thread annotate_thread(open_annotation, &annotation);

    // join threads
    align_thread.join();
    annotate_thread.join();


    // Count Reads
    std::cerr << "[Counting Reads...]\n";
   	int total_counts = getCounts(&annotation, &alignment, args.peak_detection);
   	

    // Close alignment file
    alignment.close();
   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}