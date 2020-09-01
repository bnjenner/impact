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
#include "lib/annotations.h"
#include "lib/alignments.h"

// Threads
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

    //std::cerr << annotation.get_feature("Contig0-", 500, 600) << "\n";

    std::cerr << "[Counting Reads...]\n";
    alignment.get_counts(&annotation);

    std::cerr << "Features: " << alignment.noncounts[0] << "\n";
    std::cerr << "No Features: " << alignment.noncounts[1] << "\n";
    std::cerr << "Ambiguous: " << alignment.noncounts[2] << "\n";
    std::cerr << "Low Quality: " << alignment.noncounts[3] << "\n";
    std::cerr << "Unmapped: " << alignment.noncounts[4] << "\n";


    // // Poor attempt at asynchronous counting :(
    // Count Reads
    // std::cerr << "[Counting Reads...]\n";
    
    // futures vector and iterator
    // std::vector<std::future<void>> count_futures;
    // int i = 0;
    // while (i < annotation.indicies.size() - 1) {

    //     count_futures.push_back(std::async(std::launch::async, getCounts, &annotation, &alignment, 
    //                                        (annotation.indicies[i] + 1), annotation.indicies[i+1],
    //                                        args.peak_detection));


    //     i++;
    //}
   	

    // Close alignment file
    alignment.close();
   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}