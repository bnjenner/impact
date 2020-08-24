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


    // // Count Reads
    // std::cerr << "[ Counting Reads... ]\n";
   	// int total_counts = getCounts(annotation, alignment, args.peak_detection);
   	

    // Close alignment file
    alignment.close();
   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

   
	// std::cerr << "-----------------\nProcess Complete!\n";
	// std::cerr << "Runtime: " << duration.count() << " seconds.\n";


    // std::ofstream out_file;
    // out_file.open("impact_stats.txt");

    // out_file << "----------\nIMPACT\n----------\n";
    // out_file << "Parameters:\n";
    // out_file << "  Alignment File: " << args.alignment_file << "\n";
    // out_file << "  Index File: " << args.index_file << "\n";
    // out_file << "  GFF File: " << args.gff_file << "\n";
    // out_file << "  Strandedness: " << args.strandedness << "\n";
    // out_file << "  Library: " << args.library_type << "\n";
    // out_file << "  Nonunique: " << args.nonunique_alignments << "\n";
    // out_file << "  Minimum MAPQ: " << args.mapq_min << "\n";
    // out_file << "  Peak Detection: " << args.peak_detection << "\n";
    // out_file << "  Max Components: " << args.max_components << "\n";
    // out_file << "----------------------\n";
    // out_file << "Total Counts: " << total_counts << "\n";
    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";

    // out_file.close();   

    return 0;

}