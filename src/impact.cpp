#include <iostream>
#include <chrono> 
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include "utils/parser.h"
#include "utils/peaks.h"
#include "utils/alignments.h"
#include "utils/annotations.h"
#include "utils/count.h"


// Main 
int main(int argc, char const ** argv) {

	auto start = std::chrono::high_resolution_clock::now(); 

	ImpactArguments args;
	seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);


	// Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;AlignmentFile alignment(args);
    }

    
    std::cerr << "[ IMPACT ]\n";

    // Construct alignment object
    std::cerr << "[ Parsing Alignment File... ]\n";
    AlignmentFile alignment(args);


    // Construct annotation object
    std::cerr << "[ Parsing Annotation File... ]\n";
    AnnotationFile annotation(args);


    // Count Reads
    std::cerr << "[ Counting Reads... ]\n";
   	int total_counts = getCounts(annotation, alignment, args.peak_detection);
   	

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
    std::cerr << "[ Program Complete! ]\n";
    std::cerr << "[ Runtime: " << duration.count() << " seconds ]\n";

    // out_file.close();   

    return 0;

}
