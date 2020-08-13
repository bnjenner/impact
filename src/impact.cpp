#include <iostream>
#include <chrono> 
#include <string>
#include <seqan/arg_parse.h>
#include "utils/peaks.h"
#include "utils/alignments.h"
#include "utils/count.h"

using namespace seqan;

// Argument Data Structure
struct ImpactArguments
{ 
    std::string alignment_file;    	// sam or bam file
    std::string index_file;    			// index file
    std::string gff_file;      			// gff file
    std::string strandedness;           // library strandedness
    std::string library_type;           // library type (SE or PE)
    bool nonunique_alignments;			// count primary and secondary alignments
    int mapq_min;						// minimum mapq score
    bool peak_detection; 				// enable peak detection 
    int max_components;					// max components for GMM

    CharString index_suffix = ".bai";	// Suffix for index file
};


// Argument Parser
ArgumentParser::ParseResult argparse(int argc, char const ** argv, ImpactArguments & args) {
	// Setup ArgumentParser.
    ArgumentParser parser("impact");
    addDescription(parser, "Generates read counts and identifies peaks in mapped reads from TAGseq experiments.");

	// Define Options
    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "BAM"));

    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "GFF"));

      // Strandedness 
    addOption(parser, seqan::ArgParseOption(
        "s", "strandedness", "Strandedness of library.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "strandedness", "forward");
    setValidValues(parser, "strandedness", "forward reverse unstranded");

      // Library Type
    addOption(parser, seqan::ArgParseOption(
        "l", "library-type", "Library type.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "library-type", "paired");
    setValidValues(parser, "library-type", "single paired");

      // Nonunique Alignments
	addOption(parser, seqan::ArgParseOption(
        "n", "nonunique-alignments", "Count primary and secondary read alignments."));

	  // Min Quality
	addOption(parser, ArgParseOption(
	    "q", "mapq-min",
	    "Minimum mapping quality score to consider for counts.",
	    ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "mapq-min", "-1");

	  // Find 
	addOption(parser, seqan::ArgParseOption(
        "p", "peak-detection", "Use peak detection to identify terminal exon variants."));

	  // Number of Components for GMM
	addOption(parser, ArgParseOption(
	    "m", "max-components",
	    "Maximum number of components for Gaussian Mixture Model.",
	    ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "max-components", "4");

	// Add Information 
	addUsageLine(parser, "input.bam input.gff [options]");
    setDefaultValue(parser, "version-check", "OFF");
	hideOption(parser, "version-check");
	setVersion(parser, "dev0");
    setDate(parser, "July 2020");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Check if Parse was successful
    if (res != ArgumentParser::PARSE_OK)

        return ArgumentParser::PARSE_ERROR;

    // Populate ImpactArguments 

      // Arguments (deduce file type)
    std::string file_exts[2] = {getFileExtension(getArgument(parser, 0)),
    							getFileExtension(getArgument(parser, 1))};
    
    for (int i = 0; i < 2; i++) {

	    if (file_exts[i]== "bam") {

	    	getArgumentValue(args.alignment_file, parser, i);
	    	getArgumentValue(args.index_file, parser, i);	

	    } else if (file_exts[i] == "gff" || file_exts[i] == "gtf") {

	    	getArgumentValue(args.gff_file, parser, i);	

	    } else {

	    	std::cerr << "ERROR: Unaccapetd File Format: \"." << file_exts[i] <<  "\". Accepts \".bam\", \".gtf\", or \"gff\" extension.\n";
	    	return ArgumentParser::PARSE_ERROR;
	    }
	}    
    
    // Add Index Suffix
    for (int i = 0; i < length(args.index_suffix); i++) {
		appendValue(args.index_file, args.index_suffix[i]);
	}

	  // Options
    getOptionValue(args.strandedness, parser, "strandedness");
    getOptionValue(args.library_type, parser, "library-type");
	args.nonunique_alignments = isSet(parser, "nonunique-alignments");
   	getOptionValue(args.mapq_min, parser, "mapq-min");
   	args.peak_detection = isSet(parser, "peak-detection");
   	getOptionValue(args.max_components, parser, "max-components");

    return seqan::ArgumentParser::PARSE_OK;
}


// Main 
int main(int argc, char const ** argv) {

	auto start = std::chrono::high_resolution_clock::now(); 

	ImpactArguments args;
	ArgumentParser::ParseResult res = argparse(argc, argv, args);

	// Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;


    // std::cerr << "----------------------\n";
    // std::cerr << "Parsing Alignment File\n";

    // Construct alignment object
    AlignmentFile alignment(args.alignment_file, args.index_file, 
                            args.strandedness, args.library_type,
    						args.mapq_min, args.nonunique_alignments, 
    						args.peak_detection, args.max_components);


    // std::cerr << "--------------\n";
    // std::cerr << "Counting Reads\n";

    // Count Reads
   	int total_counts = getCounts(args.gff_file, alignment, args.peak_detection);
   	

    // Close alignment file
    alignment.close();
   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

   
	// std::cerr << "-----------------\nProcess Complete!\n";
	// std::cerr << "Runtime: " << duration.count() << " seconds.\n";


    std::ofstream out_file;
    out_file.open("impact_stats.txt");

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
    out_file << "Total Counts: " << total_counts << "\n";
    out_file << "Runtime: " << duration.count() << " seconds.\n";

    out_file.close();   

    return 0;

}
