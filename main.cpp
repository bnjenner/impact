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
    addDescription(parser, "Generates read counts from TAGseq experiments and identifies terminal exon isoforms.");

	// Define Options
    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "BAM"));

    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "GFF"));

      // Nonunique Alignments
	addOption(parser, seqan::ArgParseOption(
        "n", "nonunique-alignments", "Count primary and secondary read alignments."));

	  // Min Quality
	addOption(parser, ArgParseOption(
	    "q", "mapq-min",
	    "Minimum mapping quality score to consider for counts (default is no minimum).",
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


    std::cerr << "----------\nIMPACT\n----------\n";
    std::cerr << "Parameters:\n";
    std::cerr << "  Alignment File: " << args.alignment_file << "\n";
    std::cerr << "  Index File: " << args.index_file << "\n";
    std::cerr << "  GFF File: " << args.gff_file << "\n";
    std::cerr << "  Nonunique: " << args.nonunique_alignments << "\n";
    std::cerr << "  Minimum MAPQ: " << args.mapq_min << "\n";
    std::cerr << "  Peak Detection: " << args.peak_detection << "\n";
    std::cerr << "  Max Components: " << args.max_components << "\n";
    

    std::cerr << "----------------------\n";
    std::cerr << "Parsing Alignment FIle\n";

    AlignmentFile alignment(args.alignment_file, args.index_file, 
    						args.mapq_min, args.nonunique_alignments, 
    						args.peak_detection, args.max_components);

    std::cerr << "--------------\n";
    std::cerr << "Counting Reads\n";

   	//getCounts(args.gff_file, alignment, args.peak_detection);
   	
   	
   	auto stop = std::chrono::high_resolution_clock::now(); 
   	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

   
	std::cerr << "-----------------\nProcess Complete!\n";
	std::cerr << "Runtime: " << duration.count() << " seconds.\n";

    return 0;

}
