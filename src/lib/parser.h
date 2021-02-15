#include <seqan/arg_parse.h>

using namespace seqan;

/////////////////////////////////////// 
// Arguments Data Structure
struct ImpactArguments {

    std::string alignment_file;    	    // sam or bam file
    std::string index_file;    			// index file
    std::string gff_file;      			// gff file
    int threads;                        // threads
    std::string library_type;           // library type (SE or PE)
    std::string strandedness;           // strandedness
    bool nonunique_alignments;			// count primary and secondary alignments
    int mapq_min;						// minimum mapq score
    int min_coverage;					// min coverage

};

// parameters definition for alignment and graph classes (maybea a little redundant)
struct Parameters {

    char library_type;              // library type (p, s)
    char stranded;                  // strandedness of library (f, r)
    bool nonunique_alignments;      // consider secondary alignments 
    int mapq;                       // minimum mapping quality
    int min_cov;                    // min coverage for cluster detection

};


/////////////////////////////////////// 
// Argument Parser
ArgumentParser::ParseResult argparse(int argc, char const **argv, ImpactArguments &args) {
	
    // Setup ArgumentParser.
    ArgumentParser parser("impact");
    addDescription(parser, "Generates read counts and identifies peaks in mapped reads from TAGseq experiments.");


	// Define Arguments
    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "BAM"));

    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "GFF"));


    // Define Options
      // Threads
    addOption(parser, ArgParseOption(
        "t", "threads",
        "Number of processes for multithreading.",
        ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", "1");

      // Library Type
    addOption(parser, seqan::ArgParseOption(
        "l", "library-type", "Library type. Paired end is not recommended. Only used to check proper pairing.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "library-type", "single");
    setValidValues(parser, "library-type", "single paired");

      // Strandedness
    addOption(parser, seqan::ArgParseOption(
        "s", "strandedness", "Strandedness of library.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "strandedness", "forward");
    setValidValues(parser, "strandedness", "forward reverse");

      // Nonunique Alignments
	addOption(parser, seqan::ArgParseOption(
        "n", "nonunique-alignments", "Count primary and secondary read alignments."));

	  // Min Quality
	addOption(parser, ArgParseOption(
	    "q", "mapq-min",
	    "Minimum mapping quality score to consider for counts.",
	    ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "mapq-min", "1");

	  // Min coverage for peak 
	addOption(parser, ArgParseOption(
	    "m", "min-coverage",
	    "Minimum coverage for target consideration.",
	    ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "min-coverage", "5");


	// Add Information 
	addUsageLine(parser, "input.sorted.bam input.gff [options]");
    setDefaultValue(parser, "version-check", "OFF");
	hideOption(parser, "version-check");
	setVersion(parser, "dev0");
    setDate(parser, "July 2020");


    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);


    // Check if Parse was successful
    if (res != ArgumentParser::PARSE_OK) {
        return ArgumentParser::PARSE_ERROR;
    }

    
    // Arguments (deduce file type)
    std::string file_exts[2] = {getFileExtension(getArgument(parser, 0)),
    							getFileExtension(getArgument(parser, 1))};
    
    
    // Identify input files types and index
    for (int i = 0; i < 2; i++) {

        // Assign alignment files
	    if (file_exts[i] == "bam") {

	    	getArgumentValue(args.alignment_file, parser, i);
	    	getArgumentValue(args.index_file, parser, i);
            args.index_file = args.index_file + ".bai";

        // Assign annotation files
	    } else if (file_exts[i] == "gff" || file_exts[i] == "gtf") {

	    	getArgumentValue(args.gff_file, parser, i);	

        // Throw error if unrecognized file type
	    } else {

	    	std::cerr << "ERROR: Unaccapetd File Format: \"." << file_exts[i] <<  "\". Accepts \".bam\", \".gtf\", or \"gff\" extension.\n";
	    	return ArgumentParser::PARSE_ERROR;
	    }
	}    
    

	// Populate options
    getOptionValue(args.threads, parser, "threads");
    getOptionValue(args.library_type, parser, "library-type");
    getOptionValue(args.strandedness, parser, "strandedness");
	args.nonunique_alignments = isSet(parser, "nonunique-alignments");
   	getOptionValue(args.mapq_min, parser, "mapq-min");
    getOptionValue(args.min_coverage, parser, "min-coverage");

    return seqan::ArgumentParser::PARSE_OK;
}
