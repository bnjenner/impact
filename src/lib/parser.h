#include <seqan/arg_parse.h>

using namespace seqan;

/////////////////////////////////////// 
// Arguments Data Structure

struct ImpactArguments {

    std::string alignment_file;    	// sam or bam file
    std::string index_file;    			// index file
    std::string gff_file;      			// gff file
    std::string library_type;           // library type (SE or PE)
    bool nonunique_alignments;			// count primary and secondary alignments
    int mapq_min;						// minimum mapq score
    bool peak_detection; 				// enable peak detection 
    int max_components;					// max components for GMM

};


/////////////////////////////////////// 
// Argument Parser

ArgumentParser::ParseResult argparse(int argc, char const **argv, ImpactArguments &args) {
	// Setup ArgumentParser.
    ArgumentParser parser("impact");
    addDescription(parser, "Generates read counts and identifies peaks in mapped reads from TAGseq experiments.");

	// Define Options
    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "BAM"));

    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "GFF"));


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
	setDefaultValue(parser, "mapq-min", "1");

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
    if (res != ArgumentParser::PARSE_OK) {
        return ArgumentParser::PARSE_ERROR;
    }


    // Populate ImpactArguments 

    
    // Arguments (deduce file type)
    std::string file_exts[2] = {getFileExtension(getArgument(parser, 0)),
    							getFileExtension(getArgument(parser, 1))};
    
    
    for (int i = 0; i < 2; i++) {

	    if (file_exts[i] == "bam") {

	    	getArgumentValue(args.alignment_file, parser, i);
	    	getArgumentValue(args.index_file, parser, i);	

	    } else if (file_exts[i] == "gff" || file_exts[i] == "gtf") {

	    	getArgumentValue(args.gff_file, parser, i);	

	    } else {

	    	std::cerr << "ERROR: Unaccapetd File Format: \"." << file_exts[i] <<  "\". Accepts \".bam\", \".gtf\", or \"gff\" extension.\n";
	    	return ArgumentParser::PARSE_ERROR;
	    }
	}    
    

	// Options
    getOptionValue(args.library_type, parser, "library-type");
	args.nonunique_alignments = isSet(parser, "nonunique-alignments");
   	getOptionValue(args.mapq_min, parser, "mapq-min");
   	args.peak_detection = isSet(parser, "peak-detection");
   	getOptionValue(args.max_components, parser, "max-components");

    return seqan::ArgumentParser::PARSE_OK;
}
