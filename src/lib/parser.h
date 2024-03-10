#include <seqan/arg_parse.h>

using namespace seqan;

/////////////////////////////////////// 
// Arguments Data Structure
struct ImpactArguments {

    // Files
    std::string alignment_file;    	    // sam or bam file
    std::string index_file;    			// index filed
    std::string annotation_file;      	// gff or gtf file

    // Program
    int threads;                        // threads
    std::string library_type;           // library type (SE or PE)
    std::string stranded;           // strandedness
    
    // Alignments
    bool nonunique_alignments;			// count primary and secondary alignments
    int mapq;						// minimum mapq score
    //int min_coverage;					// min coverage

    // Features
    bool isGFF = false;                         // annotaiton file is gff
    std::string feature_tag;            // name of feature tag
    std::string feature_id;             // ID of feature

    // Output
    std::string gtf_output;             // name of output gtf file

};


/////////////////////////////////////// 
// Argument Parser
ArgumentParser::ParseResult argparse(int argc, char const **argv, ImpactArguments &args) {
	
    // Setup ArgumentParser.
    ArgumentParser parser("impact");
    addDescription(parser, 
                   "Identifies expressed transcripts using clusters of mapped reads from TAGseq experiments.Generates a counts file written to stdout and optionally a GTF file of identified read clusters.");


    // Define Arguments
    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "BAM"));

    addArgument(parser, seqan::ArgParseArgument(
        ArgParseArgument::INPUT_FILE, "GTF"));


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
    // addOption(parser, ArgParseOption(
    //     "m", "min-coverage",
    //     "Minimum coverage for target consideration.",
    //     ArgParseArgument::INTEGER, "INT"));
    // setDefaultValue(parser, "min-coverage", "10");

      // Feature Tag
    addOption(parser, seqan::ArgParseOption(
        "f", "feature-tag", "Name of feature tag.",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "feature-tag", "exon");

    addOption(parser, seqan::ArgParseOption(
        "i", "feature-id", "ID of feature (use for GFFs).",
        ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "feature-id", "gene_id");

      // Output GTF
    addOption(parser, seqan::ArgParseOption(
        "o", "output-gtf", "Output read cluster GTF file and specify name.",
        ArgParseArgument::STRING, "STRING"));

    // Add Information 
    addUsageLine(parser, "input.sorted.bam annotation.gtf [options]");
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
   
   // Check file type of first positional arg
    std::string input_file_ext = getFileExtension(getArgument(parser, 0));
    if (input_file_ext != "bam") {
        std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".bam\",  extension.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // Check file type of second positional arg
    input_file_ext = getFileExtension(getArgument(parser, 1));
    if (input_file_ext != "gtf" && input_file_ext != "gff") {
        std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".gtf\" and \".gff\",  extension.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    if (input_file_ext == "gff") {
        args.isGFF = true; 
    }


    // Get arguments
    getArgumentValue(args.alignment_file, parser, 0);
    getArgumentValue(args.annotation_file, parser, 1);
    getArgumentValue(args.index_file, parser, 0);
    args.index_file = args.index_file + ".bai";
    
    // Populate options
    getOptionValue(args.threads, parser, "threads");
    getOptionValue(args.library_type, parser, "library-type");
    getOptionValue(args.stranded, parser, "strandedness");
    args.nonunique_alignments = isSet(parser, "nonunique-alignments");
    getOptionValue(args.mapq, parser, "mapq-min");
    //getOptionValue(args.min_coverage, parser, "min-coverage");
    getOptionValue(args.feature_tag, parser, "feature-tag");
    getOptionValue(args.feature_id, parser, "feature-id");
    getOptionValue(args.gtf_output, parser, "output-gtf");

    return seqan::ArgumentParser::PARSE_OK;

}
