#include <seqan/arg_parse.h>

///////////////////////////////////////
// Arguments Data Structure
struct ImpactArguments {

    // Files
    std::string alignment_file;         // sam or bam file
    std::string index_file;             // index filed
    std::string annotation_file;        // gff or gtf file

    // Program
    int threads;                        // threads
    std::string library_type;           // library type (SE or PE)
    std::string stranded;               // strandedness

    // Alignments
    bool nonunique_alignments;          // count primary and secondary alignments
    int mapq;                           // minimum mapq score
    //int min_coverage;                 // min coverage

    // Features
    bool isGFF = false;                 // annotaiton file is gff
    std::string feature_tag;            // name of feature tag
    std::string feature_id;             // ID of feature

    // Output
    std::string gtf_output;             // name of output gtf file
};


// Argument Parser
seqan::ArgumentParser::ParseResult argparse(int argc, char const **argv, ImpactArguments &args) {

    // Setup ArgumentParser.
    seqan::ArgumentParser parser("impact");
    seqan::addDescription(parser,
                          "Identifies expressed transcripts using clusters of mapped reads from TAGseq experiments.Generates a counts file written to stdout and optionally a GTF file of identified read clusters.");


    // Define Arguments
    seqan::addArgument(parser, seqan::ArgParseArgument(
                           seqan::ArgParseArgument::INPUT_FILE, "BAM"));
    seqan::addArgument(parser, seqan::ArgParseArgument(
                           seqan::ArgParseArgument::INPUT_FILE, "GTF"));


    // Define Options
    seqan::addOption(parser, seqan::ArgParseOption(
                         "t", "threads",
                         "Number of processes for multithreading.",
                         seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "threads", "1");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "l", "library-type", "Library type. Paired end is not recommended. Only used to check proper pairing.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "library-type", "single");
    seqan::setValidValues(parser, "library-type", "single paired");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "s", "strandedness", "Strandedness of library.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "strandedness", "forward");
    seqan::setValidValues(parser, "strandedness", "forward reverse");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "n", "nonunique-alignments", "Count primary and secondary read alignments."));

    seqan::addOption(parser, seqan::ArgParseOption(
                         "q", "mapq-min",
                         "Minimum mapping quality score to consider for counts.",
                         seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "mapq-min", "1");

    // addOption(parser, ArgParseOption(
    //     "m", "min-coverage",
    //     "Minimum coverage for target consideration.",
    //     ArgParseArgument::INTEGER, "INT"));
    // setDefaultValue(parser, "min-coverage", "10");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "f", "feature-tag", "Name of feature tag.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "feature-tag", "exon");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "i", "feature-id", "ID of feature (use for GFFs).",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "feature-id", "gene_id");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "o", "output-gtf", "Output read cluster GTF file and specify name.",
                         seqan::ArgParseArgument::STRING, "STRING"));

    seqan::addUsageLine(parser, "input.sorted.bam annotation.gtf [options]");
    seqan::setDefaultValue(parser, "version-check", "OFF");
    seqan::hideOption(parser, "version-check");
    seqan::setVersion(parser, "dev0");
    seqan::setDate(parser, "July 2020");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if Parse was successful
    if (res != seqan::ArgumentParser::PARSE_OK) { return seqan::ArgumentParser::PARSE_ERROR; }

    // Check file type of first positional arg
    std::string input_file_ext = seqan::getFileExtension(getArgument(parser, 0));
    if (input_file_ext != "bam") {
        std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".bam\",  extension.\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Check file type of second positional arg
    input_file_ext = seqan::getFileExtension(getArgument(parser, 1));
    if (input_file_ext != "gtf" && input_file_ext != "gff") {
        std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".gtf\" and \".gff\",  extension.\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    if (input_file_ext == "gff") { args.isGFF = true; }

    // Get arguments
    seqan::getArgumentValue(args.alignment_file, parser, 0);
    seqan::getArgumentValue(args.annotation_file, parser, 1);
    seqan::getArgumentValue(args.index_file, parser, 0);
    args.index_file = args.index_file + ".bai";

    // Populate options
    seqan::getOptionValue(args.threads, parser, "threads");
    seqan::getOptionValue(args.library_type, parser, "library-type");
    seqan::getOptionValue(args.stranded, parser, "strandedness");
    args.nonunique_alignments = seqan::isSet(parser, "nonunique-alignments");
    seqan::getOptionValue(args.mapq, parser, "mapq-min");
    //getOptionValue(args.min_coverage, parser, "min-coverage");
    seqan::getOptionValue(args.feature_tag, parser, "feature-tag");
    seqan::getOptionValue(args.feature_id, parser, "feature-id");
    seqan::getOptionValue(args.gtf_output, parser, "output-gtf");

    return seqan::ArgumentParser::PARSE_OK;
}
