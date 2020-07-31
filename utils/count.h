#include <seqan/gff_io.h>

using namespace seqan;

// Iterate Through GFF
void getCounts(CharString gff_file, AlignmentFile &alignment, bool peak_detection) {
	
	// Read in GFF
    GffFileIn gffIn;
    if (!open(gffIn, toCString(gff_file))) {
        std::cerr << "ERROR: Could not read gff file: " << gff_file << "\n";
		throw "ERROR: Could not read gff file.";
    }

    // Attach to standard output.
    GffFileOut gffOut(std::cout, Gff());

    // Copy the file record by record.
    GffRecord record;

    CharString feature_name = "";
    int num_alignments = 0;
    int i = 0;
    int start;
    int stop;

    while (!atEnd(gffIn)) {

        readRecord(record, gffIn);

        if (feature_name != record.tagValues[0]) {

        	if (feature_name != "") {

        		MappingCounts mappedCounts(feature_name, start, stop);

        		num_alignments += alignment.findAlignments(mappedCounts, record.ref, start, 
       													   stop, record.strand);


        		std::cout << feature_name << ": " << num_alignments << std::endl; 

        	}

        	num_alignments = 0;
        	feature_name = record.tagValues[0];
        	start = record.beginPos;
        	stop = record.endPos;

        } else {

        	stop = record.endPos;
        	continue;

        }

        i++;

    }

    std::cout << feature_name << ": " << num_alignments << std::endl; 
}
