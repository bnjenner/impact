#include <seqan/bam_io.h>
#include <vector>
#include <math.h>
#include "peaks.h"

using namespace seqan;

// Alignment Class
class AlignmentFile {

	public:

	// Attributes
		BamFileIn inFile;				// Bam File Object
    	BamIndex<Bai> baiIndex;			// Bam File Object
    	BamAlignmentRecord record;		// Bam File Object

    	// Parameters
    	bool nonunique_alignments;		// consider secondary alignments 
    	int mapq;						// minimum mapping quality
    	bool peak_detection;			// peak detection

    	// Counts
   		int unmapped = 0;    			// unmapped reads
    	int ambiguous = 0;   			// ambiguous reads
    	int multimapped = 0;			// multimapped reads
    	int low_qual = 0;				// low qualirt reads

    	// File Beginning Indicies
    	int ref_id;						// SOF ref name
    	int begin_pos;					// SOF start position
    	int end_pos;					// SOF end position
		bool hasAlignments; 			// initialize hasAlignment for jumpToRegion()


    // Inialize
    	AlignmentFile(CharString bam_file, CharString index_file, int pre_mapq, 
    				  bool pre_nonunique_alignments, bool pre_peak_detection) {

    		// Set Attributes
    		mapq = pre_mapq;
    		nonunique_alignments = pre_nonunique_alignments;
    		peak_detection = pre_peak_detection;

    		// Test Bam File Opening
			if (!open(inFile, toCString(bam_file))) {
		        std::cerr << "ERROR: Could not read alignment file: " << bam_file << "\n";
		        throw "ERROR: Could not read alignment file.";
		    }

		    // Test Bam Index File Opening
		    if (!open(baiIndex, toCString(index_file))) {
		        std::cerr << "ERROR: Could not read index file: " << index_file << "\n";
		        throw "ERROR: Could not read index file";
		    }

		    // Attach Bam file to standard out
		    BamFileOut out(context(inFile), std::cout, Sam());

		    // Read in Header, a necessary evil
		    BamHeader header;
		    readHeader(header, inFile);

		   	bool first = true;

		   	while (!atEnd(inFile)) {

		   		readRecord(record, inFile);

		   		// If first entry, use these stats for beginning of file
		   		if (first) {
		   			ref_id = getIdByName(ref_id, contigNamesCache(context(inFile)), record.rID);
				   	begin_pos = record.beginPos;
				   	end_pos = record.beginPos + 1;
				   	first = false;
				}


		        if (hasFlagUnmapped(record)) {
		        	unmapped++;
		        	continue;
		        } 

		        // PSEUDO CODE 
		        // if (record.qual <= mapq && countedID(record.rID) == False) {
		        // 	low_qual++;
		        // 	continue;
		        // } 

		        // HOW TO ADDRESS MULTIMAPPING AND OVERLAPS?

		        // PSEUDO CODE
		        // if (hasFlagMultiple(record) && nonunique-alignments && countedID(record.rID) == False){
		        // 	multimapped++;
		        // 	continue;
		        // }
		    }  

		    // std::cout << "Unmapped: " << unmapped << std::endl;
		    // std::cout << "Multimapped: " << multimapped << std::endl;

		   	// Jump to the first entry in file.
		    jumpToRegion(inFile, hasAlignments, ref_id, begin_pos, end_pos, baiIndex);

		}

    // Methods
    	// head Bam
		void head() {

			BamFileOut out(context(inFile), std::cout, Sam());
			BamAlignmentRecord record;

			int i = 0;
		 	while (i < 10) { // (!atEnd(bamFileIn))

		        readRecord(record, inFile);
		        writeRecord(out, record);
		        i++;

		    }  

		    // Jump the first entry in file.
		    jumpToRegion(inFile, hasAlignments, ref_id, begin_pos, end_pos, baiIndex);
		}


		// Grab Alignments within Interval Using Bam Index
		int findAlignments(CharString feature_name, CharString ref, int beginPos, int endPos, char strand) {

		    // 1-based to 0-based.
		    beginPos -= 1;
		    endPos -= 1;

		    // Translate from reference name to rID.
		    int rID = 0;
		    if (!getIdByName(rID, contigNamesCache(context(inFile)), ref)) {
		        std::cerr << "ERROR: Reference sequence named " << ref << " not known.\n";
		        return 0;
		    }


		    jumpToRegion(inFile, hasAlignments, rID, beginPos, endPos, baiIndex);

		    // Check if alignments exits
 		    if (!hasAlignments) {
		        return 0;  
		    }

		   
		    MappingCounts mappedCounts(feature_name, beginPos, endPos);


		    // Seek linearly to the selected position.
		    BamAlignmentRecord record;
		    BamFileOut out(context(inFile), std::cout, Sam());

		    int num_alignments = 0;
		    int align_end;

		    while (!atEnd(inFile)) {

		        readRecord(record, inFile);

		        // Break if not mapped, past reference id, or end of file.
		        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos) {
		            break;
		        }

		        // Check Strandedness
		        if ((hasFlagRC(record) && strand != '-') || (!hasFlagRC(record) && strand != '+')) {
		        	continue;
		        } 

		        // Check if primary alignment
		        if (hasFlagSecondary(record) && (nonunique_alignments == false)) {
		        	continue;
		        }

		        // Check if sufficient quality
		        // if (record.mapQ && mapq != -1) {
		        // 	continue;
		        // }

		        align_end = record.beginPos - 1 + getAlignmentLengthInRef(record);
		        
		        if (peak_detection) {
		        	mappedCounts.addRead(record.beginPos, align_end);
		        }

		        num_alignments++;
	
		    }


		    if (feature_name == "gl1315.NS.00574" && peak_detection) {
		   		 mappedCounts.write(); 
		    }


		    // Return to beginnig 
		    jumpToRegion(inFile, hasAlignments, ref_id, begin_pos, end_pos, baiIndex);

		   	return num_alignments;
		}

};

