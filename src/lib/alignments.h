#include "api/BamReader.h"

using namespace BamTools;

// Alignment Class
class AlignmentFile {

	public:

	// Attributes
		std::string file_name;			// Bam File Name 
		std::string index_file;			// Bam index File Name
		BamReader inFile;				// Bam File Objec
		BamAlignment alignment;			// BamAlignmentRecord record;		

    	// Program options 
    	std::string library_type;		// library type
    	bool nonunique_alignments;		// consider secondary alignments 
    	int mapq;						// minimum mapping quality

    	// Counts
   		int unmapped = 0;    			// unmapped reads
    	int ambiguous = 0;   			// ambiguous reads
    	int multimapped = 0;			// multimapped reads
    	int low_qual = 0;				// low qualirt reads

    	std::unordered_map<std::string,int> contig_cache; 	// Unordered map 


    // Constructors
    	AlignmentFile() {}

    	AlignmentFile(const ImpactArguments *args) {

    		// Set Attributes
    		file_name = args -> alignment_file;
    		index_file = args -> index_file;
    		library_type = args -> library_type;
    		nonunique_alignments = args -> nonunique_alignments;
    		mapq = args -> mapq_min - 1;

    	}


    // Methods

    	void open() {

			if (!inFile.Open(file_name)) {
			    std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
			    throw "ERROR: Could not read alignment file.";
			}

			if (!inFile.OpenIndex(index_file)) {
			    std::cerr << "ERROR: Could not read index file: " << index_file << "\n";
			    throw "ERROR: Could not read index file";
			}

			// if (!inFile.SetRegion(0, 2311417, 0, 2313233)) {
			// 	std::cerr << "ERROR: Unable to jump to region.\n";
			// }

			SamHeader head = inFile.GetHeader();
			if (head.HasSortOrder()) {

				std::string sortOrder = head.SortOrder;

				if (sortOrder.compare("coordinate") != 0) {
					std::cerr << "ERROR: Sorted alignment file required.\n";
					throw "ERROR: Could not read alignment file.";
				}

			} else {
				std::cerr << "ERROR: BAM file has no @HD SO:<SortOrder> attribute. Impossible to determine sort order.\n";
					throw "ERROR: Could determine sort status.";

			}

			// Generate Ref Map (contig indicies)
			RefVector references = inFile.GetReferenceData();
			for (int i = 0; i < references.size(); i++) {

				contig_cache[references.at(i).RefName] = i;

			}


			while (inFile.GetNextAlignment(alignment)) {

				 if (!alignment.IsMapped()) {
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

			// std::cerr << "Unmapped: " << unmapped << std::endl;
			// std::cerr << "Multimapped: " << multimapped << std::endl;

			// Jump to the first entry in file
			inFile.Rewind(); 	

		}

		// close file
		void close() {

			inFile.Close();

		}

		
	   	// head Bam
		void head() {

			int i = 0;
		 	while (inFile.GetNextAlignment(alignment) && i < 10) { // (!atEnd(bamFileIn))

		        std::cout << alignment.Name << "\n";
		        i++;

		    }  

		   	// Jump to the first entry in file.
		   	inFile.Rewind(); 	

		}


		// Grab Alignments within Interval Using Bam Index
		int findAlignments(MappingCounts &mappedCounts, std::string ref, int beginPos, int endPos, char strand,
						   bool peak_detection, int prev_stop, char prev_strand, int next_start, char next_strand) {

		    int rID = contig_cache[ref];

 			if (!inFile.SetRegion(rID, beginPos, rID, endPos)) {
 				std::cerr << mappedCounts.feature_name << " " << ref << " " << rID << beginPos << "\n";
		    	return 0;
		    }

		    int num_alignments = 0;
		    int align_end;


			while (inFile.GetNextAlignment(alignment)) {

				// Check if position in alignment file is beyond desired range
				if (alignment.RefID != rID || alignment.Position > endPos) {
		            break;
		        }

		        // Check if alignment is on the same strand
		        if ((alignment.IsReverseStrand() && strand != '-') || (!alignment.IsReverseStrand() && strand == '-')) {
					continue;
				} 

				// Check if primary alignment
				if (!alignment.IsPrimaryAlignment() && (nonunique_alignments == false)) {
					continue;
				}

		        // Check if sufficient quality
				if (alignment.MapQuality <= mapq) {
		       		continue;
				}

				if (strand == '+') {

					// Check if alignment doens't start within transcript
					if (alignment.Position < beginPos) {
						continue;
					}


				} else {

					// Check if alignment doens't start within transcript
					if (alignment.GetEndPosition() > endPos) {
						continue;
					}

				}

				if (next_start != -1 && strand == next_strand) {

					if (alignment.GetEndPosition() >= next_start) {
						continue;

					} 

				}

	        	if (prev_stop != -1 && strand == prev_strand) {

	        		if (alignment.Position <= prev_stop) {
						continue;

					} 

	        	}

		        if (peak_detection) {
		        	mappedCounts.addRead(alignment.Position, alignment.GetEndPosition());
		        }

		        num_alignments++;
	
		    }


		    // if (peak_detection && num_alignments > 50) {
		   	// 	 mappedCounts.fit(max_components); 
		    // }

		    // Return to beginnig 
		    inFile.Rewind();

		   	return num_alignments;
		}

};

