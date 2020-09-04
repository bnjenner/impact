#include "api/BamReader.h"
#include <armadillo> 

using namespace BamTools;


//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

		// Attributes
		BamReader inFile;				// Bam File Object
		BamAlignment alignment;			// BamAlignmentRecord record;	

    	// Program options 
    	std::string file_name;
    	std::string index;
    	std::string library_type;		// library type
    	std::string stranded;			// strandedness of library
    	bool nonunique_alignments;		// consider secondary alignments 
    	int mapq;						// minimum mapping quality

    	// Count Statistics
    	int noncounts[5] = {0,0,0,0,0};
    		// Order:
    		  // 0 = Features
    		  // 1 = No Feature
    		  // 2 = Ambiguous
    		  // 3 = Low Qual
    		  // 4 = Unmapped
   		
    	// Contig Cache
    	std::unordered_map<int, std::string> contig_cache; 	// Unordered map 


    // Constructors
    	AlignmentFile() {}

    	AlignmentFile(const ImpactArguments *args) {

    		// Set Attributes
    		file_name = args -> alignment_file;
    		index = args -> index_file;
    		library_type = args -> library_type;
    		stranded = args -> strandedness;
    		nonunique_alignments = args -> nonunique_alignments;
    		mapq = args -> mapq_min - 1;

    	}


    // Methods
    	void open() {

			if (!inFile.Open(file_name)) {
			    std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
			    throw "ERROR: Could not read alignment file.";
			}
			
			if (!inFile.OpenIndex(index)) {
				std::cerr << "ERROR: Could not read index file: " << index << "\n";
				throw "ERROR: Could not read index file";
		    }


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

				contig_cache[i] = references.at(i).RefName;

			}
	

		}


		// close file
		void close() {

			inFile.Close();

		}

		
	   	// head Bam
		void head() {

			int i = 0;

		 	while (inFile.GetNextAlignment(alignment) && i < 10) { // (!atEnd(bamFileIn))

		        std::cout << alignment.Name << "\t" << alignment.RefID << "\t"  << alignment.IsReverseStrand() << "\t"
		        		  << alignment.Position << "\t" << alignment.GetEndPosition() << "\n";

		        i++;

		    }  

		   	// Jump to the first entry in file.
		   	inFile.Rewind(); 	

		}


	// Grab Alignments within Interval Using Bam Index
		void get_clusters(int ref, int pos) {

		    // 1-based to 0-based.

 			if (!inFile.Jump(ref, pos)) {
		    	return;
		    }

		    int positions[1000] = {0}; 
			int end_positions[1000] = {0};
			char strands[1000];
			arma::mat adj_matrix(1000, 1000); 

			inFile.GetNextAlignment(alignment);

			end_positions[0] = alignment.GetEndPosition();
			positions[0] = round((end_positions[0] + alignment.Position) / 2);  
			strands[0] = (alignment.IsReverseStrand()) ? '-' : '+';

			int num_alignments = 1;
			int i;

			while (num_alignments < 1000) {

				inFile.GetNextAlignment(alignment);

				i = num_alignments - 1;
				end_positions[num_alignments] = alignment.GetEndPosition();
				positions[num_alignments] = round((end_positions[num_alignments] + alignment.Position) / 2); 
				strands[num_alignments] = (alignment.IsReverseStrand()) ? '-' : '+';

				while (i >= 0) {

					//std::cerr << alignment.Position << "\t" << end_positions[i] << "\n";
					if (alignment.Position > end_positions[i]) {
						break;

					} else if (alignment.Position < end_positions[i] && strands[num_alignments] == strands[i]) {
						adj_matrix(num_alignments, i) = 1;
						adj_matrix(i, num_alignments) = 1;
						i--; 

					}

					i--;

				}

				num_alignments++;

			}
			
			arma::rowvec degrees = arma::sum(adj_matrix);
			int max = arma::index_max(degrees);
			std::cerr << max << "\n";
			std::cerr << positions[max] << "\n";



		//         // Break if not mapped, past reference id, or end of file.
		//         if (alignment.RefID > rID || alignment.Position >= endPos) {
		//             break;
		//         }

		//         // Check Strandedness
		//         if (strandedness != "unstranded") {
		// 	        if ((alignment.IsReverseStrand() && strand != '-') || (!alignment.IsReverseStrand() && strand == '-')) {
		// 	        	continue;
		// 	        }
		//     	}
		     
		//         // Check if primary alignment
		//         if (!alignment.IsPrimaryAlignment() && (nonunique_alignments == false)) {
		//         	continue;
		//         }

		//         // Check if sufficient quality
		// 		if (alignment.MapQuality <= mapq) {
		//        		continue;
		// 		}
		        
		//         num_alignments++;
	
		//     }

		//     // Return to beginnig 
			inFile.Rewind();

		   	//return num_alignments;
		}			


};

