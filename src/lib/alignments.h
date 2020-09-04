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
    	int min_cov;					// min coverage for cluster detection

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
    		min_cov = args -> min_coverage - 1;

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

			end_positions[0] = alignment.GetEndPosition() - 1;
			positions[0] = round((end_positions[0] + alignment.Position) / 2);  
			strands[0] = (alignment.IsReverseStrand()) ? '-' : '+';

			int num_alignments = 1;
			int n;

			while (num_alignments < 1000) {

				inFile.GetNextAlignment(alignment);

				n = num_alignments - 1;
				end_positions[num_alignments] = alignment.GetEndPosition() - 1;
				positions[num_alignments] = round((end_positions[num_alignments] + alignment.Position) / 2); 
				strands[num_alignments] = (alignment.IsReverseStrand()) ? '-' : '+';

				while (n >= 0) {

					if (alignment.Position > end_positions[n]) {
						break;

					} else if (alignment.Position <= end_positions[n] && strands[num_alignments] == strands[n]) {
						adj_matrix(num_alignments, n) = 1;
						adj_matrix(n, num_alignments) = 1;

					} 

					n--;

				}

				num_alignments++;

			}

			arma::rowvec degrees = arma::sum(adj_matrix);
			int i = 0;
			int j;
			int increment;

			while(i < 1000) {

				increment = i + 1;

				if (degrees[i] >= min_cov) {

					int max = i;
					int counts = degrees[i];
					int nodes = 1;

					// Reverse search
					j = i - 1;
					while (j >= 0) {

						if (strands[i] == strands[j]) {

							if (adj_matrix(i, j) == 0) {
								break;
						
							} else {

								nodes++;

								if (degrees[j] > degrees[i]) {
									max = j;

								}
							
							}
						}

						j--;

					}

					// Forward search
					j = i + 1;
					while (j < 1000) {

						if (strands[i] == strands[j]) {

							if (adj_matrix(i, j) == 0) {
								break;
						
							} else {

								nodes++;

								if (degrees[j] > degrees[i]) {
									max = j;

								}
							
							}
						}

						j++;

					}

					increment = j;
				
					std::cerr << contig_cache[ref] << "\t" << positions[max] << "\t" <<  nodes << "\n";
					
				}

				i = increment;

			}

			// int max = arma::index_max(degrees);

			// int sum = degrees[max];
			// int counter = 1;

			// arma::ivec indicies = arma::conv_to<arma::ivec>::from(find(adj_matrix.col(max)));
			// arma::ivec temp;
			// arma::ivec overlap;

			// for (int i = 0; i < indicies.size(); i++) {

			// 	temp = arma::conv_to<arma::ivec>::from(find(adj_matrix.col(indicies[i])));
			// 	overlap = intersect(indicies, temp);

			// 	if (indicies.size() == overlap.size() - 1) {
			// 		break;
			// 	}



			// }


			// std::cerr << max << "\n";
			// std::cerr << positions[max] << "\n";



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

