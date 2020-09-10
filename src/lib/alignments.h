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

		    int jump = pos;
		    
		    // while reads are left on contig
		    while (jump != 0) {

		    	// Jump to desired region in bam
		    	if (!inFile.Jump(ref, jump)) {
		    		break;
		    	}
		    		
		    	std::cerr << ref << "\t" << jump << "\n";

		    	// Initialize data arrays of positions, end positions, and strands
			    int positions[10000] = {0}; 
				int end_positions[10000] = {0};
				char strands[10000];

				// Initialize adjacency matrix
				arma::mat adj_matrix(10000, 10000);
				adj_matrix.zeros(); 

				// Get First alignment
				inFile.GetNextAlignment(alignment);

				// Set values at first index
				end_positions[0] = alignment.GetEndPosition() - 1;
				positions[0] = round((end_positions[0] + alignment.Position) / 2);  
				strands[0] = (alignment.IsReverseStrand()) ? '-' : '+';


				// Initialize loop variables
				int num_alignments = 1;
				int n;
				while (num_alignments < 10000 && inFile.GetNextAlignment(alignment)) {

					// If next chromosome is reached, get out of town.
					if (alignment.RefID > ref) {
						break;
					}

					// Set values at appropriate position
					end_positions[num_alignments] = alignment.GetEndPosition() - 1;
					positions[num_alignments] = round((end_positions[num_alignments] + alignment.Position) / 2); 
					strands[num_alignments] = (alignment.IsReverseStrand()) ? '-' : '+';

					// find overlapping reads preceding this one in bam file
					n = num_alignments - 1;
					while (n >= 0) {

						// if read is too far ahead to overlap, break
						if (alignment.Position > end_positions[n]) {
							break;

						// if reads overlap
						} else if (alignment.Position <= end_positions[n] && strands[num_alignments] == strands[n]) {
							adj_matrix(num_alignments, n) = 1;
							adj_matrix(n, num_alignments) = 1;

						} 

						n--;

					}

					num_alignments++;

				}

				
				// Calculate degree of read node
				arma::rowvec degrees = arma::sum(adj_matrix);
				

				// Find approprite spot to break, we can't cut a group in half can we?
				int total = num_alignments;	

				while (degrees[total - 1] != 0) {

					total--;

					// If entirety of loop is full..... we will come back to this
					if (total < 1) {
						total = num_alignments;
						break;
					}

				}


				// Position to jump to for next loop
				if (total == num_alignments) {
					jump = positions[total - 1];
				} else {
					jump = positions[total];
				}

				std::cerr << total << "\t" << jump << "\n";

				
				// Initialize variables for counting reads by nodes in their graphs
				int i = 0;
				int x;
				int j;
				int increment;

				while(i < total) {

					increment = i + 1;

					// if nodes meet min connectedness requirement
					if (degrees[i] >= min_cov) {

						// initialize group variables. Max and counts.
						// max serves as "peak" to represent graph. It's arbitrary.
						int max = i;
						int nodes = 1;

						// Search upstream of node
						x = i;
						j = i - 1;
						while (j >= 0) {

							if (strands[x] == strands[j]) {

								if (adj_matrix(x, j) == 0) {
									break;
							
								} else {

									nodes++;
									x--;

									if (degrees[j] > degrees[x]) {
										max = j;

									}
								
								}
							}

							j--;

						}

						// Search downstream of node
						x = i;
						j = i + 1;

						while (j < total) {

							if (strands[x] == strands[j]) {

								if (adj_matrix(x, j) == 0) {
									break;
							
								} else {


									nodes++;
									x++;

									if (degrees[j] > degrees[x]) {
										max = j;

									}
								
								}
							}

							j++;

						}

						// increment increment :)
						increment = j + 1;
					
						// report read conts
						std::cout << contig_cache[ref] << "\t" << positions[max] << "\t" <<  nodes << "\n";
						
					}

					i = increment;

				}

			}

		}			


};

