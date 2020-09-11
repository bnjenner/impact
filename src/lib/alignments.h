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
		void get_counts(int ref, int pos) {

		    int jump = pos;
		    int count_overflow = 0;
		    int max_overflow = 0;
		    int peak_overflow = 0;
		    int index_overflow = 0;
		    
		    // while reads are left on contig
		    while (jump != 0) {

		    	// Jump to desired region in bam
		    	if (!inFile.Jump(ref, jump)) {
		    		break;
		    	}
			    		

		    	// Initialize data arrays of positions, end positions, and strands
			    int positions[500] = {0}; 
				int end_positions[500] = {0}; // is not actually list of end positions, just list of greatest end pos for lsit subsets
				char strands[500];

				// Initialize adjacency matrix
				arma::mat adj_matrix(500, 500);
				adj_matrix.zeros(); 


				// Initialize loop variables
				int num_alignments = 0;
				int n;
				int temp_end;
				while (num_alignments < 500) {

					inFile.GetNextAlignment(alignment);

					// If next chromosome is reached, get out of town.
					if (alignment.RefID > ref) {
						break;
					}

					// Exclude secondary alignments
					if (!alignment.IsPrimaryAlignment() && (nonunique_alignments == false)) {
						continue;
					}

					// Check if sufficient mapping quality
					if (alignment.MapQuality <= mapq) {
			       		continue;
					}

					// Set values at appropriate position
					temp_end = alignment.GetEndPosition() - 1;
					end_positions[num_alignments] = (temp_end > end_positions[num_alignments - 1]) ? temp_end : end_positions[num_alignments -1]; 
					positions[num_alignments]  = round((end_positions[num_alignments] + alignment.Position) / 2); 	
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

				if (num_alignments < 500) {
					jump = 0;

				} else {
					jump = end_positions[total - 1] + 1;
				}

				
				// Initialize variables for counting reads by nodes in their graphs
				int i = 0;
				int x;
				int j;
				int increment;		// i offset
				int nodes;			// # of counts
				int max;			// max # of nodes
				int peak;			// location of peak

				while(i < total) {

					std::cerr << i << "\n";
					increment = i + 1;

					// if nodes meet min connectedness requirement
					if (degrees[i] >= min_cov) {

						// initialize group variables. Max and counts.
						// max serves as "peak" to represent graph. It's arbitrary.
						if (i == index_overflow && count_overflow != 0) {

							nodes = count_overflow;
							max = max_overflow;
							peak = peak_overflow;

							if (nodes > degrees[i]) {
								degrees[i] = nodes;
							}

							std::cerr << nodes << "\t" << max << "\t" << degrees[i] << "\n";

						} else {
							nodes = 1;
							max = degrees[i];
						}

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


									if (degrees[j] > max) {
										max = degrees[j];
										peak = positions[j];

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

									if (end_positions[j] > end_positions[x]) {
										x++;
									}

									if (degrees[j] > max) {
										max = degrees[j];
										peak = positions[j];

										if (max_overflow != 0) {
											std::cerr << degrees[j] << "\t" << degrees[x] << "\t" << max << "\n";
										}
									}
								
								}
							}

							j++;

						}


						// increment increment :)
						increment = j + 1;
					
						// report read counts
						if (increment < total) {
							std::cout << contig_cache[ref] << "\t" << peak << "\t" <<  nodes << "\n";
						}
					}

					i = increment;


		   		}


		   		if (degrees[total - 1] != 0) {

		   			inFile.GetNextAlignment(alignment);

		   			std::cerr << "overflow detected:\n";

		   			while (true) {

			   			char next_strand = (alignment.IsReverseStrand()) ? '-' : '+';

			   			if (next_strand == strands[total - 1]) {

				   			// check if overlap with next read; 
				   			if (alignment.Position <= end_positions[total - 1]) {
								count_overflow = nodes;
								max_overflow = max;
								peak_overflow = peak;
								break;


							} else {													
								std::cout << contig_cache[ref] << "\t" << peak << "\t" <<  nodes << "\n";
								count_overflow = 0;
								index_overflow = 0;
								break;

							}

						} else {
							index_overflow ++;
							continue;

						}

					}

		   		} else {

		   			index_overflow++;

		   		}

		   	}

		}			


};

