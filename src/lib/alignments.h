#include "api/BamReader.h"

using namespace BamTools;

//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

	////////////////////////////
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


   		
    	// Contig Cache
    	std::unordered_map<int, std::string> contig_cache; 	// Unordered map 


    ////////////////////////////
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

    ////////////////////////////
    // Methods

    	///////////////////////
    	// Open files
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


		///////////////////////
		// Close files
		void close() {

			inFile.Close();

		}


		///////////////////////
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


		///////////////////////
		// Create Adjancecy Matrix and Populate Start / End / Char Arrays.
		int get_adj_matrix(std::vector<int> &start_positions, std::vector<int> &end_positions, 
						  std::vector<char> &strands, bool &contig_end, int ref, arma::mat &adj_matrix) {


			// Initialize loop variables
			int num_alignments = 0;
			int n;
			int temp_end;

			while (num_alignments < 10000) {

				inFile.GetNextAlignment(alignment);

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					contig_end = true;
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
				start_positions[num_alignments] = alignment.Position;
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


			return num_alignments;

		}


		///////////////////////
		// Grab Alignments within Interval Using Bam Index
		void get_counts(int ref, int pos) {

			// Variable accounting for overflow when group is cut off by chunking reads
		    int jump = pos;
		    int count_overflow = 0;
		    int max_overflow = 0;
		    int peak_overflow = 0;
		    int index_overflow = 0;
		    

		    // while still reading reads on desired contig
		    while (jump != 0) {


		    	// Jump to desired region in bam
		    	if (!inFile.Jump(ref, jump)) {
		    		std::cerr << "[ERROR: Could not jump to region: " << ref << ":" << jump << ".]\n";
		    		break;
		    	}
			    		

		    	// Initialize data arrays of positions, end positions, and strands
				std::vector<int> start_positions(10000, 0);
				std::vector<int> end_positions(10000, 0);
				std::vector<char> strands(10000, ' ');

				// Initialize adjacency matrix
				arma::mat adj_matrix(10000, 10000);
				adj_matrix.zeros(); 

				// Variable determining if end of contig was reached
				bool contig_end = false;

				// Create adjency matrix and get numver of aligned reads
				int num_alignments = get_adj_matrix(start_positions, end_positions,
													strands, contig_end, ref,
													adj_matrix);

				
				// Calculate degree of read node
				arma::rowvec degrees = arma::sum(adj_matrix, 0);


				// Find approprite spot to break, we can't cut a group in half can we?
				if (num_alignments < 10000 || contig_end) {
					jump = 0;

				} else {
					jump = end_positions[num_alignments - 1] + 1;
				}

				
				// Initialize variables for counting reads by nodes in their graphs
				int i = 0;
				int x;
				int j;
				int increment;		// i offset
				int nodes;			// # of counts
				int max;			// max # of nodes
				int peak;			// location of peak

				while(i < num_alignments) {

					increment = i + 1;

					// if nodes meet min connectedness requirement
					if (degrees[i] >= min_cov) {

						// initialize group variables. Max and counts.
						// max serves as "peak" to represent graph. It's arbitrary.
						if (i == index_overflow && count_overflow != 0) {

							nodes = count_overflow;
							max = max_overflow;
							peak = peak_overflow;

						} else {

							nodes = 1;
							max = degrees[i];
							peak = round((end_positions[i] + start_positions[i]) / 2); 

						}


						// Search upstream of node
						x = i;
						j = i - 1;
						while (j >= 0 && end_positions[j] > start_positions[x]) {

							if (strands[x] == strands[j]) {

								if (adj_matrix(x, j) == 0) {
									break;
							
								} else {

									nodes++;
									x = j;


									if (degrees[j] > max) {
										max = degrees[j];
										peak = round((end_positions[j] + start_positions[j]) / 2);

									}
								
								}
							}

							j--;

						}

						// Search downstream of node
						x = i;
						j = i + 1;

						while (j < total && end_positions[x] > start_positions[j]) {


							if (strands[x] == strands[j]) {

								if (adj_matrix(x, j) == 0) {
									break;
							
								} else {

									nodes++;

									if (end_positions[j] > end_positions[x]) {
										x = j;
									}

									if (degrees[j] > max) {
										max = degrees[j];
										peak = peak = round((end_positions[j] + start_positions[j]) / 2);

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


		   		// Check if the last node checked concluded a group.
		   		if (degrees[total - 1] != 0) {

		   			// Check preceding alignments to see if and where groups stops
		   			while (true) {

		   				inFile.GetNextAlignment(alignment);
			   			char next_strand = (alignment.IsReverseStrand()) ? '-' : '+';


			   			if (next_strand == strands[total - 1]) {

				   			// check if overlap with next read; 
				   			if (alignment.Position <= end_positions[total - 1]) {
								count_overflow = nodes;
								max_overflow = max;
								peak_overflow = peak;
								break;


							// If read does not overlap, overflow variables are set to 0 and counts for last group are reported.
							} else {													
								std::cout << contig_cache[ref] << "\t" << peak << "\t" <<  nodes << "\n";
								count_overflow = 0;
								index_overflow = 0;
								peak_overflow = 0;
								break;

							}

						// if next read is wrong strand, keep searching
						} else {
							index_overflow ++;
							continue;

						}

					}

				// If degree of last node is zero, no need to report
		   		} else {
		   			count_overflow = 0;
		   			index_overflow = 0;
		   			peak_overflow = 0;


		   		}

		   	} 

		}			


};

