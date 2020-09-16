using namespace BamTools;

class Overflow {

	public: 

	////////////////////////////
	// Attributes
		int count_overflow = 0;
		int max_overflow = 0;
		int peak_overflow = 0;
		int index_overflow = 0;

		int start = -1;
		int stop = -1;

	////////////////////////////
	// Methods

		///////////////////////
		// back to zero
		void reset() {

			count_overflow = 0;
			max_overflow = 0;
			peak_overflow = 0;
			index_overflow = 0;

			start = -1;
			stop = -1;

		}

		///////////////////////
		// update attributes	
		void update(int index, int count, int max, int peak, int begin, int end) {

			index_overflow = index;
			count_overflow = count;
			max_overflow = max;
			peak_overflow = peak;

			start = begin;
			stop = end;

		}

};



//////////////////////////////////////
// Graph Class
class Graph {

	public:

	////////////////////////////
	// Attributes

		// Arrays of read data
		int start_positions[5000] = {0};
		int end_positions[5000] = {0};
		char strands[5000] = {' '};


		// Initialize adjacency matrix and degree vector
		arma::mat adj_matrix;
		arma::rowvec degrees;


		// Useful global variables
		bool contig_end;
		int num_alignments;
		int ref; 
		std::string contig_name;


	////////////////////////////
	// Constructors

		Graph(int ref_num, std::string ref_name) {

			// Initialize matrix and contig number and name
			adj_matrix = arma::mat(5000, 5000, arma::fill::zeros);
			ref = ref_num;
			contig_name = ref_name;

		}


	////////////////////////////
	// Methods

		///////////////////////
		// Create Adjancecy Matrix and Populate Start / End / Char Arrays.

		void create_adjacency(BamReader &inFile, BamAlignment &alignment, Parameters &parameters) {

			// Initialize loop variables
			int n;
			int temp_end;
			num_alignments = 0;
			contig_end = false;

			while (num_alignments < 5000) {

				inFile.GetNextAlignment(alignment);

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					contig_end = true;
					break;
				}

				// Exclude secondary alignments
				if (!alignment.IsPrimaryAlignment() && (parameters.nonunique_alignments == false)) {
					continue;
				}

				// Check if sufficient mapping quality
				if (alignment.MapQuality <= parameters.mapq) {
		       		continue;
				}

				// get alignment end 
				temp_end = alignment.GetEndPosition() - 1;

				// If alignment spans double its actual length in bases, ignore it, its trash
				if ((2 * alignment.Length) <= (temp_end - alignment.Position)) {
					continue;
				}

				// read is ok, continue with populating vectors and matrix
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


		}


		///////////////////////
		// Gets degree of nodes

		void get_degrees() {

			degrees = arma::sum(adj_matrix, 0);
		}


		///////////////////////
		// gets next position to jump to in next iteration

		int get_jump() {

			int jump;

			// Find approprite spot to break, we can't cut a group in half can we?
			if (num_alignments < 5000 || contig_end) {
				jump = 0;

			} else {
				jump = end_positions[num_alignments - 1] + 1;
			}

			return jump;
		}


		///////////////////////
		// checks if last cluster is incomplete

		int check_overflow(BamReader &inFile, BamAlignment &alignment) {

			char next_strand;
			int index_overflow = 0;
			
			// Check if the last node checked concluded a group.
			if (degrees[num_alignments - 1] != 0) {

				// Check following alignments to see if and where groups stops
				while (true) {

					inFile.GetNextAlignment(alignment);
					next_strand = (alignment.IsReverseStrand()) ? '-' : '+';


					if (next_strand == strands[num_alignments - 1]) {

						// check if overlap with next read; 
						if (alignment.Position <= end_positions[num_alignments - 1]) {
							return index_overflow;


						// If read does not overlap, overflow variables are set to 0 and counts for last group are reported.
						} else {													
							return -1;

						}

					// if next read is wrong strand, keep searching
					} else {
						index_overflow ++;
						continue;

					}

				}

			// If degree of last node is zero, no need to report
			} else {
				return -1;

			}
		}


		///////////////////////
		// gets clusters by counting connected nodes and checking for overflows (this one's a bit rough....)

		void get_clusters(BamReader &inFile, BamAlignment &alignment, Parameters &parameters, Overflow &overflow) {

			// Initialize variables for counting reads by nodes in their graphs
			int i = 0;
			int x;
			int j;

			int increment;		// i offset
			int nodes;			// # of counts
			int max;			// max # of nodes
			int peak;			// location of peak
			char strand;		// group strand
			int start;			// group start
			int stop;			// group stop

			while(i < num_alignments) {

				increment = i + 1;

				// if nodes meet min connectedness requirement
				if (degrees[i] >= parameters.min_cov) {

					// initialize group variables. Max and counts.
					// max serves as "peak" to represent graph. It's arbitrary.
					if (i == overflow.index_overflow && overflow.count_overflow != 0) {

						nodes = overflow.count_overflow;
						max = overflow.max_overflow;
						peak = overflow.peak_overflow;
						start = overflow.start;
						stop = overflow.stop;

					} else {

						nodes = 1;
						max = degrees[i];
						peak = round((end_positions[i] + start_positions[i]) / 2); 
						start = start_positions[i];
						stop = end_positions[i];
					}

					strand = strands[i];

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
								start = start_positions[j];


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

					while (j < num_alignments && end_positions[x] > start_positions[j]) {


						if (strands[x] == strands[j]) {

							if (adj_matrix(x, j) == 0) {
								break;
						
							} else {

								nodes++;

								if (end_positions[j] > end_positions[x]) {
									x = j;
									stop = end_positions[j];
								}

								if (degrees[j] > max) {
									max = degrees[j];
									peak = round((end_positions[j] + start_positions[j]) / 2);

								}
							
							}
						}

						j++;

					}

					// increment increment :)
					increment = j + 1;

					// report read counts
					if (increment < num_alignments) {
						std::cout << contig_name << "\t" << strand << "\t" << start << "\t"
									<< stop << "\t" << peak << "\t" <<  nodes << "\n";
					}
				}

				i = increment;


			}
			
			int incomplete_offset = this -> check_overflow(inFile, alignment);

			if (incomplete_offset != -1) {
				overflow.update(incomplete_offset, nodes, max, peak, start, stop);
			
			} else {
				std::cout << contig_name << "\t" << strand << "\t" << start << "\t"
							<< stop << "\t" << peak << "\t" <<  nodes << "\n";
				overflow.reset();
			
			}

		}
};

