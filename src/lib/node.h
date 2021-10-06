using namespace BamTools;

//////////////////////////////////////
// Node Class (basically a node in a doubly linked list)
class Node {

	public:

	////////////////////////////
	// Attributes

		// Cluster variables
		int strand = -1;
		int read_count = -1;
		int clust_count = -1;
		int chrom_index;

		std::string gene_id = "";
		std::string chrom = "";

		// Clusters
		//    (heap allocation, minimizing can improve performance)
		std::vector<int> clust_vec{-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)

		std::vector<int> count_vec{-1}; // array of cluster start and stops (evens are starts, odds are ends)

		// Links
		Node *next = NULL;
		Node *prev = NULL;

		// Node Variables
		int ishead = 0;
		int printed = 0;
		int ambiguous = 0;


		int max_overlap = 0;

	////////////////////////////
	// Constructors

		// Empty
		Node() {}


		////////////////////////////
		// ALIGNMENT INIT 

		// Initialized (alignment)
		Node(BamAlignment &alignment, int ref_num) {

			// Get cluster properties
			strand = alignment.IsReverseStrand();
			read_count = 1;
			chrom_index = ref_num;

			// Start position
			clust_vec[0] = alignment.Position;

			// Calculate spliced alignments
			calculate_splice(alignment, clust_vec);
			clust_count = clust_vec.size() / 2;

			count_vec[0] = 1;
			// if splice, add read
			for (int i = 1; i < clust_count; i++) {
				count_vec.push_back(1);
			}

		}

		// Read cluster Initialized (region properties)
		Node(std::vector<int> &temp_vec, int temp_strand, int ref_num) {

			// Get cluster properties
			strand = temp_strand;
			read_count = 1;
			chrom_index = ref_num;
		
			// copy vector
			clust_vec = temp_vec;
			clust_count = clust_vec.size() / 2;

			count_vec[0] = 1;
			// if splice add read
			for (int i = 1; i < clust_count; i++) {
				count_vec.push_back(1);
			}
		}


		////////////////////////////
		// ANNOTATION INIT

		// Annotation Initialzed
		Node(std::string temp_gene_id, int temp_strand, std::string temp_chrom) {

			gene_id = temp_gene_id;
			strand = temp_strand;
			chrom = temp_chrom;
			read_count = 0;
			clust_count = 1;

			clust_vec = std::vector<int>{0, 0};
		}


	////////////////////////////
	// Methods

		////////////////////////////
		// Reset
		void reset() {

			// Cluster variables
			strand = -1;
			read_count = -1;

			// Clusters
			clust_vec = {-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)
			count_vec = {-1};
			clust_count = -1;
		}


		////////////////////////////
		// get first position
		int get_start() {
			return clust_vec[0];
		}

		////////////////////////////
		// get last position
		int get_stop() {
			return clust_vec[2 * (clust_count - 1) + 1];
		}

		////////////////////////////
		// get total length
		int get_total_len() {

			int len = 0;
			for (int i = 0; i < clust_count; i++) {
				len += (clust_vec[(2 * i) + 1] - clust_vec[(2 * i)]);
			}

			return len;
		}

		////////////////////////////
		// Set Next
		void set_next(Node *node) {
			next = node;
		}

		////////////////////////////
		// Set Prev
		void set_prev(Node *node) {
			prev = node;
		}


		////////////////////////////
		// Calculate splice
		void calculate_splice(BamAlignment &alignment, std::vector<int> &temp_vec) {
			
			int pos = 1;
			int inc = 0;

			// iterate through CIGAR string
			for (int i = 0; i < alignment.CigarData.size(); i++) {
			
				// If gap is encounterd, add splice, (may also need to check cigar string standards)
				if (alignment.CigarData[i].Type == 'N') {

					temp_vec[pos] = temp_vec[pos - 1] + inc - 1;

					temp_vec.reserve(temp_vec.size() + 2);

					// expand vector, slow, will improve (maybe)
					temp_vec.emplace_back(temp_vec[pos] + alignment.CigarData[i].Length);
					temp_vec.emplace_back(-1);

					
					pos += 2;
					inc = 0;

				// If not gapped, add to start position
				} else if ((alignment.CigarData[i].Type) != 'S' && (alignment.CigarData[i].Type) != 'H' &&
					       (alignment.CigarData[i].Type) != 'I') {					
					inc += alignment.CigarData[i].Length;

				}

				// if end of cigar string is reached
				if (i == alignment.CigarData.size() - 1) {
					temp_vec[pos] = temp_vec[pos - 1] + inc - 1;

				}
			
			}
		
		}


		////////////////////////////
		// Check overlap
		int check_overlap(int &temp_start, int &temp_stop, int &temp_strand) {

			// is strand correct
			if (strand != temp_strand) {
				return 0;
			}

			// iterate through all clusters
			for (int i = 0; i < clust_count; i++) {

				// check if beginning of read exists within a cluster
				if ((temp_start >= clust_vec[i * 2]) && (temp_start <= clust_vec[(i * 2) + 1])) {
					return 1;

				// check if end of read exists within a cluster				
				} else if ((temp_stop >= clust_vec[i * 2]) && (temp_stop <= clust_vec[(i * 2) + 1])) {
					return 1;

				// in read spans cluster 
				} else if ((temp_start <= clust_vec[i * 2]) && (temp_stop >= clust_vec[(i * 2) + 1])) {
					return 1;

				}

			}

			return 0;
		}


		////////////////////////////
		// Check genes
		int check_genes(int &temp_start, int &temp_stop, int &temp_strand) {

			// is strand correct
			if (strand != temp_strand) {
				return 0;
			}

			// float cov_threshold = read_count * static_cast<double>(0.50);
			// cov_threshold = 0.0;

			int overlap_score = 0;

			// iterate through all clusters
			for (int i = 1; i < clust_count; i++) {

				// if subcluster begins within exon
				if ((clust_vec[(i * 2)] <= temp_start) && (clust_vec[(i * 2) + 1] >= temp_start)) {

					// if subcluster is entirely within exon
					if (clust_vec[(i * 2) + 1] >= temp_stop) {
						overlap_score += (temp_stop - temp_start) + 1;
						// return 2;
						// max_overlap = std::max(max_overlap, 2);
					} else {
						overlap_score += (clust_vec[(i * 2) + 1] - temp_start) + 1;
						//max_overlap = std::max(max_overlap, 1);
					}
				
				// if subcluster starts before exon but overlaps with it 
				} else if ((clust_vec[(i * 2)] > temp_start) && (clust_vec[(i * 2)] <= temp_stop)) {

					if (temp_stop > clust_vec[(i * 2) + 1]) {
						overlap_score += (clust_vec[(i * 2) + 1] - clust_vec[(i * 2)]) + 1;
					
					} else {
						overlap_score += (temp_stop - clust_vec[(i * 2)]) + 1;
					}
					//max_overlap = std::max(max_overlap, 1);
					//return 1;
				}

					//////////////////////////////
					// (may reimplement)					
					// // if subcluster starts before exon but overlaps with it 
					// } else if ((clust_vec[(i * 2)] < temp_start) && (clust_vec[(i * 2) + 1] >= temp_start)) {
					// 	return 1;
					// }

					// check if beginning of read exists within a cluster
					// if ((temp_start >= clust_vec[i * 2]) && (temp_start <= clust_vec[(i * 2) + 1])) {
					// 	return 1;

					// // check if end of read exists within a cluster				
					// } else if ((temp_stop >= clust_vec[i * 2]) && (temp_stop <= clust_vec[(i * 2) + 1])) {
					// 	return 1;

					// // in read spans cluster 
					// } else if ((temp_start <= clust_vec[i * 2]) && (temp_stop >= clust_vec[(i * 2) + 1])) {
					// 	return 1;
					//////////////////////////////
				
					// }
			}

			//max_overlap = std::max(max_overlap, 1);

			return overlap_score;
		}




		////////////////////////////
		// insert another spliced region
		void insert_splice(int int_start, int int_stop, int new_pos, int new_val) {

			std::vector<int> just_a_copy = clust_vec;
			std::vector<int> just_a_copy_1 = count_vec;

			// insert new start and stop
			clust_vec.reserve(clust_vec.size() + 2);

			clust_vec.emplace_back(int_start);
			clust_vec.emplace_back(int_stop);

			std::sort(clust_vec.begin(), clust_vec.end());

			// insert new position in count vector (probably not ideal)
			std::vector<int>::iterator it;
			it = count_vec.begin() + new_pos;
			it = count_vec.insert(it, new_val);

		}

		////////////////////////////
		// Delete spliced region
		void delete_splice(int i, int int_start, int int_stop) { 

		    int close_stop = 0;

			// Iterate through remaining clusters
			for (int j = i + 1; j < clust_count; j++) {

				// Determine if clusters are joined by read
				if (clust_vec[(j * 2)] > int_stop) {
					break;
				}

				close_stop = j;	
			} 

			if (close_stop != 0 && close_stop != i) {

				int start_index = (i * 2);
				int stop_index = (2 * close_stop) + 1;

				clust_vec.erase(clust_vec.begin() + start_index + 1, clust_vec.begin() + stop_index);

				// new sum
				int new_sum = 0;
				int offset = -1;

				// create temp vector
				std::vector<int> temp_count_vec(clust_vec.size() / 2, -1);

				// populate values
				for (int x = 0; x < count_vec.size(); x++) {

					if (x < i || x > close_stop) {
						temp_count_vec[x - std::max(offset, 0)] = count_vec[x];
					} else {
						new_sum += count_vec[x];
						offset++;
					}
				
				}

				temp_count_vec[i] = new_sum;
				count_vec = temp_count_vec;
			}
		} 

		////////////////////////////
		// check how read fits into clusters
		void modify_cluster(int &temp_start, int &temp_stop, int temp_count) {


			if (temp_start == -1) {
				return;
			}

			// iterate through all clusters
			for (int i = 0; i < clust_count; i++) {

				// if read precedes cluster exit
				if (temp_stop < clust_vec[(i * 2)]) { 

					// insert region
					insert_splice(temp_start, temp_stop, i, temp_count);
					break;
				
				// if read follows cluster, go to next cluster
				} else if (temp_start > clust_vec[(i * 2) + 1]) {

					// if last cluster add
					if (i == clust_count - 1) {
						insert_splice(temp_start, temp_stop, i + 1, temp_count);
						break;
					}
				
				} else {

					// Updates beginning of cluster to longest value betweeen end of cluster and end of read
					clust_vec[(i * 2)] = (clust_vec[(i * 2)] < temp_start) ? clust_vec[(i * 2)] : temp_start;
							
					// if end of read extends past than first chunk
					if (temp_stop > clust_vec[(i * 2) + 1]) {

						// if not last read cluster
						if (i != clust_count - 1) {

							// Checks if clusters are joined by read
							delete_splice(i, temp_start, temp_stop);
						}

						// Updates end of cluster to longest value betweeen end of cluster and end of read
						clust_vec[(i * 2) + 1] = (clust_vec[(i * 2) + 1] > temp_stop) ? clust_vec[(i * 2) + 1] : temp_stop;
					}

					count_vec[i] += temp_count;
					break;

				}

			}

			clust_count = clust_vec.size() / 2;

		}

		// void filter_clusters(const Parameters &parameters) {

		// 	// // Iterate through clusters
		// 	// int i = 0;

		// 	// while (i < clust_count) {

		// 	// 	if ((count_vec[i] < (parameters.min_cov) )) {
		// 	// 		clust_vec.erase(clust_vec.begin() + (2 * i), clust_vec.begin() + (2 * i) + 2);
		// 	// 		count_vec.erase(count_vec.begin() + i, count_vec.begin() + i + 1);
		// 	// 		clust_count--;

		// 	// 	} else {
		// 	// 		i++;
				
		// 	// 	}

		// 	// 	// if cluster is too long, pass to subdivide script, 
		// 	// 	//	probably will be implemented in python
		// 	// 	// if (((clust_vec[(2 * i) + 1] - clust_vec[(2 * i)]) > (width * 2)) &&
		// 	// 	// 	  count_vec[i] > 10) {
		// 	// 	// 	std::cerr << clust_vec[(2 * i)] << "\t" << clust_vec[(2 * i) + 1] << "\n";

		// 	// 	// } 

		// 	// }
		// }

		////////////////////////////
		// report cluster and counts
		int print_cluster(const std::string &contig_name, const Parameters &parameters, int gene_count) {

			// strand character
			char s;

			// If empty or does not meet miniumum coverage
			if ((clust_vec[0] == -1)) {//|| read_count < parameters.min_cov) {
				return 0;
			}

			// Assign strand
			if (parameters.stranded == 'f') {
				s = (strand == 1) ? '-' : '+';
			} else {
				s = (strand == 1) ? '+' : '-';
			}



			// ////////////////////////////////
			// // for testing
			// if (clust_count > 1) {
			// 	return 1;
			// }
			// ////////////////////////////////


			// Print "Gene" line, not contiguous
			std::cout << contig_name << "\timpact\tcluster\t"
					  << clust_vec[0] + 1 << "\t" << clust_vec[((clust_count - 1) * 2) + 1]
					  << "\t.\t" << s << "\t.\t"
					  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
					  << "subclusters \"" << clust_count << "\";\n";
					  //<< "counts \"" << read_count << "\";\n";

			// Iterate through clusters
			for (int i = 0; i < clust_count; i++) {

				// Print name, strand, and first start
				std::cout << contig_name << "\timpact\tsubcluster\t";
				std::cout << clust_vec[(i * 2)] + 1 << "\t" << clust_vec[(i * 2) + 1];

				// Print the rest lol
				std::cout << "\t.\t" << s << "\t.\t"
						  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
						  << "subcluster_id \"" << i + 1 << "\"; " 
						  << "counts \"" << count_vec[i] <<"\";\n";

			}

			return 1;
		}

		void exon_print() {

			std::cerr << gene_id << "\n";

			for (int i = 1; i < clust_count; i++) {
				std::cerr << clust_vec[(2 * i)] << "\t" 
						  << clust_vec[(2 * i) + 1] << "\n";
			}

		}

};

