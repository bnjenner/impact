//////////////////////////////////////
// Node Class (basically a node in a doubly linked list)
class Node {

	public:

		int strand = -1;
		int read_count = -1;
		int clust_count = -1;
		int chrom_index;
		std::string gene_id = "";
		std::string chrom = "";
		std::string assigned_gene = "";

		// Clusters (heap allocation, minimizing can improve performance)
		std::vector<int> clust_vec{-1, -1};  // array of cluster start and stops (evens are starts, odds are ends)
		std::vector<int> count_vec{-1};      // array of cluster start and stops (evens are starts, odds are ends)

		// Links
		Node *next = NULL;
		Node *prev = NULL;

		// Node Variables
		int ishead = 0;
		int printed = 0;
		int assigned = 0;
		int ambiguous = 0;
		int max_overlap = 0;

		// Empty
		Node() {}

		// Initialized (alignment)
		Node(BamTools::BamAlignment &alignment, int ref_num) {

			strand = alignment.IsReverseStrand();
			read_count = 1;
			chrom_index = ref_num;
			clust_vec[0] = alignment.Position;

			calculate_splice(alignment, clust_vec); // Calculate spliced alignments

			clust_count = clust_vec.size() / 2;
			count_vec[0] = 1;

			count_vec.reserve(count_vec.size() + clust_count - 1);
			for (int i = 1; i < clust_count; i++) {
				count_vec.emplace_back(1);
			}
		}

		// Read cluster Initialized (region properties)
		Node(BamTools::BamAlignment &alignment, std::vector<int> &temp_vec, int ref_num) {
			strand = alignment.IsReverseStrand();
			read_count = 1;
			chrom_index = ref_num;
			clust_vec = temp_vec;
			clust_count = clust_vec.size() / 2;
			
			count_vec[0] = 1;
			count_vec.reserve(count_vec.size() + clust_count - 1);
			for (int i = 1; i < clust_count; i++) {
				count_vec.emplace_back(1);
			}
		}

		// Annotation Initialzed
		Node(std::string temp_gene_id, int temp_strand, int begin, int end, std::string temp_chrom) {
			gene_id = temp_gene_id;
			strand = temp_strand;
			chrom = temp_chrom;
			read_count = 0;
			clust_count = 1;
			clust_vec = std::vector<int>{begin, end};
		}


		// Reset
		void reset() {
			strand = -1;
			read_count = -1;
			clust_vec = {-1, -1}; 
			count_vec = {-1};
			clust_count = -1;
		}

		int get_start() { return clust_vec[0]; }
		int get_stop() { return clust_vec[2 * (clust_count - 1) + 1]; }

		int get_total_len() {
			int len = 0;
			for (int i = 0; i < clust_count; i++) {
				len += (clust_vec[(2 * i) + 1] - clust_vec[(2 * i)]);
			}
			return len;
		}

		void set_next(Node *node) { next = node; }
		void set_prev(Node *node) { prev = node; }

		// Calculate splice
		void calculate_splice(BamTools::BamAlignment &alignment, std::vector<int> &temp_vec) {
			
			int pos = 1;
			int inc = 0;

			for (int i = 0; i < alignment.CigarData.size(); i++) {
			
				// If gap is encounterd, add splice, (may also need to check cigar string standards)
				if (alignment.CigarData[i].Type == 'N') {

					temp_vec[pos] = temp_vec[pos - 1] + inc - 1;
					temp_vec.reserve(temp_vec.size() + 2);

					// expand vector, slow, will improve (maybe)
					temp_vec.emplace_back(temp_vec[pos] + alignment.CigarData[i].Length + 1);
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

		// Check overlap
		int check_overlap(int &temp_start, int &temp_stop, int &temp_strand) {

			if (strand != temp_strand) { return 0; }

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

		// Check genes
		int check_genes(int &temp_start, int &temp_stop, int &temp_strand) {

			// is strand correct
			if (strand != temp_strand) { return 0; }

			int max_overlap = 0;


			for (int i = 0; i < clust_count; i++) {

				// if subcluster begins within exon
				if ((clust_vec[(i * 2)] <= temp_start) && (clust_vec[(i * 2) + 1] >= temp_start)) {

					// if subcluster is entirely within exon
					if (clust_vec[(i * 2) + 1] >= temp_stop) {
						max_overlap = 2;
					} else {
						max_overlap = std::max(max_overlap, 1);
					}
				
				// if subcluster starts before exon but overlaps with it 
				} else if ((clust_vec[(i * 2)] > temp_start) && (clust_vec[(i * 2)] <= temp_stop)) {
					max_overlap = std::max(max_overlap, 1);
				}
			}

			return max_overlap;
		}

		// insert another spliced region
		void insert_splice(int int_start, int int_stop, int new_pos, int new_val) {

			std::vector<int> just_a_copy = clust_vec;
			std::vector<int> just_a_copy_1 = count_vec;

			clust_vec.reserve(clust_vec.size() + 2);  // insert new start and stop

			clust_vec.emplace_back(int_start);
			clust_vec.emplace_back(int_stop);

			std::sort(clust_vec.begin(), clust_vec.end());

			// insert new position in count vector (probably not ideal)
			std::vector<int>::iterator it;
			it = count_vec.begin() + new_pos;
			it = count_vec.insert(it, new_val);
		}

		// Delete spliced region
		void delete_splice(int i, int int_start, int int_stop) { 

		    int close_stop = 0;

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

				int new_sum = 0;
				int offset = -1;
				std::vector<int> temp_count_vec(clust_vec.size() / 2, -1);

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

		// check how read fits into clusters
		void modify_cluster(int &temp_start, int &temp_stop, int temp_count) {

			if (temp_start == -1) { return; }

			for (int i = 0; i < clust_count; i++) {

				// if read precedes cluster exit
				if (temp_stop < clust_vec[(i * 2)]) { 
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

						if (i != clust_count - 1) {
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


		// report cluster and counts
		int print_cluster(const std::string &contig_name, const ImpactArguments *parameters, int gene_count) {

			char s;
			std::string assignment;

			if ((clust_vec[0] == -1)) { return 0; }

			// Assign strand
			if (parameters -> stranded == "forward") {
				s = (strand == 1) ? '-' : '+';
			} else {
				s = (strand == 1) ? '+' : '-';
			}

			// assignment
			if (assigned == 0) {
				assignment = "__unassigned";
			} else {
				if (ambiguous == 1) {
					assignment = "__ambiguous";
				} else {
					assignment = "__assigned";
				}
			}

			std::ofstream outdata;
			outdata.open(parameters -> gtf_output, std::ios::out | std::ios::app);

			// Print "Gene" line, not contiguous
			outdata << contig_name << "\timpact\tcluster\t"
					  << clust_vec[0] + 1 << "\t" << clust_vec[((clust_count - 1) * 2) + 1] + 1
					  << "\t.\t" << s << "\t.\t"
					  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
					  << "subclusters \"" << clust_count << "\"; "
					  << "counts \"" << read_count << "\"; "
					  << "assignment \"" << assignment << "\"\n";
					  //<< "counts \"" << read_count << "\";\n";

			// Iterate through clusters
			for (int i = 0; i < clust_count; i++) {
				// Print name, strand, and first start
				outdata << contig_name << "\timpact\tsubcluster\t" << clust_vec[(i * 2)] + 1 << "\t" << clust_vec[(i * 2) + 1] + 1;
				// Print the rest lol
				outdata << "\t.\t" << s << "\t.\t"
						  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
						  << "subcluster_id \"impact." << contig_name << "." << gene_count << "." << i + 1 << "\"; " 
						  << "counts \"" << count_vec[i] <<"\";\n";
			}

			return 1;
		}
};

