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

		// Clusters
		//    (heap allocation, minimizing can improve performance)
		std::vector<int> clust_vec{-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)

		// Links
		Node *next = NULL;
		Node *prev = NULL;

		// Node Variables
		int ishead = 0;
		int printed = 0;


	////////////////////////////
	// Constructors

		// Empty
		Node() {}

		// Initialized (alignment)
		Node(BamAlignment &alignment) {

			// Get cluster properties
			strand = alignment.IsReverseStrand();
			read_count = 1;

			// Start position
			clust_vec[0] = alignment.Position;

			// Calculate spliced alignments
			calculate_splice(alignment, clust_vec);
			clust_count = clust_vec.size() / 2;

		}

		// Initialized (region properties)
		Node(std::vector<int> &temp_vec, int temp_strand) {

			// Get cluster properties
			strand = temp_strand;
			read_count = 1;

			// copy vector
			clust_vec = temp_vec;
			clust_count = clust_vec.size() / 2;
		}


	////////////////////////////
	// Methods

		////////////////////////////
		// Reset
		void reset() {

			// Cluster variables
			strand = -1;
			read_count = 1;

			// Clusters
			clust_vec = {-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)
			clust_count = 1;
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

					temp_vec[pos] = temp_vec[pos - 1] + inc;

					temp_vec.reserve(temp_vec.size() + 2);

					// expand vector, slow, will improve (maybe)
					temp_vec.emplace_back(temp_vec[pos] + alignment.CigarData[i].Length);
					temp_vec.emplace_back(-1);

					
					pos += 2;
					inc = 0;

				// If not gapped, add to start position
				} else if (alignment.CigarData[i].Type != 'S' && alignment.CigarData[i].Type != 'H') {					
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
		// insert another spliced region
		void insert_splice(int int_start, int int_stop) {

			clust_vec.reserve(clust_vec.size() + 2);

			clust_vec.emplace_back(int_start);
			clust_vec.emplace_back(int_stop);

			std::sort(clust_vec.begin(), clust_vec.end());

		}

		////////////////////////////
		// Delete spliced region
		void delete_splice(int i, int int_stop) { 

		    int factor = 0; 

			// Iterate through remaining clusters
			for (int j = i + 1; j < clust_count; j++) {

				// Determine if clusters are joined by read
				if (clust_vec[(j * 2)] > int_stop) {
					break;
				}

				factor = j;	
			} 

			if (factor != 0) {

				int start_index = (i * 2);
				int stop_index = (2 * factor) + 1;

				clust_vec.erase(clust_vec.begin() + start_index + 1, clust_vec.begin() + stop_index);
			
			}

		} 

		////////////////////////////
		// check how read fits into clusters
		void modify_cluster(int &temp_start, int &temp_stop) {

			if (temp_start == -1) {
				return;
			}

			// iterate through all clusters
			for (int i = 0; i < clust_count; i++) {

				// if read precedes cluster exit
				if (temp_stop < clust_vec[(i * 2)]) { 

					// insert region
					insert_splice(temp_start, temp_stop);
					break;
				
				// if read follows cluster, go to next cluster
				} else if (temp_start > clust_vec[(i * 2) + 1]) {

					// if last cluster add
					if (i == clust_count - 1) {
						insert_splice(temp_start, temp_stop);
						break;
					}
				
				} else {

					// Updates beginning of cluster to longest value betweeen end of cluster and end of read
					clust_vec[(i * 2)] = (clust_vec[(i * 2)] < temp_start) ? clust_vec[(i * 2)] : temp_start;
							
					// if end of read is greater than first chunk
					if (temp_stop > clust_vec[(i * 2) + 1]) {

						// Checks if clusters are joined by read
						delete_splice(i, temp_stop);

						// Updates end of cluster to longest value betweeen end of cluster and end of read
						clust_vec[(i * 2) + 1] = (clust_vec[(i * 2) + 1] > temp_stop) ? clust_vec[(i * 2) + 1] : temp_stop;


					}

					break;

				}

			}
						
			clust_count = clust_vec.size() / 2;

		}

		////////////////////////////
		// report cluster and counts
		int print_cluster(const std::string &contig_name, const Parameters &parameters, int gene_count) {

			// strand character
			char s;

			// If empty or does not meet miniumum coverage
			if ((clust_vec[0] == -1) || read_count < parameters.min_cov) {
				return 0;
			}

			// Assign strand
			if (parameters.stranded == 'f') {
				s = (strand == 1) ? '-' : '+';
			} else {
				s = (strand == 1) ? '+' : '-';
			}

			// Print "Gene" line, not contiguous
			std::cout << contig_name << "\timpact\tcluster\t"
					  << clust_vec[0] + 1 << "\t" << clust_vec[((clust_count - 1) * 2) + 1] + 1
					  << "\t.\t" << s << "\t.\t"
					  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
					  << "subclusters \"" << clust_count << "\"; "
					  << "counts \"" << read_count << "\";\n";

			// Iterate through clusters
			for (int i = 0; i < clust_count; i++) {

				// Print name, strand, and first start
				std::cout << contig_name << "\timpact\tsubcluster\t";
				std::cout << clust_vec[(i * 2)] + 1 << "\t" << clust_vec[(i * 2) + 1] + 1;

				// Print the rest lol
				std::cout << "\t.\t" << s << "\t.\t"
						  << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
						  << "subcluster_id \"" << i + 1 << "\";\n";

			}

			return 1;
		}

};



//////////////////////////////////////
// Graph Class (really just a doubly linked list)
class Graph {

	public:

	////////////////////////////
	// Attributes

		// Useful global variables
		int ref; 
		std::string contig_name;
		Parameters parameters;

		// Import node variables
		Node *head = NULL;
		Node *tail = NULL;
		Node temp;


	////////////////////////////
	// Constructors

		// Empty
		Graph() {}


	////////////////////////////
	// Methods

		///////////////////////
		// Initialize empty object
		void initialize(int ref_num, std::string ref_name, const Parameters &pre_parameters) {

			// Initialize contig number and name
			ref = ref_num;
			contig_name = ref_name;
			parameters = pre_parameters;
		}		


		///////////////////////
		// Set Head
		int set_head(BamReader &inFile, BamAlignment &alignment) {

			// Find suitable first alignment
			while (true) {

				// If no alignments
				if (!inFile.GetNextAlignment(alignment)){
					//std::cerr << "[ERROR: No Alignment]\n";
					return 0;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
			        //std::cerr << "[Finished Counting from " << contig_name << "...]\n";
					return 0;
				}

				// If alignment is a duplicate
				if (alignment.IsDuplicate()) {
					continue;
				}

				// Exclude secondary alignment 
				if (!alignment.IsPrimaryAlignment() && (parameters.nonunique_alignments == false)) {
					continue;
				} 
	
				// Check if sufficient mapping quality
				if (alignment.MapQuality < parameters.mapq) {
		       		continue;
				}

				break;
			}

			// Add alignment to head node
			temp = Node(alignment);

			temp.ishead = 1;
			head = &temp;

			return 1;

		}


		///////////////////////
		// Create clusters of overlapping reads
		void create_clusters(BamReader &inFile, BamAlignment &alignment) {

			// Initialize loop variables
			int regions;
			int temp_start; // used to kill one of loops below
			int temp_strand;
			std::vector<int> temp_vec = {-1, -1};

			// Initialize pointers
			Node *curr_node;
			tail = head;


			while (true) {

				// Start at last node
				curr_node = tail;

				// End of File Reached
				if (!inFile.GetNextAlignment(alignment)) {
					std::cerr << "[End of File Reached!]\n";
					return;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					return;
				}

				// If alignment is a duplicate
				if (alignment.IsDuplicate()) {
					continue;
				}

				// Exclude secondary alignments
				if (!alignment.IsPrimaryAlignment() && (parameters.nonunique_alignments == false)) {
					continue;
				}

				// If paired end, check propper pair
				if (!alignment.IsProperPair() && (parameters.library_type == 'p')) {
					continue;
				}

				// Check if sufficient mapping quality
				if (alignment.MapQuality < parameters.mapq) {
		       		continue;
				}

				// write to temp vector
				temp_vec = {alignment.Position, -1};

				// get alignment start
				temp_start = temp_vec[0];

				// get strand
				temp_strand = alignment.IsReverseStrand();

				// calculate splice sites
				curr_node -> calculate_splice(alignment, temp_vec);

				// check if alignment represents a new node (past first subcluster)
				if ((temp_vec[0] > curr_node -> clust_vec[1]) || (temp_strand != curr_node -> strand)) {

					// Create node
					Node *new_node = new Node(temp_vec, temp_strand);

					// link nodes within graph
					curr_node -> set_next(new_node);
					new_node -> set_prev(curr_node);
					curr_node = new_node;
					tail = curr_node;

					continue;
				
				}

				// number of aligned regions
				regions = temp_vec.size() / 2;
				
				// find overlapping region
				while ((curr_node != NULL) && (temp_start < curr_node -> get_stop()))  {

					for (int x = 0; x < regions; x++) {

						// Check if alignment overlaps with previous nodes
						if (curr_node -> check_overlap(temp_vec[(2 * x)], temp_vec[(2 * x) + 1], temp_strand)) {

							// add all clusters to vector
							for (int y = 0; y < regions; y++) {
								curr_node -> modify_cluster(temp_vec[(2 * y)], temp_vec[(2 * y) + 1]);
							}

							curr_node -> read_count++;

							// kill the loop
							temp_start = curr_node -> get_stop() + 1;
							break;
						}
					
					}

					curr_node = curr_node -> prev;

				}

			}

			return;
		}

		///////////////////////
		// print clusters in graph
		void collapse_graph() {

			// initialize pointer (I used a lot don't judge me)
			Node *curr_node = head;
			Node *next_node = NULL;
			Node *temp_node = NULL;
			Node *inter_node = NULL;


			// throw away variables, function only takes references :( will work on this
			int t_start;
			int t_stop;
			int t_strand;
			int t_overlap;
			int t_restart;

			// check if head node exists
			if (head == NULL) {
				return;
			}


			// Iterate through nodes
			while (curr_node != NULL) {

				next_node = curr_node -> next;
				temp_node = curr_node -> next;

				t_restart = 0;

				// if node hasn't already been printed / added to another cluster
				if (curr_node -> printed == 0) {
					

					while (true) {

						// checks if curr_node is last node or beyond range
						if ((temp_node == NULL) || 
							(temp_node -> get_start() > curr_node -> get_stop())) {

							if (t_restart == 0) {
								break;
							}

							t_restart = 0;
							temp_node = curr_node;


						} else {

							// populate junk vars
							t_strand = temp_node -> strand;
							t_overlap = 0;

							// iterate through every subcluster 
							for (int x = 0; x < temp_node -> clust_count; x++) {

								t_start = temp_node -> clust_vec[(x * 2)];
								t_stop = temp_node -> clust_vec[(x * 2) + 1];

								if (curr_node -> check_overlap(t_start, t_stop, t_strand) == 1) {
									t_overlap = 1;
									break;
								}
							
							}

							// if overlap detected
							if (t_overlap == 1) {

								t_restart = 1;

								// iterate through every subcluster again, and add each
								for (int x = 0; x < temp_node -> clust_count; x++) {

									t_start = temp_node -> clust_vec[(x * 2)];
									t_stop = temp_node -> clust_vec[(x * 2) + 1];

									curr_node -> modify_cluster(t_start, t_stop);

								}

								// total up read counts and mark as printed
								curr_node -> read_count += temp_node -> read_count;

								// next node needs ptr to node before temp
								if ((temp_node -> next) != NULL) {
									(temp_node -> next) -> set_prev(temp_node -> prev);
								} 

								// prev node needs ptr to node after temp
								(temp_node -> prev) -> set_next(temp_node -> next);	
								inter_node = temp_node -> prev;
								
								// if temp is next node, move on
								if (temp_node == next_node) {
									next_node = temp_node -> next;
								}

								// delete node
								delete temp_node;

								// start onto next non-deleted node
								temp_node = inter_node;

							}						
						}

						// iterate
						temp_node = temp_node -> next;	
					}
				}

				// Iterate
				curr_node = next_node;

			}
		}


		///////////////////////
		// print clusters in graph
		void print_graph() {

			// Initialize accumulator and boolean
			int gene_count = 1;
			int printed;

			// initialize pointer
			Node *curr_node = head;
			Node *next_node = NULL;
			Node *temp_node = NULL;


			// check if head node exists
			if (head == NULL) {
				return;
			}


			// Iterate through nodes
			while (curr_node != NULL) {

				next_node = curr_node -> next;
				temp_node = curr_node -> next;

				// if node hasn't already been printed / added to another cluster
				if (curr_node -> printed == 0) {

					// Print cluster (line of GFF file)
					printed = curr_node -> print_cluster(contig_name, parameters, gene_count);
					curr_node -> printed = 1;

					// If cluster printed, increment gene number
					if (printed == 1) {
						gene_count ++;
					}
				
				}

				// Iterate
				curr_node = next_node;

			}

		}

};