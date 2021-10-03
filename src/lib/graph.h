using namespace BamTools;

//////////////////////////////////////
// Graph Class (really just a doubly linked list)
class Alignmnet_Graph {

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

		// graph-specific parameters
		int num_nodes = 0;
		int multimapped_reads = 0;
		int total_reads = 0;

	////////////////////////////
	// Constructors

		// Empty
		Alignmnet_Graph() {}


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
				if (!inFile.GetNextAlignment(alignment)) { 
					return 0;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					return 0;
				}

				total_reads ++;

				// If alignment is a duplicate
				if (alignment.IsDuplicate()) {
					continue;
				}

				// If alignment is mapped
				if (!alignment.IsMapped()) {
					continue;
				}

				uint16_t NH_tag;
				alignment.GetTag("NH", NH_tag);
				if (NH_tag > 1) {
					multimapped_reads ++;
					continue;
				}

				// Exclude secondary alignment 
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

				break;
			}

			// Add alignment to head node
			temp = Node(alignment, ref);

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

			// We have a head
			int sub_total = 1;

			while (true) {

				// Start at last node
				curr_node = tail;

				// End of File Reached
				if (!inFile.GetNextAlignment(alignment)) {
					break;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					break;
				}

				total_reads ++;

				// If alignment is a duplicate
				if (alignment.IsDuplicate()) {
					continue;
				}

				// If alignment is mapped
				if (!alignment.IsMapped()) {
					continue;
				}

				uint16_t NH_tag;
				alignment.GetTag("NH", NH_tag);
				if (NH_tag > 1) {
					multimapped_reads ++;
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

				sub_total += 1;

				// write to temp vector
				temp_vec = {alignment.Position, -1};

				// get alignment start
				temp_start = temp_vec[0];

				// get strand
				temp_strand = alignment.IsReverseStrand();

				// calculate splice sites
				curr_node -> calculate_splice(alignment, temp_vec);

				// number of aligned regions
				regions = temp_vec.size() / 2;

				// check if alignment represents a new node (past first subcluster)
				if ((temp_vec[0] > curr_node -> clust_vec[1]) || (temp_strand != curr_node -> strand)) {

					// Create node
					Node *new_node = new Node(temp_vec, temp_strand, ref);

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
				while ((curr_node != NULL) && (temp_start != -1))  {

					for (int x = 0; x < regions; x++) {

						// Check if alignment overlaps with previous nodes
						if (curr_node -> check_overlap(temp_vec[(2 * x)], temp_vec[(2 * x) + 1], temp_strand)) {

							// add all clusters to vector
							for (int y = 0; y < regions; y++) {
								curr_node -> modify_cluster(temp_vec[(2 * y)], temp_vec[(2 * y) + 1], 1);
							}

							curr_node -> read_count++;
							
							// kill the loop
							temp_start = -1;
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

			// Iterate through nodes
			while (curr_node != NULL) {

				next_node = curr_node -> next;
				temp_node = curr_node -> next;

				t_restart = 0;

	
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

								curr_node -> modify_cluster(t_start, t_stop, temp_node -> count_vec[x]);

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

			// Iterate through nodes
			while (curr_node != NULL) {

				// Print cluster (line of GFF file)
				printed = curr_node -> print_cluster(contig_name, parameters, gene_count);
				curr_node -> printed = 1;

				// If cluster printed, increment gene number
				if (printed == 1) {
					gene_count ++;
				}
				
				// Iterate
				curr_node = curr_node -> next;

			}

		}

};