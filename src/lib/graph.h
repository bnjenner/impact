//////////////////////////////////////
// Graph Class (really just a doubly linked list)
class Alignmnet_Graph {

	public:

		int ref; 
		std::string contig_name;
		const ImpactArguments *parameters;

		Node *head = NULL;
		Node *tail = NULL;
		Node temp;

		// graph-specific parameters
		int num_nodes = 0;
		int multimapped_reads = 0;
		int total_reads = 0;

		// Empty
		Alignmnet_Graph() {}

		// Initialize empty object
		void initialize(int ref_num, std::string ref_name, const ImpactArguments *args) {
			ref = ref_num;
			contig_name = ref_name;
			parameters = args;
		}		

		// Set Head
		int set_head(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

			while (true) {

				if (!inFile.GetNextAlignment(alignment)) { return 0; }
				if (alignment.RefID > ref) { return 0; }

				total_reads ++;

				if (alignment.IsDuplicate()) { continue; }
				if (!alignment.IsMapped()) { continue; }

				uint16_t NH_tag;
				alignment.GetTag("NH", NH_tag);
				if ((NH_tag > 1) && (parameters -> nonunique_alignments == false)) {
					multimapped_reads ++;
					continue;
				}

				// Exclude secondary alignment 
				if (!alignment.IsPrimaryAlignment() && (parameters -> nonunique_alignments == false)) {
					continue;
				} 
	
				// If paired end, check propper pair
				if (!alignment.IsProperPair() && (parameters -> library_type == "paired")) {
					continue;
				}

				if (alignment.MapQuality < parameters -> mapq) {
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

		// Create clusters of overlapping reads
		void create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

			int regions;
			int temp_start; 	// used to kill one of loops below
			int temp_strand;
			int sub_total = 1;
			std::vector<int> temp_vec = {-1, -1};

			Node *curr_node;
			tail = head;

			while (true) {

				curr_node = tail;

				if (!inFile.GetNextAlignment(alignment)) { break; }
				if (alignment.RefID > ref) { break; }

				total_reads ++;

				if (alignment.IsDuplicate()) { continue; }
				if (!alignment.IsMapped()) { continue; }

				uint16_t NH_tag;
				alignment.GetTag("NH", NH_tag);
				if ((NH_tag > 1) && (parameters -> nonunique_alignments == false)) {
					multimapped_reads ++;
					continue;
				}

				// Exclude secondary alignments
				if (!alignment.IsPrimaryAlignment() && (parameters -> nonunique_alignments == false)) {
					continue;
				}

				// If paired end, check propper pair
				if (!alignment.IsProperPair() && (parameters -> library_type == "paired")) {
					continue;
				}

				if (alignment.MapQuality < parameters -> mapq) {
		       		continue;
				}

				sub_total += 1;

				temp_vec = {alignment.Position, -1};
				temp_start = temp_vec[0];
				temp_strand = alignment.IsReverseStrand();

				curr_node -> calculate_splice(alignment, temp_vec);
				regions = temp_vec.size() / 2;

				// check if alignment represents a new node (past first subcluster)
				if ((temp_vec[0] > curr_node -> clust_vec[1]) || (temp_strand != curr_node -> strand)) {

					Node *new_node = new Node(alignment, temp_vec, ref);

					curr_node -> set_next(new_node);
					new_node -> set_prev(curr_node);
					curr_node = new_node;
					tail = curr_node;
					continue;
				}

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
							temp_start = -1;
							break;
						}	
					}
					curr_node = curr_node -> prev;
				}
			}
			return;
		}

		// print clusters in graph
		void collapse_graph() {

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

			while (curr_node != NULL) {

				next_node = curr_node -> next;
				temp_node = curr_node -> next;
				t_restart = 0;
	
				while (true) {

					if ((temp_node == NULL) || 
						(temp_node -> get_start() > curr_node -> get_stop())) {
						if (t_restart == 0) {
							break;
						}

						t_restart = 0;
						temp_node = curr_node;

					} else {

						t_strand = temp_node -> strand;
						t_overlap = 0;

						for (int x = 0; x < temp_node -> clust_count; x++) {

							t_start = temp_node -> clust_vec[(x * 2)];
							t_stop = temp_node -> clust_vec[(x * 2) + 1];

							if (curr_node -> check_overlap(t_start, t_stop, t_strand) == 1) {
								t_overlap = 1;
								break;
							}
						
						}

						if (t_overlap == 1) {

							t_restart = 1;

							for (int x = 0; x < temp_node -> clust_count; x++) {
								t_start = temp_node -> clust_vec[(x * 2)];
								t_stop = temp_node -> clust_vec[(x * 2) + 1];
								curr_node -> modify_cluster(t_start, t_stop, temp_node -> count_vec[x]);
							}

							curr_node -> read_count += temp_node -> read_count;

							if ((temp_node -> next) != NULL) {
								(temp_node -> next) -> set_prev(temp_node -> prev);
							} 

							(temp_node -> prev) -> set_next(temp_node -> next);	
							inter_node = temp_node -> prev;
							
							if (temp_node == next_node) {
								next_node = temp_node -> next;
							}

							delete temp_node;
							temp_node = inter_node;
						}						
					}

					temp_node = temp_node -> next;	
				}

				curr_node = next_node;
			}
		}

		// print clusters in graph
		void print_graph() {

			int gene_count = 1;
			int printed;

			Node *curr_node = head;

			while (curr_node != NULL) {

				printed = curr_node -> print_cluster(contig_name, parameters, gene_count);
				curr_node -> printed = 1;

				if (printed == 1) {
					gene_count ++;
				}
				
				curr_node = curr_node -> next;
			}

		}

};