using namespace BamTools;

//////////////////////////////////////
// Graph Class
class Node {

	public:

	// Cluster variables

		int strand = -1;
		int read_count = 1;

		// Clusters
		std::vector<int> clust_vec{-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)
		int clust_count = 1;

		// Links
		Node *next = NULL;
		Node *prev = NULL;

		// Node Variables
		int ishead = 0;


	////////////////////////////
	// Constructors
		Node() {}

		Node(BamAlignment &alignment) {

			strand = alignment.IsReverseStrand();
			read_count = 1;

			// Start position
			clust_vec[0] = alignment.Position;

			// Temp vars
			int temp_end = alignment.GetEndPosition() - 1;
			int temp_junct_start = -1;
			int temp_junct_stop = -1;

			// Clusters
			calculate_splice(alignment, temp_junct_start, temp_junct_stop);

			if (temp_junct_start != -1) {

				clust_vec[1] = temp_junct_start;
				clust_vec.push_back(temp_junct_stop);
				clust_vec.push_back(temp_end);
				clust_count ++;

			} else {
				
				clust_vec[1] = temp_end;
			}

		}

		Node(int temp_start, int temp_stop, int temp_junct_start, int temp_junct_stop, int temp_strand) {

			strand = temp_strand;
			read_count = 1;

			// Start position
			clust_vec[0] = temp_start;

			if (temp_junct_start != -1) {

				clust_vec[1] = temp_junct_start;
				clust_vec.push_back(temp_junct_stop);
				clust_vec.push_back(temp_stop);
				clust_count ++;

			} else {
				
				clust_vec[1] = temp_stop;
			}

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
		void calculate_splice(BamAlignment &alignment, int &temp_junct_start, int &temp_junct_stop) {

			for (int i = 0; i < alignment.CigarData.size(); i++) {

				if (alignment.CigarData[i].Type == 'N') {
					temp_junct_stop = temp_junct_start + alignment.CigarData[i].Length;
					break;

				} else if (alignment.CigarData[i].Type == 'M' || alignment.CigarData[i].Type == 'D') {
					temp_junct_start += alignment.CigarData[i].Length;

				}
			
			}

			if (temp_junct_stop == -1) {
				temp_junct_start = -1;
			}
		
		}

		////////////////////////////
		// Check overlap
		int check_overlap(int &temp_start, int &temp_stop, int &temp_strand) {

			if (strand != temp_strand) {
				return 0;
			}

			// for all clusters
			for (int i = 0; i < clust_count; i++) {

				// check if beginning or end of read exists within a cluster
				if ((temp_start >= clust_vec[i * 2]) && (temp_start <= clust_vec[(i * 2) + 1])) {
					return 1;

				} else if ((temp_stop >= clust_vec[i * 2]) && (temp_stop <= clust_vec[(i * 2) + 1])) {
					return 1;
				}

			}

			return 0;
		}

		////////////////////////////
		// Check for insert 
		int check_insert(int &int_start, int &int_stop) {

			// for all clusters
			for (int i = 0; i < clust_count; i++) {

				if (((clust_vec[(i * 2) + 2] > int_stop) && (clust_vec[(i * 2) + 1] < int_start)) ||
					 ((clust_vec[(i * 2) + 2] == 0) && (int_start > clust_vec[(i * 2) + 1]))) {

					return i;
				}

			}

			return -1;
		}

		////////////////////////////
		// insert another spliced region
		void insert_splice(std::vector<int> &temp_vec, int int_start, int int_stop) {

			if (int_start > -1) { 

				temp_vec.push_back(int_start);
				temp_vec.push_back(int_stop);

				std::sort(temp_vec.begin(), temp_vec.end());

			}

		}

		////////////////////////////
		// Delete spliced region
		void delete_splice(std::vector<int> &temp_vec, int i, int int_stop) { 

		    int factor = 0; 

			// Iterate through remaining clusters
			for (int j = i + 1; j < clust_count; j++){
				
				// Determine if clusters are joined by read
				if (clust_vec[(j * 2)] < int_stop) {
					factor = j;	
				
				} else {
					break;
				
				}

			} 

			if (factor != 0) {

				int start_index = (i * 2);
				int stop_index = (2 * factor) + 1;

				temp_vec.erase(temp_vec.begin() + start_index + 1, temp_vec.begin() + stop_index);
			}

		} 

		////////////////////////////
		// check how read fits into clusters 
		void modify_cluster(int &temp_start, int &temp_end, int &temp_junct_start, int &temp_junct_stop) {

			std::vector<int> temp_vec(clust_vec);
			int ins_check = -1;
			int ins_start = -1;
			int ins_stop = -1;
			bool clust_del = false; 

			// for all clusters
			for (int i = 0; i < clust_count; i++) {

				// Precedes Cluster and Overlaps
				if ((temp_start < clust_vec[(i * 2)]) && (temp_end > clust_vec[(i * 2)])) {

					// is it a spliced alignment?
					//  (no)
					if (temp_junct_start == -1) {

						// 5' Extension
						temp_vec[(i * 2)] = temp_start;

						// 3' Extension
						//    (if read extends past original cluster)
						if (temp_end > clust_vec[(i * 2) + 1]) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							break;
						}

					// (yes) Check for 5' Insertion
					} else if ((temp_junct_start < clust_vec[(i * 2)])) {

						// Check for 5' Extension
						//    (if splice stops before cluster starts)
						if (temp_junct_stop < clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_junct_stop;

						}

						// Check for 3' Extension
						//    (if read end is past cluster end)
						if ((temp_end > clust_vec[(i * 2) + 1]) && (temp_junct_stop < clust_vec[(i * 2) + 1])) {

							// if read extends into next cluster and includes next cluster
							if ((temp_end > clust_vec[(i * 2) + 2]) && (temp_junct_stop < clust_vec[(i * 2) + 2])) {

								// Checks if clusters are joined by read
								delete_splice(temp_vec, i, temp_end);
								clust_del = true;

								// Updates end of cluster to longest value betweeen end of cluster and end of read
								temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							} else {
								temp_vec[(i * 2) + 1] = temp_end;							
							}

						}

						// A 5' insertion can only occur if it does not overlap with the present cluster
						if (clust_vec[(i * 2) - 1] < temp_start) {

							// Inserts new cluster
							ins_start = temp_start;
							ins_stop = temp_junct_start;
						
						}

					// (yes) Check for 3' Insertion
					} else if (temp_junct_stop > clust_vec[(i * 2) + 1]) {

						ins_check = check_insert(temp_junct_start, temp_end);
							
						if (ins_check != -1) {

							// Inserts new cluster
							ins_start = temp_junct_stop;
							ins_stop = temp_end;
							clust_vec = temp_vec;

						}

						// if read extends into next cluster and includes next cluster
						if ((temp_junct_start > clust_vec[(i * 2) + 2]) && (temp_start < clust_vec[(i * 2) + 2])) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;

						} 


						// Check for 5' Extension
						//    (splice starts after cluster starts)
						if (temp_junct_start > clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_start;
						}

						// Check for 3' Extension
						//    (if splice starts after cluster ends)
						if (temp_junct_start > clust_vec[(i * 2) + 1]) {
							temp_vec[(i * 2) + 1] = temp_junct_start;
						
						} 

					//  (yes) other
					} else {

						// Check for 5' Extension
						//    (if the splice starts after read cluster starts)
						if (temp_junct_start > clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_start;
						}

						// Check for 3' Extension
						//    (splice ends before read cluster and read end is before next cluser)
						if ((temp_junct_stop < clust_vec[(i * 2) + 1]) && (temp_end < clust_vec[(i * 2) + 2])) {
							temp_vec[(i * 2) + 1] =  (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
						}
						
						// if read extends into next cluster and includes next cluster
						if ((temp_end > clust_vec[(i * 2) + 2]) && (temp_junct_stop < clust_vec[(i * 2) + 2])) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;

						} 

					}
				
				// Follows cluster and overlaps
				} else if ((temp_start > clust_vec[(i * 2)]) && (temp_start < clust_vec[(i * 2) + 1])) {

					// is it a spliced alignment?
					//  (no)
					if (temp_junct_start == -1) {

						// 3' Extension
						//    (if end of read extends past cluster)
						if (temp_end > clust_vec[(i * 2) + 1]) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							break;	
						}

					//  (yes)
					} else {

						// Check for 3' Extension
						//    (splice start extends past cluster end and is before next cluster)
						//      (doesnt apply if current cluster is last)
						if ((temp_junct_start > clust_vec[(i * 2) + 1]) && 
							((temp_junct_start < clust_vec[(i * 2) + 2]) || (clust_vec[(i * 2) + 2] == 0))) {
							
							temp_vec[(i * 2) + 1] = temp_junct_start;

						//    (splice stop is before cluster end and read end is before next cluster)
						//      (doesnt apply if current cluster is last)
						} else if ((temp_junct_stop < clust_vec[(i * 2) + 1]) && 
								   ((temp_end < clust_vec[(i * 2) + 2]) || (clust_vec[(i * 2) + 2] == 0))) {

							temp_vec[(i * 2) + 1] = (clust_vec[(i * 2) + 1] > temp_end) ? clust_vec[(i * 2) + 1] : temp_end;
						} 


						// if read extends into next cluster and includes next cluster
						if ((temp_junct_start > clust_vec[(i * 2) + 2]) && (temp_start < clust_vec[(i * 2) + 2])) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_junct_start);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_junct_start) ? temp_vec[(i * 2) + 1] : temp_junct_start;
						}


						// Check for 3' Insertion
						//    (splice stop is after cluster)
						if (temp_junct_stop > clust_vec[(i * 2) + 1]) {

							ins_check = check_insert(temp_junct_stop, temp_end);

							if (ins_check != -1) {

								// Inserts new cluster
								ins_start = temp_junct_stop;
								ins_stop = temp_end;
								clust_vec = temp_vec;
							}
						
						}

					}

				}

				if (clust_del) {
					break;
				}

			}

			// Add new sections
			insert_splice(temp_vec, ins_start - 1, ins_stop - 1);			

			clust_vec = temp_vec;
			clust_count = temp_vec.size() / 2;

		}


		////////////////////////////
		// report cluster and counts
		void print_cluster(std::string &contig_name, Parameters &parameters) {
		
			char s;

			if ((clust_vec[0] == -1) || read_count < parameters.min_cov) {
				return;
			}

			// Assign strand
			if (parameters.stranded == 'f') {
				s = (strand == 1) ? '-' : '+';
			} else {
				s = (strand == 1) ? '+' : '-';
			}

			// Print name, strand, and first start
			std::cout << contig_name << "\t" << s << "\t" << clust_vec[0];

			// Print rest of starts
			for (int i = 1; i < clust_count; i++) {
				std::cout << "," << clust_vec[(i * 2)];
			}

			// Print first stop
			std::cout << "\t" << clust_vec[1];

			// Print rest of stops
			for (int i = 1; i < clust_count; i++) {
				std::cout << "," << clust_vec[(i * 2) + 1];
			}

			// Print counts
			std::cout << "\t" << read_count << "\n"; 

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

		// Node
		Node *head = NULL;
		Node *tail = NULL;
		Node temp;


	////////////////////////////
	// Constructors
		Graph() {}

	////////////////////////////
	// Methods

		// initialize
		void initialize(int ref_num, std::string ref_name, Parameters &pre_parameters) {

			// Initialize contig number and name
			ref = ref_num;
			contig_name = ref_name;
			parameters = pre_parameters;
		}		


		// Set Head
		int set_head(BamReader &inFile, BamAlignment &alignment) {

			while (true) {

				// If no alignments
				if (!inFile.GetNextAlignment(alignment)){
					std::cerr << "[ERROR: No Alignment]\n";
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
				if (alignment.MapQuality <= parameters.mapq) {
		       		continue;
				}

				break;
			}

			temp = Node(alignment);
			temp.ishead = 1;
			head = &temp;

			//std::cerr << "\tHEAD ALIGNMENT\t" << alignment.Name << "\n"; 
			return 1;

		}


		///////////////////////
		// Create clusters of overlapping reads
		int create_clusters(BamReader &inFile, BamAlignment &alignment) {

			// Initialize loop variables
			int temp_start;
			int temp_stop;
			int temp_strand;
			int temp_junct_start;
			int temp_junct_stop;

			Node *curr_node;
			tail = head;

			// std::cerr << alignment.Name << "\n";
			// // If no alignments
			// if (!inFile.GetNextAlignment(alignment)){
			// 	std::cerr << "[ERROR: No Alignment]\n";
			// 	return 0;
			// }

			// std::cerr << alignment.Name << "\n";

			while (true) {

				curr_node = tail;
				
				////////////////////////////////////////////
				// TEST BLOCK

				// // Interval
				// clust_count = 3;

				// clust_vec[0] = 10;
				// clust_vec[1] = 20;
				// clust_vec.push_back(30);
				// clust_vec.push_back(40);
				// clust_vec.push_back(50);
				// clust_vec.push_back(60);

				// // Test Read
				// temp_start = 1;
				// temp_junct_start = 8;
				// temp_junct_stop = 9;
				// temp_end = 78;

				// // Print
				// for (int x = 0; x < clust_vec.size(); x++){
				// 	std::cerr << clust_vec[x] << "\t"; 
				// }
				// std::cerr << "\n";
				// std::cerr << clust_count << "\n";

				// modify_cluster(temp_start, temp_end, temp_junct_start, temp_junct_stop);

				// for (int x = 0; x < clust_vec.size(); x++){
				// 	std::cerr << clust_vec[x] << "\t"; 
				// }
				// std::cerr << "\n";
				// std::cerr << clust_count << "\n";

				// break;

				/////////////////////////////////////////////

				// End of File Reached
				if (!inFile.GetNextAlignment(alignment)) {
					std::cerr << "[End of File Reached.]\n";
					return 0;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					return 0;
				}

				// If alignment is a duplicate
				if (alignment.IsDuplicate()) {
					continue;
				}

				// Exclude secondary alignments
				if (!alignment.IsPrimaryAlignment() && (parameters.nonunique_alignments == false)) {
					continue;
				}

				// Check if sufficient mapping quality
				if (alignment.MapQuality <= parameters.mapq) {
		       		continue;
				}

				// get alignment start
				temp_start = alignment.Position;

				// get alignment stop
				temp_stop = alignment.GetEndPosition() - 1;

				// get strand
				temp_strand = alignment.IsReverseStrand();

				// reset splice site variables
				temp_junct_start = alignment.Position;
				temp_junct_stop = -1;

				// calculate splice sites
				curr_node -> calculate_splice(alignment, temp_junct_start, temp_junct_stop);

				// check if alignment represents a new node
				if ((temp_start > curr_node -> get_stop()) || (temp_strand != curr_node -> strand)) {

					// Create node
					Node *new_node = new Node(temp_start, temp_stop, temp_junct_start, temp_junct_stop, temp_strand);

					// link nodes within graph
					curr_node -> set_next(new_node);
					new_node -> set_prev(curr_node);
					curr_node = new_node;
					tail = curr_node;

					continue;
				
				}

				// find overlapping region
				while ((curr_node != NULL) && (temp_start < curr_node -> get_stop()))  {

					// Check if alignment overlaps with previous nodes
					if (curr_node -> check_overlap(temp_start, temp_stop, temp_strand)) {
						curr_node -> modify_cluster(temp_start, temp_stop, temp_junct_start, temp_junct_stop);
						curr_node -> read_count ++;

						break;

					} else {
						curr_node = curr_node -> prev;

					}

				}

			}

			return 0;
		}


		///////////////////////
		// print clusters in graph
		void print_graph() {

			Node *curr_node = head;
			Node *next_node = NULL;

			if (head == NULL) {
				//std::cerr << "NULL HEAD\n";
				return;
			}

			while (curr_node != NULL) {

				next_node = curr_node -> next;

				if ((next_node != NULL) && 
					((curr_node -> get_stop() > next_node -> get_start()) && 
					 (curr_node -> strand == next_node -> strand))) {
					std::cerr << "\tCHECK\t" << curr_node -> get_start() << "\n";

				} else {
					curr_node -> print_cluster(contig_name, parameters); 
				}
				
				curr_node = next_node;

			}

		}
};

