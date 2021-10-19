// #include <armadillo>

using namespace BamTools;

//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

	////////////////////////////
	// Attributes

		BamReader inFile;				// Bam File Object
		BamAlignment alignment;			// BamAlignmentRecord record;

		// Data Structure
		Alignmnet_Graph graph;
		AnnotationFile annotation;

		// Parameters
		std::string alignment_file_name;
		std::string index;
		const ImpactArguments *parameters; 		// parameters struct (found in parser.h)

		// Thread-specific paramters
		int chr_num;
		RefVector references;
		std::unordered_map<int, std::string> contig_cache;  // Unordered map 

		// Modeling structures
		//arma::Row<int> cluster_lens;
		//arma::Row<int> cluster_exps;

		// Read Statistics
		int total_reads = 0;
		int ambiguous_reads = 0;
		int unique_reads = 0;
		int multimapped_reads = 0;
		int unassigned_reads = 0;


	////////////////////////////
	// Constructors

		// Empty
		AlignmentFile() {};

		// Initialized
		AlignmentFile(const ImpactArguments *args, int ref) {

			// Set Attributes
			alignment_file_name = args -> alignment_file;
			index = args -> index_file;
			parameters = args;

			chr_num = ref;
		}

	////////////////////////////
	// Methods

		///////////////////////
		// Open files
		void open() {

			// Open alignment file
			if (!inFile.Open(alignment_file_name)) {
				std::cerr << "ERROR: Could not read alignment file: " << alignment_file_name << "\n";
				throw "ERROR: Could not read alignment file.";
			}
			

			// Open index file
			if (!inFile.OpenIndex(index)) {
				std::cerr << "ERROR: Could not read index file: " << index << "\n";
				throw "ERROR: Could not read index file";
			}

		}


		///////////////////////
		// Close files
		void close() {
			inFile.Close();
		}


		///////////////////////
		// parse input file for contig order and jump position
		void get_order() {

			// Get header and check if file is sorted
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
			references = inFile.GetReferenceData();
			for (int i = 0; i < references.size(); i++) {

				contig_cache[i] = references.at(i).RefName;

			}

		}


		///////////////////////
		// copy contig cache for name and position in file
		void copy_order(const std::unordered_map<int, std::string> &init_contig_cache) {
			contig_cache = init_contig_cache;
		}


		///////////////////////
		// copy annotation
		void copy_annotation(const AnnotationFile &init_annotation, int temp_chr_num) {
			annotation = init_annotation;
			annotation.chrom = contig_cache[temp_chr_num];
		}


		///////////////////////
		// Grab Alignments within Interval Using Bam Index
		void get_counts() {

			// Create Graph object
			graph.initialize(chr_num, contig_cache[chr_num], parameters);

			// Jump to desired region in bam
			if (!inFile.Jump(chr_num)) {
				std::cerr << "[ERROR: Could not jump to region: " << chr_num << ".]\n";
				return;
			}

			// Set head of graph
			if (graph.set_head(inFile, alignment)) {

				// Create clusters of reads (overlap them)
				graph.create_clusters(inFile, alignment);

				// collapse overlapping clusters
				graph.collapse_graph();
			}

			multimapped_reads = graph.multimapped_reads;
			total_reads = graph.total_reads;

		}

		///////////////////////
		// Print clusters
		void print_counts() {

			// report counts for read cluster
			graph.print_graph();
		}


		///////////////////////
		// Print gemes
		void print_genes() {

			// report counts for read cluster
			annotation.print_annotation();
		}


		// ///////////////////////
		// // return average width
		// int model_cluster_width() {
		// 	if (cluster_lens.size() > 5) {
		// 		return arma::mean(cluster_lens);
		// 	}
		// 	return 0;
		// }

		// ///////////////////////
		// // return average exp
		// int model_cluster_exp() {
		// 	if (cluster_lens.size() > 5) {
		// 		return arma::mean(cluster_exps);
		// 	}
		// 	return 0;
		// }

		// void refine_clusters(float width) {

		// 	// initialize pointer
		// 	Node *curr_node = graph.head;
		// 	int i = 0;
		// 	int x = 0;

		// 	// check if head node exists
		// 	if (curr_node != NULL) {

		// 		// Iterate through nodes
		// 		while (curr_node != NULL) {

		// 			curr_node -> filter_clusters(parameters);

		// 			curr_node = curr_node -> next;
		// 			i++;
		// 		}

		// 	}

		// }	

		///////////////////////
		// // Init vectors given number of nodes
		// void init_vectors() {

		// 	cluster_lens.zeros(100);
		// 	cluster_exps.zeros(100);

		// 	// initialize pointer
		// 	Node *curr_node = graph.head;
		// 	int i = 0;
		// 	int x = 0;

		// 	// check if head node exists
		// 	if (curr_node != NULL) {

		// 		// Iterate through nodes
		// 		while (curr_node != NULL && x < 100) {



		// 			if ((curr_node -> clust_count == 1) && 
		// 				(curr_node -> get_total_len() <= 500) &&
		// 				(curr_node -> get_total_len() >= 100)) { 

		// 				cluster_lens[x] = curr_node -> get_total_len();
		// 				cluster_exps[x] = curr_node -> read_count;

		// 				x++;

		// 			}

		// 			curr_node = curr_node -> next;
		// 			i++;
		// 		}

		// 	} else {
		// 		x = 1;
		// 	}


		// 	cluster_exps.set_size(x);
		// 	cluster_lens.set_size(x);


		// }

		///////////////////////
		// overlap genes
		void overlap_genes() {

			// if no clutsers, don't bother
			if (graph.head == NULL) {
				return;
			}

			// Graph pointers
			Node *curr_clust = graph.head;
			Node *curr_gene = annotation.head;
			Node *temp_gene = NULL;

			// Vector of Pointers for start positions in list
			std::vector<Node *> prev_gene{annotation.head, annotation.head};
			
			// Temp variables for overlap type and overlapping reads
			int overlap = 0;
			int temp_overlap = 0;
			int cmp_overlap = 0;
			int overlapping_reads = 0;
			int temp_overlapping_reads = 0;

			// Temp variables 
			int temp_strand = 0;
			bool assigned = false;

			while (curr_clust != NULL) {

				assigned = false;

				// get strand, use it to index current gene to check
				temp_strand = curr_clust -> strand;
				curr_gene = prev_gene[temp_strand];

				while (curr_gene != NULL) {	


					// if on the wrong chromosome
					if (curr_gene -> chrom != contig_cache[chr_num]) {

						// iterate through and find genes on correct contig
						while (curr_gene -> chrom != contig_cache[chr_num])  {
							curr_gene = curr_gene -> next;
							if (curr_gene == NULL){
								break;
							}
						}
							
						// iterate through and find genes on same strand
						while (curr_gene != NULL) {
							
							if (curr_gene -> strand == temp_strand) {
								break;
							}
							curr_gene = curr_gene -> next;
						}

						// Set gene starting point to current gene.
						prev_gene[temp_strand] = curr_gene;	


					// gene is on correct contig
					} else {

						// If cluster is past gene, find next gene
						if (curr_clust -> get_start() > curr_gene -> get_stop()) {

							curr_gene = curr_gene -> next;

							// find next gene on same strand
							while (curr_gene != NULL) {
								if (curr_gene -> strand == temp_strand) {
									break;
								}
								curr_gene = curr_gene -> next;
							}


						// If cluster is before gene, find next cluster
						} else if (curr_clust -> get_stop() < curr_gene -> clust_vec[2]) {
							break;

						// if possibility of overlap
						} else {

							overlap = 0;
							overlapping_reads = 0;

							// iterate through subclusters. check overlap
							for (int x = 0; x < curr_clust -> clust_count; x++) {

								temp_overlap = curr_gene -> check_genes(curr_clust -> clust_vec[(2 * x)], 
																	    curr_clust -> clust_vec[(2 * x) + 1],
																	    curr_clust -> strand);

								
								if (temp_overlap != 0) {
									overlapping_reads += curr_clust -> count_vec[x];
								}

								overlap = std::max(temp_overlap, overlap);
							}


							// Check for ambigous overlaps (following genes)
							if (overlap != 0) {

								// Check forward
								temp_gene = curr_gene -> next;

								// Temp gene isn't NULL
								if (temp_gene != NULL) {

									// iterate through following genes
									while (temp_gene -> clust_vec[2] < curr_clust -> get_stop()) {

										temp_overlapping_reads = 0;
										cmp_overlap = 0;

										// check for overlap, if there is, mark cluster as ambigous
										for (int x = 0; x < curr_clust -> clust_count; x++) {

											temp_overlap = temp_gene -> check_genes(curr_clust -> clust_vec[(2 * x)], 
																				    curr_clust -> clust_vec[(2 * x) + 1],
																				    curr_clust -> strand);

											// if overlap, add to overlapping read accumulator
											if (temp_overlap > 0) {
												temp_overlapping_reads += curr_clust -> count_vec[x];
											}

											// keep track of max overlap
											cmp_overlap = std::max(temp_overlap, cmp_overlap);

										}

										// if number of overlapping reads is equal 
										if (temp_overlapping_reads == overlapping_reads) {

											// if overlapping scores are equal, assign ambiguous
											if (cmp_overlap == overlap) {
												curr_clust -> ambiguous = 1;
												break;
											
											// if new overlap is better, replace gene assignment and continue
											} else if (cmp_overlap > overlap) {
												curr_gene = temp_gene;
												overlapping_reads = temp_overlapping_reads;
											}

										// if number of overlapping reads is greater than, replace gene assignment and continue
										} else if (temp_overlapping_reads > overlapping_reads) {
											curr_gene = temp_gene;
											overlapping_reads = temp_overlapping_reads;
										}

										// move on to next gene
										temp_gene = temp_gene -> next;

										if (temp_gene == NULL) {
											break;
										}
									}
								}	
								
								assigned = true;

								// if not ambigous, assign reads to gene
								if (curr_clust -> ambiguous == 0) {
									curr_gene -> read_count += curr_clust -> read_count;
									unique_reads += curr_clust -> read_count;
								} else {
									ambiguous_reads += curr_clust -> read_count;
									
									// std::cerr << "\n";
									// std::cerr << contig_cache[chr_num] << "\t";
									// for (int i = 0; i < curr_clust -> clust_count; i++) {
									// 	std::cerr << "\t" << curr_clust -> clust_vec[(2*i)] << "\t"
									// 			  << curr_clust -> clust_vec[(2*i)+1] << "\n";
									// }
											 
									// std::cerr << curr_gene -> gene_id << "\t" 
									// 		  << temp_gene -> gene_id << "\n";
									// std::cerr << "\n";
								}

								break;

							} else {

								// Issue is that gapped clusters do not advance cluster, yet skip many genes
								// in the process, solution is to trackback to genes

								// move on to next gene
								curr_gene = curr_gene -> next;
							
								while (curr_gene != NULL) {

									if (curr_gene -> strand == temp_strand) {
										break;
									}
									curr_gene = curr_gene -> next;
								}	

							} 
							
						}

					}

				}

				if (assigned != true) {
					unassigned_reads += curr_clust -> read_count;
				}

				curr_clust = curr_clust -> next;
			}

		}



		///////////////////////
		// Close files
		void launch() {

			this -> open();

			// find clusters
			this -> get_counts();

			// build model
			//this -> init_vectors();

			// overlap genes
			this -> overlap_genes();

			this -> close();		
		}	

};

