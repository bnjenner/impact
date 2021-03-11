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
		Graph graph;

		// Program options 
		std::string file_name;
		std::string index;
		Parameters parameters; 		// parameters struct (found in parser.h)

		int chr_num;
		RefVector references;
		std::unordered_map<int, std::string> contig_cache;  // Unordered map 


	////////////////////////////
	// Constructors

		// Empty
		AlignmentFile() {};

		// Initialized
		AlignmentFile(const ImpactArguments *args, int ref) {

			// Set Attributes
			file_name = args -> alignment_file;
			index = args -> index_file;
			parameters.library_type = ((args -> library_type) == "paired") ? 'p' : 's';
			parameters.stranded = ((args -> strandedness) == "forward") ? 'f' : 'r'; 
			parameters.nonunique_alignments = args -> nonunique_alignments;
			parameters.mapq = args -> mapq_min - 1;
			parameters.min_cov = args -> min_coverage - 1;

			chr_num = ref;
		}

	////////////////////////////
	// Methods

		///////////////////////
		// Open files
		void open() {

			// Open alignment file
			if (!inFile.Open(file_name)) {
				std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
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
		// Grab Alignments within Interval Using Bam Index
		void get_counts() {

			int jump = 0;

			// Create Graph object
			graph.initialize(chr_num, contig_cache[chr_num], parameters);

			// Jump to desired region in bam
			if (!inFile.Jump(chr_num, jump)) {
				std::cerr << "[ERROR: Could not jump to region: " << chr_num << ":" << jump << ".]\n";
				return;
			}

			// Set head of graph
			if (!graph.set_head(inFile, alignment)) {
				return;
			}

			// Create clusters of reads (overlap them)
			graph.create_clusters(inFile, alignment);

			// collapse overlapping clusters
			graph.collapse_graph();

		}


		///////////////////////
		// Print clusters
		void print_counts() {

			// report counts for read cluster
			graph.print_graph();
		}


		///////////////////////
		// Close files
		void launch() {

			this -> open();
			this -> get_counts();
			this -> close();
		
		}		


};

