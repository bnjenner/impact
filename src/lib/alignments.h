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

		RefVector references;
		std::unordered_map<int, std::string> contig_cache;  // Unordered map 


	////////////////////////////
	// Constructors

		// Empty
		AlignmentFile() {}

		// Initialized
		AlignmentFile(const ImpactArguments *args) {

			// Set Attributes
			file_name = args -> alignment_file;
			index = args -> index_file;
			parameters.library_type = ((args -> library_type) == "paired") ? 'p' : 's';
			parameters.stranded = ((args -> strandedness) == "forward") ? 'f' : 'r'; 
			parameters.nonunique_alignments = args -> nonunique_alignments;
			parameters.mapq = args -> mapq_min - 1;
			parameters.min_cov = args -> min_coverage - 1;

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
		// Close files
		void close() {
			inFile.Close();
		}


		///////////////////////
		// Grab Alignments within Interval Using Bam Index
		void get_counts(int ref) {

			// Variable accounting for group cut off and name of first read in group 
			int jump = 0;
			std::string next_id = "NA";

			// Create Graph object
			graph.initialize(ref, contig_cache[ref], parameters);

			// Jump to desired region in bam
			if (!inFile.Jump(ref, jump)) {
				std::cerr << "[ERROR: Could not jump to region: " << ref << ":" << jump << ".]\n";
				return;
			}

			// Set head of graph
			if (!graph.set_head(inFile, alignment)) {
				return;
			}

			// Create adjency matrix and get number of aligned reads
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


};

