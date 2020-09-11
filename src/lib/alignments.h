using namespace BamTools;

//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

	////////////////////////////
	// Attributes

		BamReader inFile;				// Bam File Object
		BamAlignment alignment;			// BamAlignmentRecord record;	

    	// Program options 
    		std::string file_name;
	    	std::string index;
	    	Parameters parameters; 		// parameters struct (found in parser.h)

	    	std::unordered_map<int, std::string> contig_cache;  // Unordered map 
   		


    ////////////////////////////
    // Constructors

    	AlignmentFile() {}

    	AlignmentFile(const ImpactArguments *args) {

    		// Set Attributes
    		file_name = args -> alignment_file;
    		index = args -> index_file;
    		parameters.library_type = args -> library_type;
    		parameters.stranded = args -> strandedness;
    		parameters.nonunique_alignments = args -> nonunique_alignments;
    		parameters.mapq = args -> mapq_min - 1;
    		parameters.min_cov = args -> min_coverage - 1;

    	}

    ////////////////////////////
    // Methods

    	///////////////////////
    	// Open files

    	void open() {

			if (!inFile.Open(file_name)) {
			    std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
			    throw "ERROR: Could not read alignment file.";
			}
			
			if (!inFile.OpenIndex(index)) {
				std::cerr << "ERROR: Could not read index file: " << index << "\n";
				throw "ERROR: Could not read index file";
		    }


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
			RefVector references = inFile.GetReferenceData();
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
		
		void get_counts(int ref, int pos) {

			// Variable accounting for overflow when group is cut off by chunking reads
		    int jump = pos;
		    int count_overflow = 0;
		    int max_overflow = 0;
		    int peak_overflow = 0;
		    int index_overflow = 0;
		    

		    // while still reading reads on desired contig
		    while (jump != 0) {

		    	// Jump to desired region in bam
		    	if (!inFile.Jump(ref, jump)) {
		    		std::cerr << "[ERROR: Could not jump to region: " << ref << ":" << jump << ".]\n";
		    		break;
		    	}
			    		
				
				// Create Graph object
				Graph graph(ref, contig_cache[ref]);

				// Create adjency matrix and get numver of aligned reads
				graph.create_adjacency(inFile, alignment, parameters);

				// get jump for next itereation
				jump = graph.get_jump();

				// Calculate degree of read node
				graph.get_degrees();

				graph.get_clusters(inFile, alignment, parameters, index_overflow, count_overflow, max_overflow, peak_overflow);

		   	} 

		}			


};

