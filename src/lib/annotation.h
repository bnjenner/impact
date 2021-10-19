#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

//////////////////////////////////////
// Annotation Class
class AnnotationFile {

	public:

		// Useful global variables 
		std::string feature_tag;
		std::string feature_id;
		std::string annotation_file_name;
		std::string chrom;
		std::string stranded;

		// node variables
		Node *head = NULL;
		Node *tail = NULL;

		// graph-specific parameters
		int num_nodes = 0;



	////////////////////////////
	// Constructors

		// Empty
		AnnotationFile() {};

		// Proper constructor
		AnnotationFile(const ImpactArguments *args) {

			annotation_file_name = args -> annotation_file;
			feature_tag = args -> feature_tag;
			feature_id = args -> feature_id;
			stranded = args -> stranded;

		}

		// Copy constructor
	    AnnotationFile(const AnnotationFile &a1) {
	    	
	    	annotation_file_name = a1.annotation_file_name;
			feature_tag = a1.feature_tag;

			head = a1.head;
			tail = a1.tail;

	    }


	////////////////////////////
	// Methods

		///////////////////////
		// Create graph structure
		void create_gene_graph() {
		    
			// Open input file.
		    GffFileIn gffIn(toCString(annotation_file_name));
		    GffRecord record;

		    // Init temp variables
		    std::string curr_id = "";
		    std::string temp_id = "";
		    int temp_strand;
		    int begin;
		    int end;

		    // Initilize pointer
		    Node *curr_node = NULL;


			// iterate through genes
		    while (!atEnd(gffIn)) {

		    	// read next gff line
		        readRecord(record, gffIn);

		        // if the correct type
		        if (record.type == feature_tag) {

		        	// get id
		     		for (int i = 0; i < length(record.tagValues); i ++) {
		     			if (record.tagNames[i] == feature_id) {
		     				temp_id = toCString(record.tagValues[i]);
		     				break;
		     			}
		     		}

		     		// if no feature ID field, throw error
		     		if (temp_id == "") {
   		   				std::cerr << "ERROR: Could not identify Feature ID.\n";
						throw "ERROR: Could not identify Feature ID.";
	   		   		}

	   		   		// if new feature
	   		   		if (temp_id != curr_id) {

	   		   			// if (curr_id == "ENSMUSG00000027710.15") {
	   		   			// 	for (int i = 0; i < curr_node -> clust_count; i++) {
	   		   			// 		std::cerr << curr_node -> clust_vec[(2*i)] << "\t"
	   		   			// 				  << curr_node -> clust_vec[(2*i)+1] << "\n";
	   		   			// 	}
	   		   			// }

	   		   			// create new feature
	   		   			curr_id = temp_id;

	   		   			// check strandedness
	   		   			if (stranded == "reverse") {
	   		   				temp_strand = 1 - ((record.strand == '+') ? 0 : 1);
	   		   			} else {
	   		   				temp_strand = (record.strand == '+') ? 0 : 1;
	   		   			}

	   		   			// Create node
						Node *new_node = new Node(temp_id, 
												  temp_strand,
												  toCString(record.ref));

						
						// if first node
						if (curr_node == NULL) {
							curr_node = new_node;
							tail = curr_node;
							head = curr_node;
						
						} else {

							// link nodes within graph
							curr_node -> set_next(new_node);
							new_node -> set_prev(curr_node);
							curr_node = new_node;
							tail = curr_node;

						}

						// get coordinates
						begin = record.beginPos;
	        			end = record.endPos - 1;

	        			// add region
	        			curr_node -> modify_cluster(begin, end, 0);	

					
	   		   		} else {

	   		   			// get coordinates
	   		   			begin = record.beginPos;
	        			end = record.endPos - 1;

	        			// add region
	        			curr_node -> modify_cluster(begin, end, 0);


	   		   		}

		        }
			}
		}

		void print_annotation() {

			Node *curr_node = head;

			while (curr_node != NULL) {

				if (curr_node -> chrom == chrom) {
					std::cout << curr_node -> gene_id << "\t"
							  << curr_node -> read_count << "\n";
				}
			
				curr_node = curr_node -> next;
			}

		}

};