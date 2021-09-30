#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

//////////////////////////////////////
// Annotation Class
class AnnotationFile {

	public:

		// Useful global variables 
		std::string feature_tag;
		std::string annotation_file_name;
		std::string chrom; 

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
		    int begin;
		    int end;

		    int temp_count = 0;


			// Initialize first gene
			readRecord(record, gffIn);
			curr_id = toCString(record.tagValues[0]);
			Node *curr_node = new Node(curr_id, 
									   (record.strand == '+') ? 0 : 1,
									   toCString(record.ref));
			head = curr_node;
			tail = curr_node;


			// iterate through genes
		    while (!atEnd(gffIn)) {

		    	temp_count ++;

		        readRecord(record, gffIn);
		      
	     		temp_id = toCString(record.tagValues[0]);

		        if (temp_id != curr_id) {

		        	if (record.type == "gene") {

		        		curr_id = temp_id;

		        		// if (curr_node -> gene_id == "ENSMUSG00000028033.17") {

		        		// 	for (int i = 0; i < curr_node -> clust_count; i++) {
		        		// 		std::cerr << "\t" << curr_node -> clust_vec[(2 * i)]
		        		// 				  << "\t" << curr_node -> clust_vec[(2*i)+1] << "\n";
		        		// 	}
		        		// }

		        		// Create node
						Node *new_node = new Node(temp_id, 
												  (record.strand == '+') ? 0 : 1,
												  toCString(record.ref));

						// link nodes within graph
						curr_node -> set_next(new_node);
						new_node -> set_prev(curr_node);
						curr_node = new_node;
						tail = curr_node;		        	
		        	} 

		        } else if (record.type == feature_tag) {

		        	begin = record.beginPos;
		        	end = record.endPos;


		        	curr_node -> modify_cluster(begin, end, 0);

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