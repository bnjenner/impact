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
		bool isGFF;

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
			isGFF = args -> isGFF;

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
		// Print genes and counts
		void print_counts() {

			Node *curr_node = head;

			while (curr_node != NULL) {

				if (curr_node -> chrom == chrom) {
					std::cout << curr_node -> gene_id << "\t"
							  << curr_node -> read_count << "\n";
				}
			
				curr_node = curr_node -> next;
			}
		}


		///////////////////////
		// Print annotation graph
		void print_graph() {

			Node *curr_node = head;

			while (curr_node != NULL) {

				std::cout << curr_node -> gene_id << "\t";
				for (int i = 0; i < curr_node -> clust_count; i ++) {
					std::cout << curr_node -> clust_vec[(2*i)] << "\t" 
						  << curr_node -> clust_vec[(2*i) + 1] << "\t";

				}
				std::cout << "\n";
			
				curr_node = curr_node -> next;
			}

		}


		///////////////////////
		// parse lines of annotation file
		void parse_annotation_line(std::string &line, std::vector<std::string> &columns) {

			// create string stream and temp variable
			std::istringstream iss(line);
			std::string column;

			int n = 0;

			// parse string 
		    while (std::getline(iss, column, '\t')) {  // but we can specify a different one
		        columns.at(n) = column;
		        n ++;
		    }
		}


		///////////////////////
		// parse 9th column of annotation files (tags)
		void parse_annotation_tags(std::string &tag_column, std::vector<std::string> &tags, bool isGFF) {

			std::istringstream iss(tag_column);
			std::string tag;

			int n = 0;

			if (isGFF) {
				
				std::string subtag;

				// parse tags
			    while(std::getline(iss, tag, ';')) {  // but we can specify a different one

			    	std::istringstream iss2(tag);

			    	while (std::getline(iss2, subtag, '=')) {
				        tags.push_back(subtag);
			    	}

			    }

			} else {

				// parse tags
			    while(std::getline(iss, tag, ' ')) {  // but we can specify a different one
			        
			        // n is even
			    	if (n % 2 == 1) {
			    		tag.resize(tag.size() - 2);
			    		tag.erase(0,1);
			    	}

			        tags.push_back(tag);
			        n ++;
			    }
			}

			
		}


		///////////////////////
		// Create graph structure
		void create_gene_graph() {

			// open annotation file
			std::ifstream infile(annotation_file_name);

			// create col and line variables
			std::vector<std::string> columns{9, ""};
			std::string line;

			int i = 0;

			 // Init temp variables
		    std::string curr_id = "";
		    std::string temp_id = "";
		    int temp_strand;
		    int begin;
		    int end;

		    // Initilize pointer
		    Node *curr_node = NULL;

			// Iterate through lines in file
			while (std::getline(infile, line)) {

				// if not header line
				if (line.substr(0,1) != "#") {

					// parse line
					parse_annotation_line(line, columns);

					// if correct type
		       		if (columns[2] == feature_tag) {

		       			// parse tags
						std::vector<std::string> tags;
						parse_annotation_tags(columns[8], tags, isGFF);

						// get id
			     		for (int i = 0; i < tags.size(); i++) {
			     			if (tags[i] == feature_id) {
			     				temp_id = tags[i + 1];
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


		   		   			// create new feature
		   		   			curr_id = temp_id;

		   		   			// check strandedness
		   		   			if (stranded == "reverse") {
		   		   				temp_strand = 1 - ((columns[6] == "+") ? 0 : 1);
		   		   			} else {
		   		   				temp_strand = (columns[6] == "+") ? 0 : 1;
		   		   			}

		   		   			// get coordinates
							begin = std::stoi(columns[3]) - 1;
		        			end = std::stoi(columns[4]) - 1;

		   		   			// Create node
							Node *new_node = new Node(temp_id, 
													  temp_strand,
													  begin, end,
													  columns[0]);

							
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

						
		   		   		} else {

		   		   			// get coordinates
		   		   			begin = std::stoi(columns[3]) - 1;
		        			end = std::stoi(columns[4]) - 1;

		        			// add region
		        			curr_node -> modify_cluster(begin, end, 0);


		   		   		}

				 	}
				}
			}

		}

};