#include <seqan/basic.h>
#include <seqan/gff_io.h>

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
	AnnotationFile(const AnnotationFile &anno) {
		annotation_file_name = anno.annotation_file_name;
		feature_tag = anno.feature_tag;
		head = anno.head;
		tail = anno.tail;
	}


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

	// Print annotation graph
	void print_graph() {
		Node *curr_node = head;
		while (curr_node != NULL) {
			std::cout << curr_node -> gene_id << "\t";
			for (int i = 0; i < curr_node -> clust_count; i ++) {
				std::cout << curr_node -> clust_vec[(2 * i)] << "\t"
				          << curr_node -> clust_vec[(2 * i) + 1] << "\t";
			}
			std::cout << "\n";
			curr_node = curr_node -> next;
		}

	}

	// parse lines of annotation file
	void parse_annotation_line(std::string &line, std::vector<std::string> &columns) {

		size_t n = 0;
		std::istringstream iss(line);
		std::string column;

		while (std::getline(iss, column, '\t')) {  // but we can specify a different one
			columns.at(n) = column;
			n++;
		}
	}

	// parse 9th column of annotation files (tags)
	void parse_annotation_tags(std::string &tag_column, std::vector<std::string> &tags, bool isGFF) {

		size_t n = 0;
		std::istringstream iss(tag_column);
		std::string tag;

		// GTF
		if (isGFF) {

			std::string subtag;

			while (std::getline(iss, tag, ';')) { // but we can specify a different one

				std::istringstream iss2(tag);
				while (std::getline(iss2, subtag, '=')) {
					tags.push_back(subtag);
				}
			}

			// GFF
		} else {

			while (std::getline(iss, tag, ' ')) { // but we can specify a different onew

				if (n % 2 == 1) {
					tag.resize(tag.size() - 2);
					tag.erase(0, 1);
				}
				tags.push_back(tag);
				n++;
			}
		}


	}

	// Create graph structure
	void create_gene_graph() {

		std::ifstream infile(annotation_file_name);

		if (!infile) {
			std::cerr << "ERROR: Could not read annotation file: " << annotation_file_name << "\n";
			throw "ERROR: Could not read annotation file.";
		}

		// create col and line variables
		std::vector<std::string> tags;
		std::vector<std::string> columns{9, ""};
		std::string line;

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

			// if (line.substr(0, 1) != "#") {   // if not header line
			if (line.find_first_of('#') != 0) { 
				parse_annotation_line(line, columns);

				if (columns[2] == feature_tag) {	// if correct type

					tags.clear();
					parse_annotation_tags(columns[8], tags, isGFF);

					// get id (index is useful here)
					for (int i = 0; i < tags.size(); i++) {
						if (tags[i] == feature_id) {
							temp_id = tags[i + 1];
							break;
						}
					}

					// if no feature ID field, throw error
					if (temp_id.empty()) {
						std::cerr << "ERROR: Could not identify Feature ID.\n";
						throw "ERROR: Could not identify Feature ID.";
					}

					// if new feature
					if (temp_id.compare(curr_id) != 0) {

						curr_id = temp_id;

						if (stranded.compare("reverse") == 0) {
							temp_strand = 1 - ((columns[6] == "+") ? 0 : 1);
						} else {
							temp_strand = (columns[6] == "+") ? 0 : 1;
						}

						begin = std::stoi(columns[3]) - 1;
						end = std::stoi(columns[4]) - 1;

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
							curr_node -> set_next(new_node);
							new_node -> set_prev(curr_node);
							curr_node = new_node;
							tail = curr_node;
						}


					} else {
						begin = std::stoi(columns[3]) - 1;
						end = std::stoi(columns[4]) - 1;
						curr_node -> modify_cluster(begin, end, 0);
					}

				}
			}
		}

	}
};