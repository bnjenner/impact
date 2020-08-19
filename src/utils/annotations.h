// Run of the mill split function
std::vector<std::string> split_gff(std::string line, char d) {

	std::vector<std::string> cols;
	std::string col = "";

	for (char c: line) {

		if (c == d) { 

	    	cols.push_back(col); 
	    	col = ""; 

		} else { 

	    	col = col + c;

		} 
	}  

	if (col != "") {
		cols.push_back(col); 
	}

	return cols;

}


std::string get_id(std::string line, std::string file_suffix) {

	std::string id = "";
	std::string tag;
	std::string prev;
	char sep;
	int n;

	if (file_suffix == ".gtf") {

		tag = "gene_id";
		prev = "XXXXXXX";
		sep = ' ';
		n = 6;

	} else {

		tag = "ID";
		prev = "XX";
		sep = '=';
		n = 2;

	}


	bool record = false;

	for (char c: line) {

		if (c == sep && prev == tag) { 

			//std::cout << "got em\n";

			record = true;

		} else if (record == true && c != '\'' && c != '\"' && c != ' '){

			if (c == ';' && record == true) {
				break;
			} 
			
			id = id + c;

		} else {

			prev = prev.substr(1,n) + c; 
		
		}
	}  

	return id;

}

class Feature {

	public:

		std::string name;
		std::string type;
		std::string contig;
		char strand;

		int start;
		int stop;


		Feature() {}

		Feature(std::string line, std::string file_suffix) {

			std::vector<std::string> split = split_gff(line, '\t');

			name = get_id(split[8], file_suffix);
	    	type = split[3];
	    	contig = split[0];
	    	strand = split[6][0];
	    	start = std::stoi(split[3]);
	    	stop = std::stoi(split[4]);

	    	split.clear();


		}

		void print() {

			std::cout << name << "\t";
			std::cout << type << "\t";
			std::cout << contig << "\t";
			std::cout << strand << "\t";
			std::cout << start << "\t";
			std::cout << stop << "\n";

		}

};


// Alignment Class
class  AnnotationFile {

	public:

	// Attributes
		std::string file_name;					// Name of annotation file
		std::string file_suffix;					// Type of annotation file
		int total_features = 0;					// Number of features (lines)		
    	std::vector<Feature> feature_cache; 	// Unordered map 


    // Initialize
    	AnnotationFile(const ImpactArguments args) {

    		// Set Attributes
    		file_name = args.gff_file;
    		file_suffix = file_name.substr(file_name.length() - 4, 4);

    		// make lowercase
    		std::for_each(file_suffix.begin(), file_suffix.end(), [](char & c) {
		        c = ::tolower(c);
		    });

    		// Get # Lines 
			std::string line;
			std::ifstream gff_file (args.gff_file);

			if (gff_file.is_open()) {

				while ( getline(gff_file, line) ) {

					if (line.substr(0,2) != "##") {

						Feature feature_obj(line, file_suffix);
						feature_cache.push_back(feature_obj);
						total_features++;

					}
				}

				gff_file.close();
			}

		}

};