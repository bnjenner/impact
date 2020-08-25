class Feature {

	public:

		std::string name;
		std::string type;
		std::string contig;
		char strand;

		int start;
		int stop;


		Feature() {}

		Feature(const std::string *line, const std::string *file_suffix) {

			// Found in lib/utils.H
			std::vector<std::string> split = split_gff(line, '\t');

			name = get_id(&split[8], file_suffix);
	    	type = split[3];
	    	contig = split[0];
	    	strand = split[6][0];

	    	// 1 to 0 based indexing
	    	start = std::stoi(split[3]) - 1;
	    	stop = std::stoi(split[4]) - 1;

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


// Annotation Class
class  AnnotationFile {

	public:

	// Attributes
		std::string file_name;					// Name of annotation file
		std::string file_suffix;					// Type of annotation file
		int total_features = 0;					// Number of features (lines)		
    	std::vector<Feature> feature_cache; 	// Unordered map 


    // Constructor
    	AnnotationFile() {}

    	AnnotationFile(const ImpactArguments *args) {

    		// set attributes
    		file_name = args -> gff_file;
    		file_suffix = file_name.substr(file_name.length() - 4, 4);

    		// make lowercase
    		std::for_each(file_suffix.begin(), file_suffix.end(), [](char & c) {
		        c = ::tolower(c);
		    });

		}


	// Methods
		void open() {

			// parse annotation file 
			std::string line;
			std::ifstream gff_file (file_name);

			if (gff_file.is_open()) {

				while ( getline(gff_file, line) ) {

					if (line[0] != '#') {

						feature_cache.emplace_back(&line, &file_suffix);
						total_features++;

					}
				}

				gff_file.close();
			}


		}

};