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

    	AnnotationFile(const ImpactArguments args) {

    		// Set Attributes
    		file_name = args.gff_file;
    		file_suffix = file_name.substr(file_name.length() - 4, 4);

    		// make lowercase
    		std::for_each(file_suffix.begin(), file_suffix.end(), [](char & c) {
		        c = ::tolower(c);
		    });

		}


	// Methods
		void open() {

			// Get # Lines 
			std::string line;
			std::ifstream gff_file (file_name);

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