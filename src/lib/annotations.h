#include <functional>
#include <istream>

// Split by Delinator Template and Class
template<char deliminator>
class Deliminator: public std::string {};

// Feature Class 
class Feature {

	public:

		std::string name;
		std::string type;
		std::string contig;
		char strand;

		int start;
		int stop;

		int hash;
		std::string hash_string;

		Feature() {}

		Feature(const std::vector<std::string> *split, const std::string *file_suffix) {

			// Found in lib/utils.H
			name = get_id(&(*split)[8], file_suffix);

	    	type = (*split)[3];
	    	contig = (*split)[0];
	    	strand = (*split)[6][0];

	    	// 1 to 0 based indexing
	    	start = std::stoi((*split)[3]) - 1;
	    	stop = std::stoi((*split)[4]) - 1;

	    	hash_string = contig + strand;
	    	hash = std::hash<std::string>{}(hash_string) % 10000;

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


// Contig Bin Class 
class ContigBin {

	public:

		std::string hash_string;
		std::vector<Feature> feature_cache; 	// Unordered map 

	  // Constructor Classes
		ContigBin() {}

		ContigBin(std::string str) {
			hash_string = str;

		}

	  // Methods
		void add(Feature feat_obj) {
			feature_cache.push_back(feat_obj);
		
		}

};


// Annotation Class
class  AnnotationFile {

	public:

	// Attributes
		std::string file_name;					// Name of annotation file
		std::string file_suffix;					// Type of annotation file
		int total_features = 0;					// Number of features (lines)			
    	ContigBin contig_bin[10000];


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
		int key_exists(int hash, std::string hash_string, int increment) {

			// Check if hash is too large
			if (hash >= 10000) {
				return key_exists(hash, hash_string, -hash);
			}	

			// Check if Hash exists
			if (contig_bin[hash + increment].hash_string.empty()) {
				return increment; 
			}

			// Check if strings match
			if (contig_bin[hash].hash_string != hash_string) {
				return key_exists(hash, hash_string, increment + 1);
			
			}

			return -1;
		}

		void open() {

			// parse annotation file 
			std::string line;
			std::ifstream gff_file (file_name);
			std::string prev_name = "";
			int counter = 0;
			int gate;

			if (gff_file.is_open()) {

				while ( getline(gff_file, line) ) {

					if (line[0] != '#') {

						std::istringstream iss(line);
						std::vector<std::string> results((std::istream_iterator<Deliminator<'\t'>>(iss)),
														  std::istream_iterator<Deliminator<'\t'>>());


						Feature feat_obj(&results, &file_suffix);

						gate = key_exists(feat_obj.hash, feat_obj.hash_string, 0);

						if (gate == 0) {
							contig_bin[feat_obj.hash] = ContigBin(feat_obj.contig + feat_obj.strand);

						} else if (gate > 0) {
							contig_bin[feat_obj.hash + gate] = ContigBin(feat_obj.contig + feat_obj.strand);
						
						}

						contig_bin[feat_obj.hash].add(feat_obj); 

						total_features++;

					}

				}

				gff_file.close();
			}

			std::cerr << total_features << "\n";
			
	}

};