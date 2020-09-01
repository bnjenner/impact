#include <functional>
#include <istream>

//////////////////////////////////////
// Split by Delinator Template and Class
template<char deliminator>
class Deliminator: public std::string {};


//////////////////////////////////////
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

		int count = 0;


	// Constructors
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

	// Methods
		void update(int end) {

			if (end > stop) {
				stop = end;
			}
		
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


//////////////////////////////////////
// Contig Bin Class 
class ContigBin {

	public:

		std::string hash_string;
		std::vector<Feature> feature_cache; 
		int num_features = 0;

	// Constructors
		ContigBin() {}

		ContigBin(std::string str) {
			hash_string = str;

		}

	// Methods
		void add(Feature feat_obj) {
			feature_cache.push_back(feat_obj);
			num_features++;
		
		}

		int get_interval(int start, int stop) {

			// Result Codes:
		 	  // 0 = Feature
		 	  // 1 = No Feature
		 	  // 2 = Ambiguous

			
			// This function is essentially a binary search for item completely encompassing
			//	a reads alignment start. It also checks for neighboring overlaps

			int found = false;
			int i = 0;
			int j = num_features - 1;
			int x;

			while (found != true) {

				x = floor((i + j) / 2);

				// std::cerr << feature_cache[x].name << "\t" << i << "\t" << j
				//  		  << "\t" << x << "\n"; 
				
				// If interval is between features or 
				//	overlaping with front of feature (5' end) or
				//  interval is 0
				if ((feature_cache[x-1].stop < start && feature_cache[x].start > stop) ||
					(feature_cache[x].start > start && feature_cache[x].start <= stop) || 
					(j < i)) {
					
					return 1;
				
				} 

				// If interval matches
				if (feature_cache[x].start <= start && feature_cache[x].stop >= start ) {
					found = true;
				
				  // If interval is greater than feature
				} else if (feature_cache[x].stop < start) {
					i = x + 1;
				
				  // If interval is less than feature
				} else if (feature_cache[x].start > stop) {
					j = x - 1;
				
				}

			}
			
			// Check for Next Overlap
			if ((feature_cache[x+1].start <= stop && feature_cache[x].name != feature_cache[x+1].name) && 
				(x != num_features - 1)) {
				
				return 2;
			} 

			// Check for Previous Overlap
			if ((feature_cache[x-1].stop >= start && feature_cache[x].name != feature_cache[x-1].name) && 
			    (x != 0)) {

				return 2;
			}

			feature_cache[x].count++;
			return 0;

		}

};


//////////////////////////////////////
// Annotation Class
class  AnnotationFile {

	public:

	// Attributes
		std::string file_name;					// Name of annotation file
		std::string file_suffix;				// Type of annotation file
		int total_features = 0;					// Number of features (lines)			
    	ContigBin contig_bin[100000];			// Bin for contigs
    	std::vector<int> order;


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
			if (hash >= 99999) {
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

			return 100000;
		}


		int get_feature(std::string hash_string, int start, int stop) {

			// Result Codes:
		 	  // 0 = Feature
		 	  // 1 = No Feature
		 	  // 2 = Ambiguous

			int hash = std::hash<std::string>{}(hash_string) % 100000;
			int gate = key_exists(hash, hash_string, 0);

			if (gate == 0) {
				return 1;
			} 

			while (contig_bin[hash].hash_string != hash_string) {

				if (contig_bin[hash].hash_string.empty()) {
					return 1;
				}

				hash++;

				if (hash >= 99999) {
					hash = 0;
				}
			
			}

			return contig_bin[hash].get_interval(start, stop);

		}

		void open() {

			// open annotation file 
			std::string line;
			std::ifstream gff_file(file_name);

			// Gate Serves as hash offset in case of collisions
			int gate; 
			int final_hash;
			bool push;
			int bin_num_features;

			// Iterate through annotation file
			if (gff_file.is_open()) {

				while (getline(gff_file, line)) {

					// Ignore header lines
					if (line[0] != '#') {

						// Split string by tabs
						std::istringstream iss(line);
						std::vector<std::string> results((std::istream_iterator<Deliminator<'\t'>>(iss)),
														  std::istream_iterator<Deliminator<'\t'>>());

						// Create Feature object using split line, get hash of object in contig_bin
						Feature feat_obj(&results, &file_suffix);
						final_hash = feat_obj.hash;

						// Check if key exists and add offset 
						gate = key_exists(final_hash, feat_obj.hash_string, 0);


						// If key doesn't exist, create new ContigBin object
						if (gate == 0) {
							contig_bin[final_hash] = ContigBin(feat_obj.contig + feat_obj.strand);

						  // If Key has offset 
						} else if (gate != 0 && gate != 100000) {
							final_hash = feat_obj.hash + gate;
							contig_bin[final_hash] = ContigBin(feat_obj.contig + feat_obj.strand);
						
						}


						push = true;
						bin_num_features = contig_bin[final_hash].num_features;


						// If object isnt first and key exists (gate == 10000)
						if (bin_num_features != 0 && gate == 100000) {

							if (contig_bin[final_hash].feature_cache[bin_num_features - 1].name == feat_obj.name) {
								contig_bin[final_hash].feature_cache[bin_num_features - 1].update(feat_obj.stop);
								push = false;
							
							} 

						} 

						if (push == true) {
							order.push_back(final_hash);
							contig_bin[final_hash].add(feat_obj); 
							total_features++;
						
						}

					}

				}

				gff_file.close();
			}

		}

};