#include <armadillo>
#include <fstream>

using namespace arma;

// Mapping Class
class MappingCounts {

	public:

	// Attributes
		seqan::CharString feature_name;
		int start;
		int stop;
		int length;

		mat counts;

    // Inialize
    	MappingCounts(seqan::CharString feat_name, int interval_start, int interval_stop) {

    		// Set Attributes
    		feature_name = feat_name;
    		start = interval_start;
    		stop = interval_stop;
    		length = stop - start + 1;

    		counts.set_size(1, length);
		} 

    // Methods
    	void addRead(int read_start, int read_stop) {

    		int begin;
    		int end;

    		if (read_start < start) {
    			begin = 0;
    			read_start = start;
    		} else {
    			begin = read_start - start;
    		}

    		if (read_stop > stop) {
    			end = 0;
    			read_stop = stop;
    		} else {
    			end = stop - read_stop;
    		}


    		mat unmapped_pre(1, begin, fill::zeros);
			mat mapped(1, (read_stop - read_start + 1), fill::ones);
			mat unmapped_post(1, end, fill::zeros);

			mat read = join_rows(join_rows(unmapped_pre, mapped), unmapped_post);

    		counts = counts + read;

    	}

    	void print() {

    		std::cout << counts << std::endl;

    	}

    	void write() {

    		std::ofstream out_file;
			out_file.open ("gl1315.NS.00574.txt");

			for (int i = 0; i <= counts.n_cols; i++) {
				out_file << counts.at(0,i) << std::endl;
			}

			out_file.close();	

    	}

    	//int gmmfit

};