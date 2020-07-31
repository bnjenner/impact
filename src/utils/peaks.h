#include <fstream>
#include <armadillo>
#include <math.h>

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
    		counts.zeros();

    		//std::cerr << "Interval: " << start << "-" << stop
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

			// std::cerr << "___\n";
			// std::cerr << counts.n_cols << "\n";
			// std::cerr << unmapped_pre.n_cols << " " << mapped.n_cols << " " << unmapped_post.n_cols << "\n";

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

    	void fit(int max_components) {

    		std::vector<double> bic_scores;
    		bool status;
    		double likelihood;
    		double min_bic = 10000000; // arbitrarily large number
    		double bic;
    		int best_k;
    		mat mean;

    		gmm_diag model;

    		for (int k = 1; k <= max_components; k++ ) {

    			status = model.learn(counts, k, eucl_dist, random_subset, 10, 5, 1e-10, false);

    			if (!status) {
    				std::cerr << counts << "\n";
    			}

				likelihood = model.sum_log_p(counts);
				bic = likelihood - (2.5 * log(length)); // L - (1/2) p log(n)

				if (min_bic > bic) {
					min_bic = bic;
					best_k = k;
					mean = model.means;
				}

    		}

    		if (best_k > 1) {

	    		std::cerr << "------\nName: " << feature_name << "\n";
	    		std::cerr << "Mean : " << mean.at(0,0);
	    		for (int i = 1; i < mean.n_cols; i++) {
	    			std::cerr << "\t" << mean.at(0,i);
				}
				std::cerr << "\n";
				std::cerr << "K: " << best_k << "\n";
	    		std::cerr << "BIC: " << min_bic << "\n";
    		}
    	}


    	//int gmmfit

};