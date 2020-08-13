#include <fstream>
#include <armadillo>
#include <math.h>

using namespace arma;

// Mapping Class
class MappingCounts {

	public:

	// Attributes
		std::string feature_name;
		int start;
		int stop;
		int length;

		mat counts;
        int num_counts = 0;

    // Inialize
    	MappingCounts(std::string feat_name, int interval_start, int interval_stop) {

    		// Set Attributes
    		feature_name = feat_name;
    		start = interval_start;
    		stop = interval_stop;
    		length = stop - start + 1;

    		counts.set_size(1, length);
    		counts.zeros();

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
            num_counts ++;

    	}

    	void print() {

    		std::cout << counts << std::endl;

    	}

    	void write() {

    		std::ofstream out_file;
			out_file.open(feature_name + ".txt");

			for (int i = 0; i <= counts.n_cols; i++) {

                int temp = round(counts.at(0,i));
				out_file << temp << std::endl;
			}

			out_file.close();	

    	}

    	void fit(int max_components) {

    		bool status;
    		double likelihood;
    		double max_bic;
    		double bic;
    		int best_k = 0;
    		int p;
    		mat mean;
            mat cov;

            mat data;
            mat temp;
            data.set_size(1, 0);

            double cutoff = (num_counts * 0.05);

            for (uword i = 0; i < counts.n_cols; i++) {

                if (counts[i] >= 10) {
                    mat temp(1, counts[i], fill::ones);
                    temp = temp * (i + 1);

                    data = join_rows(data, temp);
                }
            }

            if (data.n_cols > 0) {

                //double data_norm = norm(data, 2);
                //data = data / data_norm;

        		gmm_diag model;


        		for (int k = 1; k <= max_components; k++ ) {

        			status = model.learn(data, k, eucl_dist, random_spread, 0, 10, 1e-10, false);

        			if (!status) {
        				continue;
        			}
                    
        			p = (2 * k) + k - 1;

    				likelihood = model.sum_log_p(data);

    				bic = (2 * likelihood) - (p  * log(data.n_cols)); // Liklihood penalized by complexity

    				if (max_bic < bic || k == 1) {
    					max_bic = bic;
    					best_k = k;
    					mean = model.means;
                        cov = model.dcovs;
    				}

        		}

                // Check if any models ran.
                if (best_k == 0) {
                    std::cerr << "ERROR: No model could be fit to reads for gene: " << feature_name << ".\n";
                    return;
                }



                for (int i = 1; i < mean.n_cols; i++) {

                    if (mean.at(0,i) == 0) {
                        continue;
                    }

                    for (int j = 0; j < mean.n_cols; j++) {


                        if (abs(mean.at(0,i) - mean.at(0,j)) <= 250 && i != j && mean.at(0,j) != 0) {

                            if (counts.at(0, round(mean.at(0,j))) < counts.at(0, round(mean.at(0,i)))) {

                                mean.replace(mean.at(0,j), 0);

                            } else {

                                mean.replace(mean.at(0,i), 0);

                            }

                            best_k -= 1;

                        }

                    }

                }

                // std::cerr << mean << "\n";


        		if (best_k > 0) {
                    
    	    		for (int i = 0; i < mean.n_cols; i++) {

                        if (mean.at(0,i) > cutoff) {

                            int pos = round(mean.at(0,i));
                            int intensity = counts.at(0, pos);
                            double var = sqrt(cov.at(0, i) / best_k);

                            std::cerr << feature_name << "\t";
                            std::cerr << std::fixed << pos + start << "\t";
                            std::cerr << std::fixed << intensity << "\t";
                            std::cerr << std::setprecision(2) << var << "\n";


                        }
                    }

        		}
            }
    	}

};