using namespace BamTools;

//////////////////////////////////////
// Graph Class
class Graph {

	public:

	////////////////////////////
	// Attributes

		// Cluster variables
		int start = -1;
		int end = -1;
		int strand = -1;
		int nodes = 1;
		int junct_start = -1; // Will be implemented later to account for splicing
		int junct_stop = -1;

		// Jump position
		int next_pos = -1;

		// Useful global variables
		int ref; 
		std::string contig_name;


	////////////////////////////
	// Constructors

		Graph(int ref_num, std::string ref_name) {

			// Initialize contig number and name
			ref = ref_num;
			contig_name = ref_name;

		}


	////////////////////////////
	// Methods


		void reset() {

			start = -1;
			end = -1;
			strand = -1;
			nodes = 0;
			next_pos = -1;
			junct_start = -1; // Will be implemented later to account for splicing
			junct_stop = -1;

		}

		void calculate_splice(BamAlignment &alignment, int &temp_junct_start, int &temp_junct_stop) {

			for (int i = 0; i < alignment.CigarData.size(); i++) {

				if (alignment.CigarData[i].Type == 'N') {
					temp_junct_stop = temp_junct_start + alignment.CigarData[i].Length;
					break;

				} else if (alignment.CigarData[i].Type == 'M' || alignment.CigarData[i].Type == 'D') {
					temp_junct_start += alignment.CigarData[i].Length;

				}
			
			}
		}


		///////////////////////
		// Create Adjancecy Matrix and Populate Start / End / Char Arrays.

		void create_adjacency(BamReader &inFile, BamAlignment &alignment, Parameters &parameters, 
							  std::string &next_id, int *dead_end) {

			// Initialize loop variables
			int temp_end;
			int temp_strand;
			int temp_junct_start;
			int temp_junct_stop;

			while (true) {

				// If end of file is reached
				if (!inFile.GetNextAlignment(alignment)){
					next_pos = 0;
					break;
				}

				// If next chromosome is reached, get out of town.
				if (alignment.RefID > ref) {
					next_pos = 0;
					break;
				}

				if (alignment.IsDuplicate()) {
					continue;
				}

				// Exclude secondary alignments
				if (!alignment.IsPrimaryAlignment() && (parameters.nonunique_alignments == false)) {
					continue;
				}

				// Check if sufficient mapping quality
				if (alignment.MapQuality <= parameters.mapq) {
		       		continue;
				}

				// get alignment end 
				temp_end = alignment.GetEndPosition() - 1;

				// // If alignment spans double its actual length in bases, ignore it, its trash
				// if ((50000) <= (temp_end - alignment.Position)) {
				// 	continue;
				// }

				// get strand
				temp_strand = alignment.IsReverseStrand();


				temp_junct_start = alignment.Position;
				temp_junct_stop = -1;

				calculate_splice(alignment, temp_junct_start, temp_junct_stop);

				// Check if cluster has been initialized
				if (start == -1) {

					// If new cluster or next id matches
					if (next_id == alignment.Name || next_id == "NA") {
						start = alignment.Position;
						end = temp_end;
						strand = temp_strand;
						next_id = "READY";
						nodes = 1;

						if (temp_junct_stop != -1) {

							junct_start = temp_junct_start;
							junct_stop = temp_junct_stop;

						}

					}

					continue;
				}

				// Check strand
				if (temp_strand != strand) {

					if (alignment.Name == "A01102:83:HCWVFDSXY:2:2321:17924:30686") {
						std::cerr << next_id << "\t" 
								  << alignment.Position << "\t" << dead_end[temp_strand] << "\n";
					}
					
					if ((next_id == "READY") && (alignment.Position > dead_end[temp_strand])) {
						next_id = alignment.Name;
						next_pos = alignment.Position;
					}

					continue;
				}


				// if reads overlap
				if ((alignment.Position <= end) && 
					(((alignment.Position < junct_start) || (junct_start == -1)) || 
						((temp_end > junct_stop) && (temp_junct_stop < junct_stop) && (junct_stop != -1)))) {

					if (start == 25892048) {
						std::cerr << alignment.Name << "\t" << alignment.Position << "\t" << temp_end << "\n"
								  << "\t" << junct_start << "\t" << junct_stop << "\n";
					}

					nodes ++;

					if (temp_end > end) {
						end = temp_end;

					}

					if (junct_stop > alignment.Position) {
						junct_stop = alignment.Position;
					}

					if (temp_junct_stop != -1) {

						junct_start = (junct_start > temp_junct_start) ? junct_start : temp_junct_start;
						junct_stop = (junct_stop < temp_junct_stop && junct_stop != -1) ? junct_stop : temp_junct_stop;

					} else if ((end > junct_start) && (junct_start != -1)) { 
						junct_start = temp_end;

					} 

				} else {

					if (nodes == 1) {
						start = alignment.Position;
						end = temp_end;
						strand = temp_strand;

						if (temp_junct_stop != -1) {

							junct_start = temp_junct_start;
							junct_stop = temp_junct_stop;

						} else {
							junct_start = -1;
							junct_stop = -1;
						}

						continue;

					} else {

					// Hopefully prevents repeating read clusters
						dead_end[strand] = end;

						if (next_id == "READY") {
							next_id = alignment.Name;
							next_pos = alignment.Position;
						}
						break;
					}

				}

			}

		}

		///////////////////////
		// gets next position to jump to in next iteration

		int get_jump() {
			return next_pos;
		}

		///////////////////////
		// report count of recent cluster

		void print_counts(char stranded) {
		
			char s;

			if (start == -1) {
				return;
			}

			if (stranded == 'f') {
				s = (strand == 1) ? '-' : '+';
			} else {
				s = (strand == 1) ? '+' : '-';
			}

			/// SOMETHING IN HERE REPORTS A JUNCT START > JUNCT STOP, NEEDS A FIX

			std::cout << contig_name << "\t" << s << "\t" << start + 1 << "\t" 
					  << end + 1 << "\t" << junct_start << "\t" 
					  << junct_stop << "\t" << nodes << "\n"; 

		
			// if (junct_start != -1 && junct_stop != -1) {
			// 	std::cerr << junct_start << "\t" << junct_stop << "\n";
			// }
		}

};

