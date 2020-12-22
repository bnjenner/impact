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

			if (temp_junct_stop == -1) {
				temp_junct_start = -1;
			}
		}


		///////////////////////
		// Create Adjancecy Matrix and Populate Start / End / Char Arrays.

		void create_adjacency(BamReader &inFile, BamAlignment &alignment, Parameters &parameters, 
							  std::string &next_id, int *dead_zone) {

			// Initialize loop variables
			int temp_end;
			int temp_strand;
			int temp_junct_start;
			int temp_junct_stop;

			while (true) {

				// If end of file is reached
				if (!inFile.GetNextAlignment(alignment)){
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

				// get strand
				temp_strand = alignment.IsReverseStrand();

				// reset splice site variables
				temp_junct_start = alignment.Position;
				temp_junct_stop = -1;

				// calculate splice sites
				calculate_splice(alignment, temp_junct_start, temp_junct_stop);
						

				// Check if new cluster hasnt been initialized
				if (start == -1) {

					// If new cluster or next id matches
					if (next_id == alignment.Name || next_id == "NA") {

						// begin new gene cluster
						start = alignment.Position;
						end = temp_end;
						strand = temp_strand;
						next_id = "READY";
						next_pos = 0;
						nodes = 1;

						// determine if splice junction is present
						if (temp_junct_stop != -1) {

							junct_start = temp_junct_start;
							junct_stop = temp_junct_stop;

						}
					}

					continue;
				}

				// Check if strand is the same as read cluster
				if (temp_strand != strand) {
					
					// if not, and no next ID is found, assign net
					if ((next_id == "READY") && 
						((alignment.Position > dead_zone[temp_strand + 2]) || (temp_end < dead_zone[temp_strand]))) {
						next_id = alignment.Name;
						next_pos = alignment.Position;

					}
				}


				// if read overlaps wit hread cluster (excluding splice junction)
				if ((alignment.Position <= end) && (temp_strand == strand) &&
					(((alignment.Position < junct_start) || (junct_start == -1)) || 
						((temp_end > junct_stop) && (temp_junct_stop < junct_stop) && (junct_stop != -1)))) {

					// add read to group
					nodes ++;

					// elongate read cluster end
					if (temp_end > end) {
						end = temp_end;
					}

					////////////////////////////////////////

					// FOLLOWING CODE NEEDS TO BE EDITTED TO ACCOUNT FOR MULTIPLE SPLICE SITES

					// reassign slice junction bounds if read is within splice bounds
					if ((junct_stop > alignment.Position) && (alignment.Position > junct_start)) {
						junct_stop = alignment.Position;
					}

					if ((temp_end > junct_start) && (junct_start != -1)) { 
						junct_start = temp_end;
					} 

					junct_start = (junct_start > temp_junct_start) ? junct_start : temp_junct_start;
					junct_stop = ((junct_stop < temp_junct_stop) && (junct_stop != -1)) ? junct_stop : temp_junct_stop;

					////////////////////////////////////////


				} else {

					// If read has no overlapp and new read is passed
					if ((nodes == 1) && (alignment.Position > end)) {
					
						// if no new cluster, assign new cluster 
						if (next_id == "READY") { 
							start = alignment.Position;
							end = temp_end;
							strand = temp_strand;
							next_id = "READY";
							next_pos = 0;

							junct_start = temp_junct_start;
							junct_stop = temp_junct_stop;

						// or go to new cluster
						} else {
							break;
						}


					} else {

						// if new no new 
						if (next_id == "READY") {
							next_id = alignment.Name;
							next_pos = alignment.Position;
						}

						// Hopefully prevents repeating read clusters
						if (alignment.Position > end) {
							
							dead_zone[strand] = (junct_stop == -1) ? start : junct_stop;
							dead_zone[strand + 2] = end;
							break;
						}
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

			if ((start == -1) || nodes == 1) {
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

