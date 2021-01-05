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

		// Clusters
		std::vector<int> clust_vec{-1, -1}; // array of cluster start and stops (evens are starts, odds are ends)
		int clust_count = 1;
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

		int check_overlap(BamAlignment &alignment, int &temp_end) {

			// for all clusters
			for (int i = 0; i < clust_count; i++) {

				// check if beginning or end of read exists within a cluster
				if ((alignment.Position > clust_vec[i * 2]) && (alignment.Position < clust_vec[(i * 2) + 1])) {
					return 1;

				} else if ((temp_end > clust_vec[i * 2]) && (temp_end < clust_vec[(i * 2) + 1])) {
					return 1;
				}

			}

			return 0;
		}


		void insert_splice(std::vector<int> &temp_vec, int start_sw, int end_sw) {

			if (start_sw != -1) { 

				temp_vec.push_back(start_sw);
				temp_vec.push_back(end_sw);

				std::sort(temp_vec.begin(), temp_vec.end());

			}

		}

		void delete_splice(std::vector<int> &temp_vec, int i, int temp_end) { 

		    int factor = 0; 

			// Iterate through remaining clusters
			for (int j = i + 1; j < clust_count; j++){
				
				// Determine if clusters are joined by read
				if (clust_vec[(j * 2)] < temp_end) {
					factor = j;	
				
				} else {
					break;
				
				}

			} 

			if (factor != 0) {

				int start_index = (i * 2);
				int stop_index = (2 * factor) + 1;

				temp_vec.erase(temp_vec.begin() + start_index + 1, temp_vec.begin() + stop_index);
			}

		} 

		// BamAlignment &alignment, 
		void modify_cluster(int &temp_start, int &temp_end, int &temp_junct_start, int &temp_junct_stop) {

			std::vector<int> temp_vec(clust_vec);
			int count_buffer = 0;
			int ins_start = -1;
			int ins_stop = -1;
			bool clust_del = false; 

			// for all clusters
			for (int i = 0; i < clust_count; i++) {

				// Precedes Cluster and Overlaps
				if ((temp_start < clust_vec[(i * 2)]) && (temp_end > clust_vec[(i * 2)])) {

					// is it a spliced alignment?
					//  (no)
					if (temp_junct_start == -1) {

						// 5' Extension
						temp_vec[(i * 2)] = temp_start;

						// 3' Extension
						//    (if read extends past original cluster)
						if (temp_end > clust_vec[(i * 2) + 1]) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							break;
						}

					// (yes) Check for 5' Insertion
					} else if (temp_junct_start < clust_vec[(i * 2)]) {

						// Check for 5' Extension
						//    (if splice stops before cluster starts)
						if (temp_junct_stop < clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_junct_stop;

						}

						// Check for 3' Extension
						//    (if read end is past cluster end)
						if ((temp_end > clust_vec[(i * 2) + 1]) && (temp_junct_stop < clust_vec[(i * 2) + 1])) {

							// if read extends into next cluster and includes next cluster
							if ((temp_end > clust_vec[(i * 2) + 2]) && (temp_junct_stop < clust_vec[(i * 2) + 2])) {

								// Checks if clusters are joined by read
								delete_splice(temp_vec, i, temp_end);
								clust_del = true;

								// Updates end of cluster to longest value betweeen end of cluster and end of read
								temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							} else {
								temp_vec[(i * 2) + 1] = temp_end;							
							}

						}

						// A 5' insertion can only occur if it does not overlap with the present cluster
						if (clust_vec[(i * 2) - 1] < temp_start) {

							// Inserts new cluster
							ins_start = temp_start;
							ins_stop = temp_junct_start;
							count_buffer ++;
						
						}

					// (yes) Check for 3' Insertion
					} else if (temp_junct_stop > clust_vec[(i * 2) + 1]) {

						// (if read ends before next cluster or current cluster is the last)
						if ((clust_vec[(i * 2) + 2] > temp_end) || (clust_vec[(i * 2) + 2] == 0)) {

							// Inserts new cluster
							ins_start = temp_junct_stop;
							ins_stop = temp_end;
							count_buffer ++;

						}

						// Check for 5' Extension
						//    (splice starts after cluster starts)
						if (temp_junct_start > clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_start;
						}

						// Check for 3' Extension
						//    (if splice starts after cluster ends)
						if (temp_junct_start > clust_vec[(i * 2) + 1]) {
							temp_vec[(i * 2) + 1] = temp_junct_start;
						
						} 

					//  (yes) other
					} else {

						// Check for 5' Extension
						//    (if the splice starts after read cluster starts)
						if (temp_junct_start > clust_vec[(i * 2)]) {
							temp_vec[(i * 2)] = temp_start;
						}

						// Check for 3' Extension
						//    (splice ends before read cluster and read end is before next cluser)
						if ((temp_junct_stop < clust_vec[(i * 2) + 1]) && (temp_end < clust_vec[(i * 2) + 2])) {
							temp_vec[(i * 2) + 1] =  (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
						}
						
						// if read extends into next cluster and includes next cluster
						if ((temp_end > clust_vec[(i * 2) + 2]) && (temp_junct_stop < clust_vec[(i * 2) + 2])) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
						
						} else {
							temp_vec[(i * 2) + 1] = temp_end;							
						}

					}
				
				// Follows cluster and overlaps
				} else if ((temp_start > clust_vec[(i * 2)]) && (temp_start < clust_vec[(i * 2) + 1])) {

					// is it a spliced alignment?
					//  (no)
					if (temp_junct_start == -1) {

						// 3' Extension
						//    (if end of read extends past cluster)
						if (temp_end > clust_vec[(i * 2) + 1]) {

							// Checks if clusters are joined by read
							delete_splice(temp_vec, i, temp_end);
							clust_del = true;

							// Updates end of cluster to longest value betweeen end of cluster and end of read
							temp_vec[(i * 2) + 1] = (temp_vec[(i * 2) + 1] > temp_end) ? temp_vec[(i * 2) + 1] : temp_end;
							
							break;	
						}

					//  (yes)
					} else { 

						// Check for 3' Extension
						//    (splice start extends past cluster end and is before next cluster)
						//      (doesnt apply if current cluster is last)
						if ((temp_junct_start > clust_vec[(i * 2) + 1]) && 
							((temp_junct_start < clust_vec[(i * 2) + 2]) || (clust_vec[(i * 2) + 2] == 0))) {
							
							temp_vec[(i * 2) + 1] = temp_junct_start;

						//    (splice stop is before cluster end and read end is before next cluster)
						} else if ((temp_junct_stop < clust_vec[(i * 2) + 1]) && 
								   (temp_end < clust_vec[(i * 2) + 2])) {

							temp_vec[(i * 2) + 1] = (clust_vec[(i * 2) + 1] > temp_end) ? clust_vec[(i * 2) + 1] : temp_end;
						}


						// Check for 3' Insertion
						//    (splice stop is after cluster)
						if (temp_junct_stop > clust_vec[(i * 2) + 1]) {

							//    (read ends before next cluster and splice ends after current cluster)
							//      (doesnt apply if current cluster is last)
							if (((clust_vec[(i * 2) + 2] > temp_end) && (clust_vec[(i * 2) + 1] < temp_junct_stop)) ||
							     (clust_vec[(i * 2) + 2] == 0)) {

								// Inserts new cluster
								ins_start = temp_junct_stop;
								ins_stop = temp_end;
								clust_vec = temp_vec;
								count_buffer ++;
							}

						}

					}

				}

				if (clust_del) {
					break;
				}

				// if (i == 0) {
				// 	break;
				// }

			}

			// Add new sections
			insert_splice(temp_vec, ins_start, ins_stop);			

			clust_vec = temp_vec;
			clust_count += count_buffer;

		}


		///////////////////////
		// Create Adjancecy Matrix and Populate Start / End / Char Arrays.

		void create_adjacency(BamReader &inFile, BamAlignment &alignment, Parameters &parameters, 
							  std::string &next_id, int *dead_zone) {

			// Initialize loop variables
			int temp_start;
			int temp_end;
			int temp_strand;
			int temp_junct_start;
			int temp_junct_stop;

			while (true) {

				////////////////////////////////////////////
				// TEST BLOCK

				// Interval
				clust_count = 3;

				clust_vec[0] = 10;
				clust_vec[1] = 20;
				clust_vec.push_back(30);
				clust_vec.push_back(40);
				clust_vec.push_back(50);
				clust_vec.push_back(60);

				// Test Read
				temp_start = 2;
				temp_end = 75; 

				temp_junct_start = 8;
				temp_junct_stop = 25;


				// Print
				for (int x = 0; x < clust_vec.size(); x++){
					std::cerr << clust_vec[x] << "\t"; 
				}
				std::cerr << "\n";

				modify_cluster(temp_start, temp_end, temp_junct_start, temp_junct_stop);

				for (int x = 0; x < clust_vec.size(); x++){
					std::cerr << clust_vec[x] << "\t"; 
				}
				std::cerr << "\n";

				break;


				/////////////////////////////////////////////

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

						junct_start = temp_junct_start;
						junct_stop = temp_junct_stop;

						if (junct_start != -1) {

							clust_vec[0] = start;
							clust_vec[1] = junct_start;

							clust_vec.push_back(junct_stop);
							clust_vec.push_back(end);

							clust_count = 2;
						} else {
							clust_vec[0] = start;
							clust_vec[1] = end;
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
				if ((temp_strand == strand) && check_overlap(alignment, temp_end)) {

					// add read to group
					nodes ++;

					//modify_cluster(alignment, temp_end, temp_junct_start, temp_junct_stop);

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


							if (junct_start != -1) {
								clust_vec[0] = start;
								clust_vec[1] = junct_start;
								clust_vec.push_back(junct_stop);
								clust_vec.push_back(end);
								clust_count = 2;
							} else {
								clust_vec[0] = start;
								clust_vec[1] = end;
								clust_count = 1;
							}

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

