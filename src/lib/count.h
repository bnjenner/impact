// // Iterate Through GFF
int getCounts(AnnotationFile *annotation, AlignmentFile *alignment, bool peak_detection) {
	
    int total_features = annotation -> total_features;
	Feature feat_obj;
    bool first = true;
    int prev_index = 0;
    int prev_stop;
    char prev_strand;
    int next_start;
    char next_strand;

	std::string feature_name = "";
	std::string contig;
    char strand;
    int num_alignments = 0;
    int i = 0;
    int start;
    int stop;

    int total_counts = 0;


    while (i <= total_features) {

		feat_obj = annotation -> feature_cache[i];

        if (feature_name != feat_obj.name || i == total_features) {

        	if (feature_name != "" || i == total_features) {

                MappingCounts mappedCounts(feature_name, start, stop);

               // Overlap cases
                // Previous 
                if (first == true) {
                    prev_stop = -1;
                    prev_strand = '=';
                    first = false;

                } else {
                    prev_stop = annotation -> feature_cache[i - (prev_index + 1)].stop + 1;
                    prev_strand =  annotation -> feature_cache[i - (prev_index + 1)].strand;
                }

                // Next 
                if (i == total_features) {
                    next_start = -1;
                    next_strand = '=';

                } else {
                    next_start = annotation -> feature_cache[i].start;
                    next_strand =  annotation -> feature_cache[i].strand;
                }


                // // Find alignments
                num_alignments += alignment -> findAlignments(mappedCounts, contig, start, stop, strand, peak_detection,
                                                              prev_stop, prev_strand, next_start, next_strand);

               
                // Print counts
                std::cout << feature_name << ": " << num_alignments << std::endl;

                total_counts += num_alignments;
                num_alignments = 0;
                prev_index = 0;
                
        	}

        	feature_name = feat_obj.name;
        	contig = feat_obj.contig;
        	strand = feat_obj.strand;
        	start = feat_obj.start;
        	stop = feat_obj.stop;

        } else {

        	if (feat_obj.stop > stop) {
	        	stop = feat_obj.stop;
	        }

        }

        prev_index++;
        i++;

    } 

    return total_counts;
}
