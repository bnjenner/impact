

// // Iterate Through GFF
int getCounts(AnnotationFile &annotation, AlignmentFile &alignment, bool peak_detection) {
	
	Feature feat_obj;
	std::string feature_name = "";
	std::string contig;
    char strand;
    int num_alignments = 0;
    int i = 0;
    int start;
    int stop;

    int total_counts = 0;

    while (i <= annotation.total_features) {

		feat_obj = annotation.feature_cache[i];

        if (feature_name != feat_obj.name || i == annotation.total_features) {

        	if (feature_name != "" || i == annotation.total_features) {


        		MappingCounts mappedCounts(feature_name, start, stop);


        		num_alignments += alignment.findAlignments(mappedCounts, contig, start, stop, strand);


        		std::cout << feature_name << ": " << num_alignments << std::endl;

                total_counts += num_alignments;

        	}

        	num_alignments = 0;
        	feature_name = feat_obj.name;
        	contig = feat_obj.contig;
        	strand = feat_obj.strand;
        	start = feat_obj.start;
        	stop = feat_obj.stop;

        } else {

        	stop = feat_obj.stop;

        }

        i++;

    } 

    return total_counts;
}
