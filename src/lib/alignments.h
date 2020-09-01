#include "api/BamReader.h"

using namespace BamTools;

//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

	// Attributes

		std::string file_name;			// Bam File Name 
		std::string index;			// Bam Index File 
		BamReader inFile;				// Bam File Object
		BamAlignment alignment;			// BamAlignmentRecord record;		

    	// Program options 
    	std::string library_type;		// library type
    	std::string stranded;			// strandedness of library
    	bool nonunique_alignments;		// consider secondary alignments 
    	int mapq;						// minimum mapping quality

    	// Count Statistics
    	int noncounts[5] = {0,0,0,0,0};
    		// Order:
    		  // 0 = Features
    		  // 1 = No Feature
    		  // 2 = Ambiguous
    		  // 3 = Low Qual
    		  // 4 = Unmapped
   		
    	// Contig Cache
    	std::unordered_map<int, std::string> contig_cache; 	// Unordered map 


    // Constructors
    	AlignmentFile() {}

    	AlignmentFile(const ImpactArguments *args) {

    		// Set Attributes
    		file_name = args -> alignment_file;
    		index = args -> index_file
    		library_type = args -> library_type;
    		stranded = args -> strandedness;
    		nonunique_alignments = args -> nonunique_alignments;
    		mapq = args -> mapq_min - 1;

    	}


    // Methods
    	void open() {

			if (!inFile.Open(file_name)) {
			    std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
			    throw "ERROR: Could not read alignment file.";

			if (!inFile.OpenIndex(index)) {
				std::cerr << "ERROR: Could not read index file: " << index << "\n";
				throw "ERROR: Could not read index file";
		    }


			SamHeader head = inFile.GetHeader();
			if (head.HasSortOrder()) {

				std::string sortOrder = head.SortOrder;

				if (sortOrder.compare("coordinate") != 0) {
					std::cerr << "ERROR: Sorted alignment file required.\n";
					throw "ERROR: Could not read alignment file.";
				}

			} else {
				std::cerr << "ERROR: BAM file has no @HD SO:<SortOrder> attribute. Impossible to determine sort order.\n";
					throw "ERROR: Could determine sort status.";

			}

			// Generate Ref Map (contig indicies)
			RefVector references = inFile.GetReferenceData();
			for (int i = 0; i < references.size(); i++) {

				contig_cache[i] = references.at(i).RefName;

			}

			//inFile.Rewind(); 	

		}


		// close file
		void close() {

			inFile.Close();

		}

		
	   	// head Bam
		void head() {

			int i = 0;
		 	while (inFile.GetNextAlignment(alignment) && i < 10) { // (!atEnd(bamFileIn))

		        std::cout << alignment.Name << "\n";
		        i++;

		    }  

		   	// Jump to the first entry in file.
		   	inFile.Rewind(); 	

		}


		void get_counts(AnnotationFile *annotation) {

			std::string contig;
			char strand;

			int result;

			while (inFile.GetNextAlignment(alignment)) {


				if (!alignment.IsMapped()) {
					noncounts[4]++;
					continue;
				}


				contig = contig_cache[alignment.RefID];
				
				if (library_type == "paired" && !alignment.IsProperPair()) {
					noncounts[2]++;
					continue;
				}


				if (alignment.IsReverseStrand()) {
					strand = '-';
				
				} else {
					strand = '+';
					
				}

				// Check if proper strand
				if (stranded == "forward") {

					if (library_type == "paired") {

						if ((strand == '+' && !alignment.IsFirstMate()) || 
							(strand == '-' && alignment.IsFirstMate())) {
							noncounts[2]++;
							continue;
						}
					
					}

				} else if (stranded == "reverse") {

					if (library_type == "paired") {

						if ((strand == '+' && alignment.IsFirstMate()) || 
							(strand == '-' && !alignment.IsFirstMate())) {
							noncounts[2]++;
							continue;
						}
					
					}

				}

				// Check if primary alignment
				if (!alignment.IsPrimaryAlignment() && (nonunique_alignments == false)) {
					continue;
				}

				// Check if sufficient quality
				if (alignment.MapQuality <= mapq) {
					noncounts[3]++;
		       		continue;
				}

				result = annotation -> get_feature(contig + strand, alignment.Position, alignment.GetEndPosition());						

				noncounts[result]++; 
			}

		}


};

