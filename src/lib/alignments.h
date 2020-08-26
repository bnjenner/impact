#include "api/BamReader.h"

using namespace BamTools;

//////////////////////////////////////
// Alignment Class
class AlignmentFile {

	public:

	// Attributes

		std::string file_name;			// Bam File Name 
		BamReader inFile;				// Bam File Objec
		BamAlignment alignment;			// BamAlignmentRecord record;		

    	// Program options 
    	std::string library_type;		// library type
    	bool nonunique_alignments;		// consider secondary alignments 
    	int mapq;						// minimum mapping quality

    	// Count Statistics
    	int noncounts[4] = {0,0,0,0};
    		// Order:
    		  // 0 = Features
    		  // 1 = No Feature
    		  // 2 = Ambiguous
    		  // 3 = Low Qual
   		
    	// Contig Cache
    	std::unordered_map<int, std::string> contig_cache; 	// Unordered map 


    // Constructors
    	AlignmentFile() {}

    	AlignmentFile(const ImpactArguments *args) {

    		// Set Attributes
    		file_name = args -> alignment_file;
    		library_type = args -> library_type;
    		nonunique_alignments = args -> nonunique_alignments;
    		mapq = args -> mapq_min - 1;

    	}


    // Methods
    	void open() {

			if (!inFile.Open(file_name)) {
			    std::cerr << "ERROR: Could not read alignment file: " << file_name << "\n";
			    throw "ERROR: Could not read alignment file.";
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

				//std::cerr << alignment.Name << "\t";

				contig = contig_cache[alignment.RefID];
				
				if (alignment.IsReverseStrand()) {
					strand = '-';
				
				} else {
					strand = '+';
					
				}

				// Check if primary alignment
				if (!alignment.IsPrimaryAlignment() && (nonunique_alignments == false)) {
					//std::cerr << "secondary\n";
					continue;
				}

				// Check if sufficient quality
				if (alignment.MapQuality <= mapq) {
					noncounts[3]++;
		       		continue;
				}

				//if (alignment.Name == "A00351:130:HHGJFDSXX:3:1623:3703:36620_GTGTAG") {
				result = annotation -> get_feature(contig + strand, alignment.Position, alignment.GetEndPosition());				
				//}

				//std::cerr << result << "\n";
				

				noncounts[result]++; 
			}

		}


};

