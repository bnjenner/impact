// Run of the mill split function
std::vector<std::string> split_gff(const std::string *line, char d) {


	std::vector<std::string> cols;
	std::string col = "";

	for (char c: *line) {

		if (c == d) { 

	    	cols.push_back(col); 
	    	col = ""; 

		} else { 

	    	col = col + c;

		} 
	}  

	if (col != "") {
		cols.push_back(col); 
	}

	return cols;

}


// Parse GFF annotation for ID 
std::string get_id(const std::string *line, const std::string *file_suffix) {

	std::string id = "";
	std::string tag;
	std::string prev;
	char sep;
	int n;

	// Check for gff of gtf
	if (*file_suffix == ".gtf") {

		tag = "gene_id";
		prev = "XXXXXXX";
		sep = ' ';
		n = 6;

	} else {

		tag = "ID";
		prev = "XX";
		sep = '=';
		n = 2;

	}


	bool record = false;

	for (char c: *line) {

		if (c == sep && prev == tag) { 

			//std::cout << "got em\n";

			record = true;

		} else if (record == true && c != '\'' && c != '\"' && c != ' '){

			if (c == ';' && record == true) {
				break;
			} 
			
			id = id + c;

		} else {

			prev = prev.substr(1,n) + c; 
		
		}
	}  

	return id;

}
