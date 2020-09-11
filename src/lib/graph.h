class Graph {

	public:

		std::vector<int> start_positions(10000, 0);
		std::vector<int> end_positions(10000, 0);
		std::vector<char> strands(10000, ' ');

		// Initialize adjacency matrix
		arma::mat adj_matrix(10000, 10000);
		adj_matrix.zeros(); 

		arma::rowvec degrees;


	// Constructors
		Graph() {}
