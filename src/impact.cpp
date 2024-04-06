#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <math.h>
#include <condition_variable>
#include <mutex>
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "lib/parser.h"
#include "lib/node.h"
#include "lib/graph.h"
#include "lib/annotation.h"
#include "lib/alignments.h"
#include "lib/queue.h"

// Mutex and Conditinoal vars
std::mutex main_mut;
std::condition_variable main_cv;
bool MAIN_THREAD = false;

//////////////////////////////////////
// Main 
int main(int argc, char const ** argv) {

    // initialize clock
    auto start = std::chrono::high_resolution_clock::now(); 

    std::cerr << "[IMPACT]\n";

    // Parse arguments
    ImpactArguments args;
    seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);

    // Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;
    }

    
    std::cerr << "[Parsing Input Files...]\n";

    // Parse annoatation file
    std::cerr << "[...Annotation File...]\n";
    AnnotationFile init_annotation(&args);
    init_annotation.create_gene_graph();

    // Parse alignment file    
    std::cerr << "[...Alignment File...]\n";
    AlignmentFile init_alignment(&args, 0);
    init_alignment.open();
    init_alignment.get_order();
    init_alignment.close();


    // Number of contigs for subdividing work
    //  across multiple threads
    int n = init_alignment.references.size();

    // Create vector of pointers for multithreading
    std::vector<AlignmentFile*> alignments;
    alignments.reserve(n);

    // load 
    for (int i = 0; i < n; i++) {
        alignments.emplace_back(new AlignmentFile(&args, i));
        alignments[i] -> copy_order(init_alignment.contig_cache);
        alignments[i] -> copy_annotation(init_annotation, i);
    }


    /*
    TO DO:
        - add mutext lock for writing to stdeff
    */

    std::cerr << "[Processing Data...]\n";

    int i = 0;
    const int proc = std::max(args.threads - 1, 1);
    {
        thread_queue call_queue(proc);  // initialize dispatch queue with n threads
        do {
            // populate dispatch queue with necessary jobs
            while (i < n) {
                call_queue.dispatch([&, i]{alignments[i] -> launch();}); // dispatch job
                i++;
            }

        // Wait for queue to be emptied
        } while (!call_queue.finished());
    }


    std::unique_lock<std::mutex> main_lock(main_mut);   // lock main thread
    main_cv.wait(main_lock, []{return MAIN_THREAD;});   // wait for thread_queue destructor to let us go
    main_lock.unlock();                                 // unlock thread

    // replace with struct for christ sake
    int total_ambiguous = 0;
    int total_multimapping = 0;
    int total_no_feature = 0;
    int total_low_quality = 0;
    int total_unique = 0;
    int total_reads = 0;


    std::cerr << "[Writing Results...]\n";
    std::cerr << "[...Counts Data...]\n";

    for (int i = 0; i < n; i++) {
        alignments[i] -> print_genes();
        total_ambiguous += alignments[i] -> ambiguous_reads; 
        total_unique += alignments[i] -> unique_reads; 
        total_multimapping += alignments[i] -> multimapped_reads;
        total_no_feature += alignments[i] -> unassigned_reads;
        total_reads += alignments[i] -> total_reads;     
    }
    
    std::cout << "__unique\t" << total_unique << "\n";
    std::cout << "__ambiguous\t" << total_ambiguous << "\n";
    std::cout << "__multimapping\t" << total_multimapping << "\n";
    std::cout << "__unassigned\t" << total_no_feature << "\n";
    std::cout << "__total\t" << total_reads << "\n";


    // Output read cluster gtf if specified
    if (args.gtf_output != "") {
        std::cerr << "[...Output GTFs...]\n";
        std::ofstream newFile;
        newFile.open(args.gtf_output);
        for (const auto &align: alignments) {
            align -> print_gtf();
        }
    }


	// The most complicated line of "get the time" I have ever seen. 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 

    // Say goodbye :)
    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;
}
