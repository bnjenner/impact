#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <math.h>
#include <condition_variable>
#include <mutex>
//#include <armadillo>
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "lib/parser.h"
#include "lib/model.h"
#include "lib/node.h"
#include "lib/graph.h"
#include "lib/annotation.h"
#include "lib/alignments.h"
#include "lib/queue.h"


// Mutex and Conditinoal vars
std::mutex main_mut;
std::condition_variable main_cv;
bool MAIN_THREAD = false; // found in queue.h


// Main 
int main(int argc, char const ** argv) {


    std::cerr << "[IMPACT]\n";
    auto start = std::chrono::high_resolution_clock::now(); 


    ////////////////////////////
    // Parse arguments
    ImpactArguments args;
    seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);


    // Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;
    }

    ////////////////////////////
    // Parse input files
    std::cerr << "[Parsing Input Files...]\n";

    std::cerr << "[...Annotation File...]\n";
    AnnotationFile init_annotation(&args);
    init_annotation.create_gene_graph();


    std::cerr << "[...Alignment File...]\n";
    AlignmentFile init_alignment(&args, 0);
    init_alignment.open();
    init_alignment.get_order();
    init_alignment.close();



    ////////////////////////////
    // Generate clusters

    // Number of contigs for subdividing work
    int n = init_alignment.references.size();

    // Create vector of objects for multithreading
    std::vector<AlignmentFile*> alignments;
    alignments.reserve(n);

    for (int i = 0; i < n; i++) {
        alignments.emplace_back(new AlignmentFile(&args, i));
        alignments[i] -> copy_order(init_alignment.contig_cache);
        alignments[i] -> copy_annotation(init_annotation, i);
    }

    ///////////
    // TO DO
    // - add mutex lock for writing to stderr
    ///////////

    // initialize increment and process number
    int i = 0;
    int proc = std::max(args.threads - 1, 1);


    // establish scope for queue, once it leaves scope, it will call 
    //  destructor which joins the remaining threads and modifies 
    //  MAIN_THREAD var and allows main thread to continue
    std::cerr << "[Processing Data...]\n";
    {

        // initialize dispatch queue with n threads
        thread_queue call_queue(proc); 
        do {

            // populate dispatch queue with necessary jobs
            while (i < n) {

                // dispatch job
                call_queue.dispatch([&, i]{alignments[i] -> launch();});
                i++;
        
            }

         // Wait for queue to be emptied
        } while (!call_queue.finished());
    
    }


    // lock main thread
    std::unique_lock<std::mutex> main_lock(main_mut);
    // wait for thread_queue destructor to let us go
    main_cv.wait(main_lock, []{return MAIN_THREAD;});
    // unlock thread
    main_lock.unlock(); 
    //std::cerr << "[Processing Complete!]\n";


    // // Start modeling
    // std::cerr << "[Modeling Peak Width...]\n";
    // // Create model for width
    // Model model;

    // for (int i = 0; i < n; i++) {
    //     if (alignments[i] -> model_cluster_width() > 0) {
    //         model.tot_width += alignments[i] -> model_cluster_width(); 
    //         model.tot_exp += alignments[i] -> model_cluster_exp(); 
    //         model.n += 1;
    //     }
    // }

    // model.width = model.tot_width / model.n;
    // model.exp = model.tot_exp / model.n;

    

    ////////////////////////////
    // Writing results

    int total_ambiguous = 0;
    int total_multimapping = 0;
    int total_no_feature = 0;
    int total_low_quality = 0;
    int total_unique = 0;
    int total_reads = 0;


    // Report counts (this is single threaded for order reasons)
    std::cerr << "[Writing Results...]\n";
    std::cerr << "[...Counts Data...]\n";
    for (int i = 0; i < n; i++) {
        //alignments[i] -> refine_clusters(model.width);
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

        // Overwrite file
        std::ofstream newFile;
        newFile.open(args.gtf_output);

        for (int i = 0; i < n; i++) {
            alignments[i] -> print_gtf();
        }
    }


	// The most complicated line of "get the time" I have ever seen. 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 



    ////////////////////////////
    // Say goodbye :)
    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}
