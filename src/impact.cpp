#include <iostream>
#include <chrono> 
#include <string>
#include <vector>
#include <thread>
#include <math.h>
#include <condition_variable>
#include <mutex>
#include <armadillo>
#include "api/BamReader.h" 
#include "api/BamAux.h"
#include "lib/parser.h"
#include "lib/graph.h"
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

    // Parse arguments
    ImpactArguments args;
    seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);


    // Return Error if Parsing Error
    if (res != seqan::ArgumentParser::PARSE_OK) {
            return res;
    }


    // Parse input files
    std::cerr << "[Parsing Input Files...]\n";
    AlignmentFile alignment(&args, 0);
    alignment.open();
    alignment.close(); 


    // Number of contigs for subdividing work
    int n = alignment.references.size();

    // Create vector of objects for multithreading
    std::vector<AlignmentFile*> alignments;
    alignments.reserve(n);

    for (int i = 0; i < n; i++) {
        alignments.emplace_back(new AlignmentFile(&args, i));
    }

    ////////////////////////////////////////
    // TO DO

    // - add mutex lock for writing to stderr

    ////////////////////////////////////////

    // initialize increment and process number
    int i = 0;
    int proc = std::max(args.threads - 1, 1);

    // establish scope for queue, once it leaves scope, it will call 
    //  destructor which joins the remaining threads and modifies 
    //  MAIN_THREAD var and allows main thread to continue
    std::cerr << "[Processing Alignments...]\n";
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

    std::cerr << "[Processing Complete!]\n";

    // Report counts (this is single threaded for order reasons)
    std::cerr << "[Writing Results...]\n";
    for (int i = 0; i < n; i++) {
        alignments[i] -> print_counts();
    }


	// The most complicated line of "get the time" I have ever seen. 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 


    // Say goodbye :)
    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]\n";  

    return 0;

}