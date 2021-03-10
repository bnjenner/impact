#include <functional>
#include <cstdint>
#include <queue>

// Special thanks to https://github.com/embeddedartistry/embedded-resources/blob/master/examples/cpp/dispatch.cpp 
//	for the code inspiration :)
extern bool MAIN_THREAD;

class thread_queue {

	typedef std::function<void(void)> call;

	public:

	////////////////////////////
	// Attributes

		std::mutex mlock;
		std::vector<std::thread> threads;
		std::queue<call> call_queue;
		std::condition_variable cv;
		bool quit = false;

	////////////////////////////
	// Constructors

		thread_queue(size_t thread_cnt) : threads(thread_cnt) {
			
			for (size_t i = 0; i < threads.size(); i++) {
				threads[i] = std::thread(&thread_queue::thread_handler, this);
			}

		}

	////////////////////////////
	// Destructors

		~thread_queue() {

			// creat lock and set quit to true
			std::unique_lock<std::mutex> lock(mlock);
			quit = true;

			// unlock, notify all threads
			lock.unlock();
			cv.notify_all();


			// wait for threads to finish
			for (size_t i = 0; i < threads.size(); i++) {
				
				if (threads[i].joinable()) {
					threads[i].join();
				}

			}

			::MAIN_THREAD = true;
		}

	////////////////////////////
	// Methods

		////////////////////////////
		// dispatch (basically enqueue)
		void dispatch(const call &job) {
	
			// create lock
			std::unique_lock<std::mutex> lock(mlock);
			// enqueue
			call_queue.push(std::move(job));
			// unlock
			lock.unlock();
			// notify conditional var
			cv.notify_one();
		}


		////////////////////////////
		// thread handler
		void thread_handler(void) {

			// create lock
			std::unique_lock<std::mutex> lock(mlock);

			do {

				// wait for data (still need to think about this one)
				cv.wait(lock, [this]{
					return (call_queue.size() || quit);
				});

				// after wait, we own the lock
				if (!quit && call_queue.size()) {

					// create job
					auto job = std::move(call_queue.front());
					call_queue.pop();

					//unlock now that we're done messing with the queue
					lock.unlock();

					////////////////////////////////////////

					// HERE IS THE ERROR, SOMETHING WITH HOW IT IS CALLED IDK
					// call job
					//job();

					///////////////////////////////////////

					std::cerr << "calling\n";

					// lock 
					lock.lock();
				}

			} while (!quit);
		
		}

		bool finished(){
			return (call_queue.empty());
		}

		int size(){
			return call_queue.size();
		}
};
