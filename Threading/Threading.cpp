// condition_variable example
#include <iostream>           // std::cout
#include <thread>             // std::thread
#include <mutex>              // std::mutex, std::unique_lock
#include <condition_variable> // std::condition_variable

const int NumberOfThreads = 4;

std::mutex mutThread[NumberOfThreads];
std::condition_variable cvsThread[NumberOfThreads];

std::mutex mutMain[NumberOfThreads];
std::condition_variable cvsMain[NumberOfThreads];

int isMainReady = 0;
int areThreadsReady = 0;

void print_id(int id) {

	// keep doing this
	while (true) {

		// attempt to capture the lock
		std::unique_lock<std::mutex> lck(mutThread[id]);
		
		// wait until we have the lock
		while (!isMainReady) cvsThread[id].wait(lck);

		// do work
		std::cout << "thread " << id << '\n';
	}
}

void go() {

}

void StartThreads()
{
	for (int id = 0; id < NumberOfThreads; id++)
	{
		std::unique_lock<std::mutex> lck(mutThread[j]);
		isMainReady = 1;
		cv.notify_all();
	}
}

void WaitThreads()
{
	for (int j = 0; j < NumberOfThreads; j++)
	{

	}
}

int main()
{
	std::thread threads[10];
	// spawn 10 threads:
	for (int i = 0; i < 10; ++i)
		threads[i] = std::thread(print_id, i);

	std::cout << "10 threads ready to race...\n";
	go();                       // go!

	for (auto& th : threads) th.join();

	std::cout << "Threads are done.";

	return 0;
}