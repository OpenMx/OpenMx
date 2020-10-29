#ifndef COVENTRYWISEPARALLEL_H_
#define COVENTRYWISEPARALLEL_H_

#include <queue>
#include "omxState.h"
#include "ConcurrentQueue.h"

template <class T> struct pairMaxCompare {
	bool operator()(const T &x, const T &y){
		return (x.first + x.second) > (y.first + y.second);
	}
};

template <typename CalcEntry>
void CovEntrywiseParallel(int numThreads, CalcEntry &ce)
{
	const bool debug = false;
	std::mutex delayedMutex;
	std::priority_queue< std::pair<int,int>,
		std::vector< std::pair<int,int> >,
		pairMaxCompare< std::pair<int,int> > > delayed;
	ConcurrentDeque< std::pair<int,int> > todo;
	int numCols = ce.getNumCols();
	int numColsStar = triangleLoc1(numCols);
	if (debug) mxLog("CovEntrywiseParallel %d", numThreads);
	// Use type 'long' to ensure that writes don't interfere with each other
	Eigen::Array<long, Eigen::Dynamic, 1> thrDone(numThreads);
	thrDone.setZero();
	// isDone might behave badly so we copy it and track status ourselves
	Eigen::Array<long, Eigen::Dynamic, 1> diagDone(numCols);
	for (int cx=0; cx < numCols; ++cx) {
		diagDone[cx] = ce.isDone(cx,cx);
		todo.push_nolock(std::make_pair(cx,cx));
	}
#pragma omp parallel num_threads(numThreads)
	while (1) {
		int tid = omp_get_thread_num();
		auto t1 = todo.pop();
		if (debug) mxLog("todo.pop -> %d,%d", t1.first, t1.second);
		if (t1.first == -1 && t1.second == -1) {
			int progress = 0;
			while (delayed.size()) {
				std::pair<int,int> top;
				{
					std::unique_lock<std::mutex> mlock(delayedMutex);
					if (delayed.size() == 0) break;
					top = delayed.top();
					if (debug) mxLog("delayed.top is (%d,%d) pending %d",
							 top.first, top.second, int(delayed.size()));
					if (diagDone[top.first] && diagDone[top.second]) {
						delayed.pop();
						if (debug) mxLog("todo.push(%d,%d) %d left",
								 top.first, top.second, int(delayed.size()));
					} else break;
				}
				todo.push_back(top);
				progress += 1;
			}
			if (debug) mxLog("sentinal queued %d (%d/%d)", progress,
			      int(progress + thrDone.sum()), numColsStar);
			if (int(progress + thrDone.sum()) < numColsStar) {
				todo.push_back(std::make_pair(-1, -1));
			} else {
				for (int tx=0; tx < numThreads; ++tx) {
					todo.push_back(std::make_pair(-1, tx));
				}
			}
			continue;
		}
		if (t1.first < 0) break;
		if (t1.first == t1.second) {
			if (!diagDone[t1.first]) {
				try {
					ce.onDiag(t1.first);
				} catch (const std::exception& e) {
					omxRaiseErrorf("%s", e.what());
					if (debug) mxLog("CATCH: %s", e.what());
				} catch (...) {
					omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
					if (debug) mxLog("CATCH: %s line %d: unknown exception", __FILE__, __LINE__);
				}
				diagDone[t1.first] = 1; // regardless of ce.isDone
				if (debug) mxLog("diag %d done", t1.first);
			}
			for (int rx=0; rx < t1.first; ++rx) {
				if (ce.isDone(rx, t1.first)) {
					thrDone[tid] += 1;
					continue;
				}
				if (diagDone(rx)) {
					todo.push_back(std::make_pair(rx, t1.first));
				} else {
					std::unique_lock<std::mutex> mlock(delayedMutex);
					delayed.push(std::make_pair(rx, t1.first));
					if (debug) mxLog("delay %d %d", rx, t1.second);
				}
			}
			if (t1.first == numCols-1) {
				todo.push_back(std::make_pair(-1, -1));
			}
		} else {
			try {
				if (!isErrorRaised()) ce.offDiag(t1.first, t1.second);
			} catch (const std::exception& e) {
				omxRaiseErrorf("%s", e.what());
				if (debug) mxLog("CATCH: %s", e.what());
			} catch (...) {
				omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
				if (debug) mxLog("CATCH: %s line %d: unknown exception", __FILE__, __LINE__);
			}
		}
		thrDone[tid] += 1;
		if (tid == 0) {
			try {
				ce.reportProgress(thrDone.sum());
			} catch (const std::exception& e) {
				omxRaiseErrorf("%s", e.what());
				if (debug) mxLog("CATCH: %s", e.what());
			} catch (...) {
				omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
				if (debug) mxLog("CATCH: %s line %d: unknown exception", __FILE__, __LINE__);
			}
			bool gotInt = Global->interrupted();
			if (gotInt && debug) mxLog("interrupt");
		}
		if (isErrorRaised()) {
			// push_front to stop ASAP
			for (int tx=0; tx < numThreads; ++tx) {
				todo.push_front(std::make_pair(-1, tx));
			}
		}
	}
}

#endif
