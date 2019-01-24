#ifndef COVENTRYWISEPARALLEL_H_
#define COVENTRYWISEPARALLEL_H_

#include "omxState.h"
#include "ConcurrentQueue.h"

template <typename CalcEntry>
void CovEntrywiseParallel(int numThreads, CalcEntry &ce)
{
	const bool debug = false;
	ConcurrentDeque< std::pair<int,int> > todo;
	// Use type 'long' to ensure that writes don't interfere with each other
	Eigen::Array<long, Eigen::Dynamic, 1> thrDone(numThreads);
	thrDone.setZero();
	int numCols = ce.getNumCols();
	int numColsStar = triangleLoc1(numCols);
	for (int cx=0; cx < numCols; ++cx) {
		todo.push_nolock(std::make_pair(cx,cx));
	}
#pragma omp parallel num_threads(numThreads)
	while (1) {
		int tid = omp_get_thread_num();
		auto t1 = todo.pop();
		if (t1.first < 0) break;
		if (debug) mxLog("todo.pop -> %d,%d", t1.first, t1.second);
		if (t1.first == t1.second) {
			if (!ce.isDone(t1.first, t1.first)) {
				try {
					ce.var(t1.first);
				} catch (const std::exception& e) {
					omxRaiseErrorf("%s", e.what());
				} catch (...) {
					omxRaiseErrorf("CovEntrywiseParallel: unknown exception");
				}
			}
			thrDone[tid] += 1;
			for (int rx=0; rx < t1.first; ++rx) {
				if (ce.isDone(rx, t1.first)) { thrDone[tid] += 1; continue; }
				if (debug) mxLog("todo.push_front(%d,%d)", rx, t1.first);
				if (ce.isDone(t1.first, t1.first)) {
					todo.push_front(std::make_pair(rx, t1.first));
				} else {
					todo.push_back(std::make_pair(rx, t1.first));
				}
			}
		} else {
			if (!ce.isDone(t1.first, t1.first)) {
				if (debug) mxLog("requeue %d %d", t1.first, t1.second);
				todo.push_back(t1);
			} else {
				try {
					ce.cov(t1.first, t1.second);
				} catch (const std::exception& e) {
					omxRaiseErrorf("%s", e.what());
				} catch (...) {
					omxRaiseErrorf("CovEntrywiseParallel: unknown exception");
				}
				thrDone[tid] += 1;
			}
		}
		int thrDoneSum = thrDone.sum();
		if (tid == 0) {
			ce.reportProgress(thrDoneSum);
			bool gotInt = omxGlobal::interrupted();
			if (gotInt && debug) mxLog("interrupt");
		}
		if (thrDoneSum == numColsStar || isErrorRaised()) {
			for (int tx=0; tx < numThreads; ++tx) {
				todo.push_front(std::make_pair(-1, tx));
			}
		}
	}
}

#endif
