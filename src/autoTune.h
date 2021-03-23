#ifndef u_AutoTune_H_
#define u_AutoTune_H_

#include "omxDefines.h"

template <typename T>
class AutoTune {
  const char *name;
	const int ELAPSED_HISTORY_SIZE;
  bool used;
  nanotime_t startTime;
	std::vector<nanotime_t> elapsed0;
	std::vector<nanotime_t> elapsed1;
	int curElapsed;
	int numThreadsBookmark;
	int maxAvailThreads;
	int verbose;
  int curNumThreads;
  std::unique_ptr<T> workPtr;

  void start()
  {
    used = true;
		startTime = get_nanotime();
		int th = std::max(1, numThreadsBookmark - curElapsed % 2);
    if (th != curNumThreads) {
      workPtr->setNumThreads(th);
      curNumThreads = th;
    }
  }

  void finish()
  {
		double el1 = (get_nanotime() - startTime);
    if (curElapsed < ELAPSED_HISTORY_SIZE * 2) {
      if (verbose >= 2) {
        mxLog("%s: test[%d] curNumThreads=%d %fms",
              name, curElapsed, curNumThreads, el1/1000000.0);
      }
			int repl = curElapsed / 2;
			if (curElapsed % 2) {
				elapsed1[repl] = el1;
			} else {
				elapsed0[repl] = el1;
			}
			curElapsed += 1;
			if (curElapsed == ELAPSED_HISTORY_SIZE * 2) {
				std::sort(elapsed0.begin(), elapsed0.end());
				std::sort(elapsed1.begin(), elapsed1.end());
				double e0 = elapsed0[elapsed0.size()/2];
				double e1 = elapsed1[elapsed1.size()/2];
        if (verbose) {
          mxLog("%s: took %fms with %d threads and %fms with %d threads",
                name, e0/1000000.0, numThreadsBookmark, e1/1000000.0, std::max(1, numThreadsBookmark-1));
        }
				if (e0 > e1 && numThreadsBookmark > 1) {
					numThreadsBookmark = std::max(numThreadsBookmark - 1, 1);
					if (numThreadsBookmark > 1) curElapsed = 0;
				}
        if (verbose && curElapsed > 0) {
          mxLog("%s: looks like %d threads offer the best performance",
                name, numThreadsBookmark);
        }
			}
		}
  }

 public:
  AutoTune(const char *u_name) :
  name(u_name), ELAPSED_HISTORY_SIZE(3), used(false) {
    setMaxThreads(1);
  }

  ~AutoTune() {
    if (used) {
      diagParallel(OMX_DEBUG, "%s: used %d/%d threads",
                   name, numThreadsBookmark, maxAvailThreads);
    } else {
      diagParallel(OMX_DEBUG, "%s: not used", name);
    }
  }

  T &work() { return *workPtr; }

  void setWork(std::unique_ptr<T> w) {
    curNumThreads = -1;
    workPtr = std::move(w);
  }

  void setMaxThreads(int th)
  {
    if (used) mxThrow("%s: already used", name);
    th = std::max(th, 1);
    maxAvailThreads = th;
		verbose = th>1 && Global->parallelDiag;

    if (workPtr) {
      numThreadsBookmark = std::min(th, workPtr->getMaxUsableThreads());
    } else {
      numThreadsBookmark = 1;
    }
    if (numThreadsBookmark < 1) numThreadsBookmark = 1;
		if (numThreadsBookmark == 1) {
			curElapsed = ELAPSED_HISTORY_SIZE * 2;
      elapsed0.clear();
      elapsed1.clear();
		} else {
      curElapsed = 0;
			elapsed0.resize(ELAPSED_HISTORY_SIZE);
			elapsed1.resize(ELAPSED_HISTORY_SIZE);
		}
  }

  template<typename ...Args> void operator()(Args&&... args)
    {
      start();
      (*workPtr)(std::forward<Args>(args)...);
      finish();
    }
};

#endif
