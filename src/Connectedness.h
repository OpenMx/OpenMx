/*
 * Copyright 2016 The OpenMx Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _CONNECTEDNESS_HPP_
#define _CONNECTEDNESS_HPP_

#include <set>

namespace UndirectedGraph {

	// Finds connected subgraphs of a unidirected graph.  Nodes
	// are labeled 0 to (# of nodes)-1.  The region vector holds
	// the assignment of nodes to connected regions. Region -1 is
	// a special region that indicates that that node is not
	// connected to any other nodes. The nodes that comprise a
	// region are given in connected[ region[ax] ]. The purpose of
	// region -1 is to avoid allocation of a std::set for every
	// node. It is assumed that there are many fewer subgraphs
	// than there are nodes.

	class Connectedness {
	public:
		typedef std::vector< std::set<int> > SubgraphType;
	private:
		std::vector<int> &region;
		SubgraphType &connected;
		bool verbose;
		int subgraphs;

		void preconnect(int ax)
		{
			if (region[ax] != -1) return;
			region[ax] = connected.size();
			connected.resize(connected.size() + 1);
			connected[ region[ax] ].insert(ax);
			if (verbose) {
				mxLog("preconnect %d to region %d", ax, region[ax]);
			}
		}
	public:
		Connectedness(std::vector<int> &region, SubgraphType &connected,
			      int size, bool verbose) :
			region(region), connected(connected), verbose(verbose)
		{
			region.assign(size, -1);
			connected.clear();
			subgraphs = size;
		}

		void log()
		{
			if (!verbose) return;
			Eigen::Map< Eigen::VectorXi > regionMap(region.data(), region.size());
			mxPrintMat("region", regionMap);
		}

		void connect(int ax, int bx)
		{
			if (ax == bx) Rf_error("Cannot connect %d to itself", ax);
			if (region[ax] == -1) preconnect(ax);
			if (region[ax] == region[bx]) return; //already connected
			subgraphs -= 1;
			if (region[bx] == -1) {
				region[bx] = region[ax];
				connected[ region[ax] ].insert(bx);
				if (verbose) {
					mxLog("add %d to region %d", bx, region[ax]);
				}
			} else {
				if (region[bx] > region[ax]) std::swap(region[bx], region[ax]);
				// as1 > as2
				if (region[bx] != region[ax]) {
					if (verbose) {
						mxLog("merge region %d (%d elem) to region %d (%d elem)",
						      region[ax], (int)connected[region[ax]].size(),
						      region[bx], (int)connected[region[bx]].size());
					}
					// merge to as2
					std::set<int> &as1set = connected[region[ax]];
					std::set<int> &as2set = connected[region[bx]];
					for (std::set<int>::iterator it2 = as1set.begin(); it2 != as1set.end(); ++it2) {
						region[*it2] = region[bx];
						as2set.insert(*it2);
					}
					as1set.clear();
				}
			}
		}

		int numSubgraphs() const { return subgraphs; }
	};
};

#endif
