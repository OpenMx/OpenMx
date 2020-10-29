/*
 * Copyright 2016-2019 by the individuals mentioned in the source code history
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
		Connectedness(std::vector<int> &_region, SubgraphType &_connected,
			      int size, bool _verbose) :
			region(_region), connected(_connected), verbose(_verbose)
		{
			region.assign(size, -1);
			connected.clear();
			subgraphs = size;
		}

		void log()
		{
			if (!verbose) return;
			mxLog("subgraph count = %d", subgraphs);
			Eigen::Map< Eigen::VectorXi > regionMap(region.data(), region.size());
			mxPrintMat("region", regionMap);
			for (int cx=0; cx < (int) connected.size(); ++cx) {
				if (!connected[cx].size()) continue;
				std::set<int> &as1set = connected[cx];
				std::string str = string_snprintf("group %d:", cx);
				for (std::set<int>::iterator it2 = as1set.begin(); it2 != as1set.end(); ++it2) {
					str += string_snprintf(" %d", *it2);
				}
				str += "\n";
				mxLogBig(str);
			}
		}

		int getSizeIfConnected(int ax, int bx)
		{
			if (region[ax] == -1 && region[bx] == -1) return 2;
			if (region[ax] == region[bx]) return connected[ region[ax] ].size();
			if (region[ax] == -1) return connected[ region[bx] ].size() + 1;
			if (region[bx] == -1) return connected[ region[ax] ].size() + 1;
			return (connected[ region[ax] ].size() +
				connected[ region[bx] ].size());
		}

		void connect(int ax, int bx)
		{
			if (ax == bx) mxThrow("Cannot connect %d to itself", ax);
			if (region[ax] == -1) preconnect(ax);
			if (region[ax] == region[bx]) return; //already connected
			subgraphs -= 1;
			//if (subgraphs < 1) mxThrow("problem");
			if (region[bx] == -1) {
				region[bx] = region[ax];
				connected[ region[ax] ].insert(bx);
				if (verbose) {
					mxLog("add %d to region %d", bx, region[ax]);
				}
			} else {
				int rax = region[ax];
				int rbx = region[bx];
				if (rbx > rax) std::swap(rbx, rax);
				// as1 > as2
				if (rbx != rax) {
					// merge to as2
					if (verbose) {
						mxLog("merge region %d (%d elem) to region %d (%d elem)",
						      rax, (int)connected[rax].size(),
						      rbx, (int)connected[rbx].size());
					}
					std::set<int> &as1set = connected[rax];
					std::set<int> &as2set = connected[rbx];
					for (std::set<int>::iterator it2 = as1set.begin(); it2 != as1set.end(); ++it2) {
						region[*it2] = rbx;
						as2set.insert(*it2);
					}
					as1set.clear();
				}
			}
			if (verbose) log();
		}

		int numSubgraphs() const { return subgraphs; }
	};
}

#endif
