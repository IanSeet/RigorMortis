vector<pair<int, int> > bondSpecList;
map<int, int> specMap;
vector<int> colour;
void bondGen()
{
	for (int i = 0; i < bondList.size(); i++)
	{
		pair<int, int> p1, p2;
		p1.first = bondList[i].rbx[0] - RBlist; p1.second = bondList[i].px[0].second;
		p2.first = bondList[i].rbx[1] - RBlist; p2.second = bondList[i].px[1].second;
		int hashf1 = ((p1.first + 1) << 16) + p1.second, hashf2 = ((p2.first + 1) << 16) + p2.second;
		if (specMap.count(hashf1) == 0)
		{
			bondSpecList.push_back(p1);
			adjList.resize(adjList.size() + 1);
			specMap[hashf1] = adjList.size() - 1;
		}
		if (specMap.count(hashf2) == 0)
		{
			bondSpecList.push_back(p2);
			adjList.resize(adjList.size() + 1);
			specMap[hashf2] = adjList.size() - 1;
		}
		int small = min(specMap[hashf1], specMap[hashf2]), big = max(specMap[hashf1], specMap[hashf2]);
		//cout << i << ' ' << bondList.size() << ' ' << small << ' ' << p1.first << ' ' << p1.second << ' ' << big << ' ' << p2.first << ' ' << p2.second << '\n';
		pair <int, int> p(big, 1);
		if (bondList[i].isChain) p.second = 4;
		adjList[small].push_back(p);
	}
}

void dfsPriority(vector<vector<int> > &adjList, vector<int> &priority, int x, int &curr)
{
	if (priority[x] != -1) return;
	priority[x] = curr;
	for (int i = 0; i < adjList[x].size(); i++)
	{
		if (priority[adjList[x][i]] == -1)
		{
			curr++;
			dfsPriority(adjList, priority, adjList[x][i], curr);
		}
	}
}

bool colourAux(int n, vector<int> &colour, int x, vector<vector<int> > &adjList, vector<int> &priority)
{
	//cout << "colourAux: " << x << '\n';
	for (int c = x%n; c < n + x%n; c++)
	{
		c %= n;
		bool clash = 0, succ = 1;
		for (int i = 0; i < adjList[x].size(); i++)
		{
			if (priority[adjList[x][i]] < priority[x] && colour[adjList[x][i]] == c)
			{
				clash = 1;
				break;
			}
		}
		if (!clash)
		{
			colour[x] = c;
			for (int i = 0; i < adjList[x].size(); i++)
			{
				if (priority[adjList[x][i]] > priority[x])
				{
					succ = colourAux(n, colour, adjList[x][i], adjList, priority);
					if (!succ) break;
				}
			}
			if (succ) return true;
		}
	}
	colour[x] = -1;
	return false;
}

bool colouring(int n, vector<vector<int> > &adjList)
{
	/*for (int i = 0; i < adjList.size(); i++)
	{
		for (int j = 0; j < adjList[i].size(); j++) cout << adjList[i][j] << ' ';
		cout << '\n';
	}*/
	vector<int> colourNat, priority; colourNat.assign(adjList.size(), -1); priority.assign(adjList.size(), -1);
	int curr = 0;
	for (int i = 0; i < adjList.size(); i++) dfsPriority(adjList, priority, i, curr);
	//for (int i = 0; i < priority.size(); i++) cout << priority[i] << '\n';
	bool b = 1;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (colourNat[i] == -1) b &= colourAux(n, colourNat, i, adjList, priority);
	}
	if (b)
	{
		cout << "Colouring found\n";
		colour.resize(colourNat.size());
		for (int i = 0; i < colourNat.size(); i++)
		{
			colour[i] = colourNat[i];
		}
		return true;
	}
	else
	{
		cout << "No colouring found for " << n << '\n';
		return false;
	}
}

void topdb(ofstream &ofs)
{
	for (int i = 0; i < RBlistSize; i++) for (int j = 0; j < RBlist[i].massListSize; j++)
	{
		ofs << "HETATM" << setw(4) << i + j + 1 << setw(3) << "C" << setw(6) << "R" << i << setw(6) << 1 << setw(12) << fixed << setprecision(3) << RBlist[i].massList[j].decomp.i << setw(8) << RBlist[i].massList[j].decomp.j << setw(8) << RBlist[i].massList[j].decomp.k << setw(6) << "1.00" << setw(6) << "0.00" << setw(12) << "C\n"; 
	}
}

void togjf(ofstream &ofs)
{
	ofs << "# opt pm6\n\ntitle\n\n0 1\n";
	for (int i = 0; i < RBlistSize; i++) 
	{
		for (int j = 0; j < RBlist[i].massListSize; j++) ofs << " " << 'H' << right << setw(24) << fixed << setprecision(8) << RBlist[i].massList[j].decomp.i << setw(16) << RBlist[i].massList[j].decomp.j << setw(16) << RBlist[i].massList[j].decomp.k << '\n';
		for (int j = 1; j < RBlist[i].specListSize; j++) ofs << " " << 'H' << right << setw(24) << fixed << setprecision(8) << RBlist[i].specList[j].decomp.i << setw(16) << RBlist[i].specList[j].decomp.j << setw(16) << RBlist[i].specList[j].decomp.k << '\n';
	}
}

void printConfigGjf(int n, vector<config> &trajectory, ofstream &ofs)
{
	config& cf = trajectory[n];
	ofs << "# opt pm6\n\ntitle\n\n0 1\n";
	for (int i = 0; i < RBlistSize; i++)
	{
		RBlist[i] = cf.vs[i]; RBlist[i].decompose();
		for (int j = 0; j < RBlist[i].massListSize; j++) ofs << " " << 'H' << right << setw(24) << fixed << setprecision(8) << RBlist[i].massList[j].decomp.i << setw(16) << RBlist[i].massList[j].decomp.j << setw(16) << RBlist[i].massList[j].decomp.k << '\n';
		for (int j = 1; j < RBlist[i].specListSize; j++) ofs << " " << 'H' << right << setw(24) << fixed << setprecision(8) << RBlist[i].specList[j].decomp.i << setw(16) << RBlist[i].specList[j].decomp.j << setw(16) << RBlist[i].specList[j].decomp.k << '\n';
	}
}

void printConfigxyz(int n, vector<config> &trajectory, ofstream &ofs)
{
	bool changeType = 1;
	config& cf = trajectory[n];
	int i, j, count = 0;
	for (i = 0; i < RBlistSize; i++)
	{
		for (j = 0; j < RBlist[i].massListSize; j++) count++;
		//for (j = 1; j < RBlist[i].special.size(); j++) count++;
	}
	count += bondSpecList.size();
	ofs << count << "\n\n";
	for (i = 0; i < RBlistSize; i++)
	{
		RBlist[i] = cf.vs[i]; RBlist[i].decompose();
		for (int j = 0; j < RBlist[i].massListSize; j++)
		{
			string displayAs;
			if (!changeType) displayAs = "H";
			else displayAs = paramTable[RBlist[i].massList[j].type].displayAs;
			ofs << " " << displayAs << right << setw(24) << fixed << setprecision(8) << RBlist[i].massList[j].decomp.i << setw(16) << RBlist[i].massList[j].decomp.j << setw(16) << RBlist[i].massList[j].decomp.k << '\n';
		}
		//for (int j = 1; j < RBlist[i].specListSize; j++) ofs << " " << colour[i] << right << setw(24) << fixed << setprecision(8) << RBlist[i].specList[j].decomp.i << setw(16) << RBlist[i].specList[j].decomp.j << setw(16) << RBlist[i].specList[j].decomp.k << '\n';
	}
	for (int u = 0; u < bondSpecList.size(); u++)
	{
		i = bondSpecList[u].first; j = bondSpecList[u].second;
		ofs << " " << 'H' << right << setw(24) << fixed << setprecision(8) << RBlist[i].specList[j].decomp.i << setw(16) << RBlist[i].specList[j].decomp.j << setw(16) << RBlist[i].specList[j].decomp.k << '\n';
	}
}

void printTopolMol2(ofstream &ofs, char idx, bool changeType)
{
	string res[25] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLY", "GLU", "GLN", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ADE",
	"CYT", "GUA", "THY", "URA"};
	int i, j, k, l = 0, m = 0, n = 0, count = 0, bonds = 0;
	for (i = 0; i < RBlistSize; i++)
	{
		for (j = 0; j < RBlist[i].massListSize; j++) count++;
		//for (j = 1; j < RBlist[i].special.size(); j++) count++;
		if (idx < 2) for (j = 0; j < RBlist[i].adjList.size(); j++) bonds += RBlist[i].adjList[j].size();
	}
	if (idx == 0)
	{
		for (i = 0; i < adjList.size(); i++) bonds += adjList[i].size();
	}
	else if (idx == 2)
	{
		for (i = 0; i < adjList.size(); i++) for (j = 0; j < adjList[i].size(); j++) if (adjList[i][j].second == 1) bonds++;
	}
	else if (idx == 3)
	{
		for (i = 0; i < adjList.size(); i++) for (j = 0; j < adjList[i].size(); j++) if (adjList[i][j].second == 4) bonds++;
	}
	count += bondSpecList.size();
	ofs << "#\n#\n#\n\n#\n#\n\n@<TRIPOS>MOLECULE\nPanda\n";
	ofs << count << ' ' << bonds << '\n';
	ofs << "SMALL\nNO CHARGES\n\n\n@<TRIPOS>ATOM\n";
	for (i = 0; i < RBlistSize; i++) for (j = 0; j < RBlist[i].massListSize; j++)
	{
		string displayAs;
		if (!changeType) displayAs = "H";
		else displayAs = paramTable[RBlist[i].massList[j].type].displayAs;
		n++;
		string res2 = res[colour[i]%25]; if (idx == 2) res2 = "GLY"; if (idx == 3) res2 = "GLU";
		ofs << n << ' ' << displayAs << n << ' ' << RBlist[i].massList[j].decomp.i << ' ' << RBlist[i].massList[j].decomp.j << ' ' << RBlist[i].massList[j].decomp.k << ' ' << displayAs << ' ' << i  << ' ' << res2 << '\n'; 
	}
	k = n;
	for (i = 0; i < bondSpecList.size(); i++)
	{
		k++;
		int f = bondSpecList[i].first, s = bondSpecList[i].second;
		string res2 = res[colour[f]%25]; if (idx == 2) res2 = "GLY"; if (idx == 3) res2 = "GLU";
		ofs << k << ' ' << 'H' << k << ' ' << RBlist[f].specList[s].decomp.i << ' ' << RBlist[f].specList[s].decomp.j << ' ' << RBlist[f].specList[s].decomp.k << ' ' << 'H' << ' ' << f  << ' ' << res2 << '\n';
	}
	ofs << "@<TRIPOS>BOND";
	if (idx < 2)
	{
		for (i = 0; i < RBlistSize; i++) 
		{
			for (j = 0; j < RBlist[i].adjList.size(); j++) for (k = 0; k < RBlist[i].adjList[j].size(); k++)
			{
				l++;
				ofs << '\n' << l << ' ' << m + j + 1 << ' ' << m + RBlist[i].adjList[j][k] + 1 << " 2";  
			}
			m += RBlist[i].massListSize;
		}
	}
	if (idx == 0)
	{
		for (i = 0; i < adjList.size(); i++) for (j = 0; j < adjList[i].size(); j++)
		{
			l++;
			ofs << '\n' << l << ' ' << n + i + 1 << ' ' << n + adjList[i][j].first + 1 << " 1";
		}
	}
	else if (idx == 2)
	{
		for (i = 0; i < adjList.size(); i++) for (j = 0; j < adjList[i].size(); j++)
		{
			if (adjList[i][j].second == 1)
			{
				l++;
				ofs << '\n' << l << ' ' << n + i + 1 << ' ' << n + adjList[i][j].first + 1 << " 1";
			}
		}
	}
	else if (idx == 3)
	{
		for (i = 0; i < adjList.size(); i++) for (j = 0; j < adjList[i].size(); j++)
		{
			if (adjList[i][j].second == 4)
			{
				l++;
				ofs << '\n' << l << ' ' << n + i + 1 << ' ' << n + adjList[i][j].first + 1 << " 4";
			}
		}
	}
	ofs << '\n';
}

void printCurrConfig(ofstream &ofs)
{
	for (int i = 0; i < RBlistSize; i++)
	{
		RBlist[i].center.print(ofs);
		RBlist[i].momentum.print(ofs);
		RBlist[i].orient.print(ofs);
		RBlist[i].angmom.print(ofs);
	}
}

void printConfigCout(int n, vector<config> &trajectory, ofstream &ofs)
{
	config& cf = trajectory[n];
	for (int i = 0; i < RBlistSize; i++)
	{
		cf.vs[i].center.printtocout();
		cf.vs[i].momentum.printtocout();
		cf.vs[i].orient.printtocout();
		cf.vs[i].angmom.printtocout();
	}
}

void printObservables(ofstream &ofs, bool writePOL)
{
	if (!writePOL)
	{
		map<string, pair<double, double> >::iterator it;
		for (it = observableMatrix.begin(); it != observableMatrix.end(); ++it)
		{
			for (int j = 0; j < it->first.length(); j++) ofs << (int)((unsigned char)it->first[j]) << ' ';
			ofs << it->second.first << ' ' << it->second.first / it->second.second << '\n';
		}
	}
	else
	{
		//cout << "polsize: " << printObservableList.size() << '\n';
		for (int i = 0; i < printObservableList.size(); i++)
		{
			ofs << "Timestep No: " << i << '\n';
			for (int j = 0; j < printObservableList[i].size(); j++)
			{
				//cout << "pol: " << i << ' ' << j << '\n';
				ofs << "Observable No: " << j << '\n';
				map<float, int>::iterator it;
				for (it = printObservableList[i][j].begin(); it != printObservableList[i][j].end(); ++it)
				{
					ofs << it->first << ": " << it->second << '\n';
				}
				ofs << '\n';
			}
			ofs << '\n';
		}
	}
}

/*inline void printCost(ofstream &ofs)
{
	for (int i = 0; i < printObservableList.size(); i++)
	{
		for (int j = 0; j < printObservableList[i].size(); j++) ofs << (int)printObservableList[i][j] << ' ';
		ofs << '\n';
	}
}*/

void printExpenditures(ofstream &ofs)
{
	double mean = 0, stdev = 0;
	//cout << "exp1\n";
	for (int i = 0; i < expenditureList.size(); i++)
	{
		mean += expenditureList[i];
		ofs << expenditureList[i];
		for (int j = 0; j < finalObservableList[i].size(); j++)
		{
			ofs << ' ';// << j << ":";
			for (int k = 0; k < finalObservableList[i][j].size(); k++) ofs << ' ' << (int)finalObservableList[i][j][k];
			ofs << ' ';
		}
		ofs << '\n';
	}
	mean /= expenditureList.size();
	//cout << "exp2\n";
	for (int i = 0; i < expenditureList.size(); i++)
	{
		stdev += fabs(expenditureList[i] - mean);
	}
	stdev /= expenditureList.size();
	ofs << "Mu: " << mean << " Sigma: "  << stdev << '\n';
}

void printExtDipVectors(ofstream &ofs, bool extend)
{
	for (int i = 0; i < intListSize; i++)
	{
		if (interactionList[i]->type == 5)
		{
			vector3d v1 = interactionList[i]->rbx[0]->specList[0].decomp, v2 = *(interactionList[i]->vx[0]);
			if (extend)
			{
				vector3d v3 = v2 - v1;
				v2 += v3; v1 -= v3;
			}
			v1.printNoNewline(ofs);
			v2.print(ofs);
		}
	}
}

void printExtDipVectorsFrame(ofstream &ofs, bool extend, vector<config> &traj)
{
	for (int frame = 0; frame < traj.size(); frame++)
	{
		ofs << "Frame " << frame << '\n';
		config &cf = traj[frame];
		for (int i = 0; i < RBlistSize; i++)
		{
			RBlist[i] = cf.vs[i]; RBlist[i].decompose();
		}
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->type == 5)
			{
				if (isActive && interactionList[i]->active != 0 && active.count(interactionList[i]->active) == 0) {}
				else
				{
					vector3d v1 = interactionList[i]->rbx[0]->specList[0].decomp, v2 = *(interactionList[i]->vx[0]);
					if (extend)
					{
						vector3d v3 = v2 - v1;
						v2 += v3; v1 -= v3;
					}
					v1.printNoNewline(ofs);
					v2.print(ofs);
				}
			}
		}
		for (int i = 0; i < constrainAxisTimeList.size(); i++)
		{
			vector3d v1 = constrainAxisTimeList[i].rbx[0]->specList[0].decomp, v2 = *(constrainAxisTimeList[i].vx[0]);
			if (extend)
			{
				vector3d v3 = v2 - v1;
				v2 += v3; v1 -= v3;
			}
			v1.printNoNewline(ofs);
			v2.print(ofs);
		}
	}
}
