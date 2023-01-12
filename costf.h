double lossTraj(double start, double end, double depth, double equil, int target, ofstream &ofs)
{
	vector<vector<map<float, int> > > & pol = printObservableList;
	int a = start * pol.size(), b = end * pol.size();
	double loss = 0.0;
	for (int i = 0; i < pol.size() - 1; i++)
	{
		ofs << "Step no: " << i;
		map<float, int>::iterator it;
		if (i < a || i > b)
		{
			ofs << " Relaxed\n";
			for (it = pol[i][target].begin(); it != pol[i][target].end(); ++it)
			{
				float f = (it->first - equil)*(it->first - equil)*it->second;
				loss += f;
				ofs << it->first << ' '<< it->second << ' ' << f << '\n';
			}
		}
		else
		{
			ofs << " Engaged\n";
			for (it = pol[i][target].begin(); it != pol[i][target].end(); ++it)
			{
				if ((depth > 0 && it->first < equil + depth) || (depth < 0 && it->first > equil + depth))
				{
					float f = (it->first - (equil + depth))*(it->first - (equil + depth))*it->second;
					loss += f;
					ofs << it->first << ' ' << it->second << ' ' << f << '\n';
				}
			}
		}
		//cout << i << ' ' << pol[i][target] << ' ' << loss << '\n';
	}
	return loss;
}

double fitness = 0;

double lossFunction(double start, double end, double depth, double equil, double trajFrac, int target, int curr, ofstream & ofs)
{
	//cout << "expcost: " << (1 - trajFrac)*expenditureList[0] << '\n';
	fitness = trajFrac*lossTraj(start, end, depth, equil, target, ofs)/100000 + (1 - trajFrac)*expenditureList[curr];
	return fitness;
}
