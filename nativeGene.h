void particleMaker()
{
	for (int i = 0; i < RBlistSize; i++)
	{
		unordered_map<vector3d, int> um;
		for (int j = 0; j < RBlist[i].specListSize; j++)
		{
			um[RBlist[i].specList[j].initial] = j + 1;
		}
		for (int j = 0; j < RBlist[i].massListSize; j++)
		{
			int a = um[RBlist[i].massList[j].initial];
			if (a > 0)
			{
				particle p(&RBlist[i].massList[j], &RBlist[i].specList[a - 1]);
				RBlist[i].particleList.push_back(p);
			}
		}
	}
}

vector<gene> geneList;

void mutate(float rate)
{
	int a = rate*1000, k = 3;
	double unitDisp = 0.02;
	for (int i = 0; i < geneList.size(); i++)
	{
		int c = rand()%1000;
		if (a < c)
		{
			vector3d v(rand()%k*unitDisp, rand()%k*unitDisp, rand()%k*unitDisp);
			geneList[i].disp += v;
			geneList[i].disp.printtocout();
		}
	}
}

void makeGenes()
{
	for (int i = 0; i < intListSize; i++)
	{
		if (interactionList[i]->type && !interactionList[i]->obs)
		{
			gene g;
			g.type = 0;
			g.vecInt.push_back(i);
			geneList.push_back(g);
		}
	}
	for (int i = 0; i < RBlistSize; i++)
	{
		for (int j = 0; j < RBlist[i].particleList.size(); j++)
		{
			gene g;
			g.type = 1;
			pair<int, int> p(i, j);
			g.vecParticle.push_back(p);
			geneList.push_back(g);
		}
	}
}

void mutate(vector<gene> &gList, float rate)
{
	int a = rate*1000, k = 3;
	double unitDisp = 0.02;
	gList.clear(); gList.resize(geneList.size());
	for (int i = 0; i < gList.size(); i++)
	{
		int c = rand()%1000;
		if (a < c)
		{
			vector3d v(rand()%k*unitDisp, rand()%k*unitDisp, rand()%k*unitDisp);
			gList[i].disp += v;
			gList[i].disp.printtocout();
		}
	}
}

void printGenes(ofstream &ofs)
{
	ofs << geneList.size() << '\n';
	for (int i = 0; i < geneList.size(); i++)
	{
		ofs << geneList[i].type << ' ' <<  geneList[i].shift << '\n';
		geneList[i].disp.print(ofs);
		ofs << geneList[i].vecParticle.size() << '\n';
		for (int j = 0; j < geneList[i].vecParticle.size(); j++)
		{
			ofs << geneList[i].vecParticle[j].first << ' ' << geneList[i].vecParticle[j].second << '\n';
		}
		ofs << geneList[i].vecInt.size() << '\n';
		for (int j = 0; j < geneList[i].vecInt.size(); j++)
		{
			ofs << geneList[i].vecInt[j] << '\n';
		}
	}
}

void phenotype()
{
	for (int i = 0; i < geneList.size(); i++)
	{
		if (geneList[i].type == 0)
		{
			gene &g = geneList[i];
			for (int i = 0; i < g.vecInt.size(); i++)
			{
				interactionList[g.vecInt[i]]->d *= g.shift;
			}
		}
		else
		{
			gene &g = geneList[i];
			for (int i = 0; i < g.vecParticle.size(); i++)
			{
				int a = g.vecParticle[i].first, b = g.vecParticle[i].second;
				rigidBody &rb = RBlist[a]; particle &p = rb.particleList[b];
				p.mass->initial += g.disp;
				p.spec->initial += g.disp;
				if (p.charge != NULL) p.charge->initial += g.disp;
			}
		}
	}
	for (int i = 0; i < intListSize; i++)
	{
		if (!interactionList[i]->type)
		{
			constraint* c = (constraint*)interactionList[i];
			c->locus = *(c->vx[0]);
		}
	}
	for (int i = 0; i < RBlistSize; i++)
	{
		RBlist[i].recenter();
	}
}

void parseGenes(ifstream &ifs)
{
	int a;
	ifs >> a; geneList.clear(); geneList.resize(a);
	for (int i = 0; i < geneList.size(); i++)
	{
		ifs >> geneList[i].type >>  geneList[i].shift;
		geneList[i].disp.read(ifs);
		ifs >> a; geneList[i].vecParticle.resize(a);
		for (int j = 0; j < geneList[i].vecParticle.size(); j++)
		{
			ifs >> geneList[i].vecParticle[j].first >> geneList[i].vecParticle[j].second;
		}
		ifs >> a; geneList[i].vecInt.resize(a);
		for (int j = 0; j < geneList[i].vecInt.size(); j++)
		{
			ifs >> geneList[i].vecInt[j];
		}
	}
}
