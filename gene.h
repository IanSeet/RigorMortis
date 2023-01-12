struct parent
{
	vector<gene> geneList;
	double fitness;

	bool operator < (const parent &otr) const 
	{
	        return fitness < otr.fitness;
	}
};

void crossover(parent &s1, parent &s2, vector<vector<gene> > &offspringList, int target)
{
	offspringList[target].resize(s1.geneList.size());
	for (int j = 0; j < s1.geneList.size(); j++)
	{
		char c = rand()%2;
		if (c) offspringList[target][j] = s1.geneList[j];
		else offspringList[target][j] = s2.geneList[j];
	}
}

struct lossParam
{
	double start, end, depth, equil;
	int target;
};

struct genParam
{
	vector<simulation*> vecsim;
	vector<lossParam*> vecloss;
	int steps, runs, roteqsteps, eqsteps;
	genParam(int _steps, int _runs) : steps(_steps), runs(_runs){eqsteps = 50000, roteqsteps = 0;}
	genParam(int _steps, int _runs, int _eqsteps, int _roteqsteps):
	steps(_steps), runs(_runs), eqsteps(_eqsteps), roteqsteps(_roteqsteps){}
};

static void *genF(void *p)
{
	genParam *sp = (genParam*)p;
	int steps = sp->steps, runs = sp->runs;
	const int minsteps = 20000, eqsteps = sp->eqsteps, areq = sp->eqsteps, roteqsteps = sp->roteqsteps, write = steps/10;
	bool roteq = 0; if (sp->roteqsteps > 0) roteq = 1;
	string s, fout, ftype = "top";
	for (int j = 0; j < sp->vecsim.size(); j++)
	{
		simulation &sim = *(sp->vecsim[j]);
		cout << "starting simulation " << sim.name << '\n';
		fout = sim.name;
		//cout << sim.SP() << '\n';
		//sim.printContributions(0);
		sim.SDmin(minsteps);
		cout << "min complete\n";
		sim.equilibriate(eqsteps, eqsteps/10);
		if (roteq)
		{
			sim.roteq(roteqsteps, roteqsteps/10);
		}
		sim.optConfig.init(sim.RBlist, sim.RBlistSize);
		sim.optConfig.assign();
		sim.combinedEnergies.resize(steps/write + 1);
		//parseInputCommon(ftype, fout, sim);
		sim.maxRuns = runs;
		s = fout + "Lossf";
		ofstream ofsl(s.c_str());
		for (int i = 0; i < runs; i++)
		{
			sim.runNo = i + 1;
			sim.optConfig.overwrite(); sim.catpInit();
			sim.areq(eqsteps, eqsteps/10);
			sim.langevinIntegrator(steps, write, 0);
			for (int j = 0; j < sim.MDtrajectory.size(); j++) sim.combinedEnergies[j].add(sim.MDtrajectory[j]);
			sim.expenditureList.push_back(sim.expenditure);
			//ofsl << "runNo: " << i << '\n'; 
			sim.findObservables();
		}
		cout << "lossF: " << sim.lossFunction(0.25, 0.75, -30, 90.0, 0.1, 0, sim.expenditureList.size() - 1, ofsl) << '\n';
		s = fout + ".xyz";
		ofstream ofs(s.c_str());
		for (int i = 0; i < sim.MDtrajectory.size(); i++)
		{
			sim.printConfigxyz(i, sim.MDtrajectory, ofs); 
		}
		s = fout + "Energy";
		ofstream ofs2(s.c_str());
		for (int i = 0; i < sim.combinedEnergies.size(); i++)
		{
			sim.combinedEnergies[i].divide(runs);
			ofs2 << i << ' ' << sim.combinedEnergies[i].TKE << ' ' << sim.combinedEnergies[i].RKE << ' ' << sim.combinedEnergies[i].V << ' ' <<
			sim.combinedEnergies[i].TKE + sim.combinedEnergies[i].RKE + sim.combinedEnergies[i].V << ' ' << sim.combinedEnergies[i].expenditure << '\n';
		}
		s = fout + "Observables";
		ofstream ofs3(s.c_str());
		sim.printObservables(ofs3, writePOL);
		s = fout + "Expenditures";
		ofstream ofs4(s.c_str());
		sim.printExpenditures(ofs4);
	}
}

void genAlgo(int generations, int popSize, float mutRate, simulation &seed, genAlgoParam &gap)
{
	vector<vector<parent> > parentList(generations);
	vector<vector<vector<gene> > > offspringList(generations);
	for (int i = 0; i < generations; i++)
	{
		offspringList[i].resize(popSize);
	}
	vector<vector<simulation> > vs(generations);
	for (int i = 0; i < generations; i++) vs[i].resize(popSize);
	
	const int maxthreads = gap.maxthreads;
	vector<pthread_t> genThreads(maxthreads);
	vector<genParam> vp;
	const double threshold = 0.25;
	int limit = threshold * popSize; if (limit == 0) limit = 1;
	for (int j = 0; j < generations; j++)
	{
		for (int k = 0; k < maxthreads; k++)
		{
			genParam gp(gap.steps, gap.runs, gap.eqsteps, gap.roteqsteps);
			for (int i = k; i < popSize; i += maxthreads)
			{
				seed.deepCopy(vs[j][i]);
				if (j > 0)
				{
					vs[j][i].geneList = offspringList[j - 1][i];
				}
				vs[j][i].initDecompAll();
				vs[j][i].particleMaker(); vs[j][i].makeGenes();
				vs[j][i].name = seed.name + "_c" + to_string(i);
				//cout << vs[i].name << '\n';
				//vs[i].mutate(gp.mutrate);
				//vs[i].phenotype();
				cout << "mutation complete\n";
				gp.vecsim.push_back(&vs[j][i]);
			}
			vp.push_back(gp);
		}
		cout << "initalising threads\n";
		for (int i = 0; i < maxthreads; i++)
		{
			int rc = pthread_create(&genThreads[i], NULL, genF, &vp[i]);
			//cout << &vs[i] << '\n';
	  		if (rc)
			{
				cout << "Error: unable to create thread " << rc << '\n';
				exit(-1);
	  		}
		}
		for (int i = 0; i < maxthreads; i++) 
		{
			pthread_join(genThreads[i], NULL);
		}
		cout << "threads joined\n";
		
		vector<parent> temp(popSize);
		for (int i = 0; i < popSize; i++)
		{
			temp[i].geneList = vs[j][i].geneList;
			temp[i].fitness = vs[j][i].fitness;
		}
		sort(temp.begin(), temp.end());
		cout << "parents sorted\n";
		parentList[j].resize(limit);
		
		for (int i = 0; i < limit; i++) parentList[j][i] = temp[i];
		vector<vector<gene> > &g = offspringList[j];
		for (int i = 0; i < popSize; i++)
		{
			int parent1 = rand()%limit, parent2 = rand()%limit;
			crossover(parentList[j][parent1], parentList[j][parent2], g, i);
		}
		cout << "crossover done\n";
	}
}
