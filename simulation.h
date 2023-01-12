class simulation
{
	//Simulation units: Energy - kT at 298K, Temperature - kelvin, Distance - 0.5nm, Charge - elementary charge
	public:
	int currStep, RBlistSize, dummyListSize, intListSize, maxTimeStep, chain;
	double contribs[6]; double LJcontrib;
	bool printStep, updateCellList; //If a particle breaches the boundaries of a cell, trigger update
	vector3d minDim, maxDim, boxDim; //Smallest and largest values of x, y and z reached. Needed for grid cell calculations.
	rigidBody * RBlist;
	inline void haltAll()
	{
		for (int i = 0; i < RBlistSize; i++) RBlist[i].halt();
	}
	config optConfig;
	vector<double> expenditureList;
	vector<config> SDtrajectory, MDtrajectory, eqtrajectory; //Stores trajectory
	vector<energyConfig> combinedEnergies; //Averaged energies from multiple runs
	double expenditure; //Energy expenditure
	vector<constrainAxisTime> constrainAxisTimeList;
	vector<bond> bondList;
	
	struct trobs
	{
		vector<int> vi;
		int last;
		trobs()
		{
			last = -1;
		}
	};
	
	vector<interaction*> observables; unordered_map<float, trobs> triggerObs;
	vector< vector<unsigned char> > finalObservables;
	vector<vector<vector<unsigned char> > > finalObservableList;
	vector<vector<map<float, int> > > printObservableList;
	vector<vector<float> > polAverage;
	vector<map<float, int> > poltemp;
	map<string, pair<double, double> > observableMatrix;
	double interactionTime, integratorTime, liTime, cellListTime;
	interaction ** interactionList;
	inline int min(int a, int b) {if (a < b) return a; else return b;}
	inline int max(int a, int b) {if (a > b) return a; else return b;}
	double offScale, offPhaseFactor;
	double LJumbrellaEnergy;
	vector3d externalDipole, rotationAxis, initialDipole, prevDipole, offsetDipole, initOffDipole, prevOffDipole; bool clockwise;
	bool isEquil, isMin, isMD;
	string name; int runNo, maxRuns;
	vector<pair<int, int> > topSpecList;
	vector<vector<pair<int, int> > > mtopSpecList;
	unordered_set<int> active; bool isActive;
	double sinDipoleAngle; bool circular;
	double loss;
	int intervalCount;
	double sumTKE = 0, sumRKE = 0;
	int totalKEsample = 0;
	double expectedTKE, expectedRKE;
	
	inline void rotateED() {rotateEDPrime(); rotateOD();}

	inline void rotateEDPrime()
	{
		double timeFrac = ((double)currStep)/maxTimeStep + offPhaseFactor;
		if (!sinDipole)
		{
			if (clockwise)
			{
				quaternion q(rotationAxis, (timeFrac - offPhaseFactor)*2*chain*PI*stepMult);
				sinDipoleAngle = (timeFrac - offPhaseFactor)*2*chain*PI*stepMult; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
				prevDipole = externalDipole;
				externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
			}
			else
			{
				quaternion q(rotationAxis, -(timeFrac - offPhaseFactor)*2*chain*PI*stepMult);
				sinDipoleAngle = -(timeFrac - offPhaseFactor)*2*chain*PI*stepMult; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
				prevDipole = externalDipole;
				externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
			}
		}
		else
		{
			double angleConst = PI*sinDipoleFactor, angle = timeFrac*2*chain*PI;
			if (circular)
			{
				double x = timeFrac; if (x > 0.5 && x <= 1) x = 1 - x; if (x > 1) x -= 1;
				//if (printStep) cout << timeFrac << ' ' << x << '\n';
				if (clockwise)
				{
					quaternion q(rotationAxis, (x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor);
					sinDipoleAngle = (x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
					prevDipole = externalDipole;
					externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
				}
				else
				{
					quaternion q(rotationAxis, -(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor);
					sinDipoleAngle = -(x - offPhaseFactor)*chain*PI*stepMult*2*sinDipoleFactor; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
					//if (printStep) cout << "sda: " << sinDipoleAngle << '\n';
					prevDipole = externalDipole;
					externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
				}
			}
			else
			{
				if (clockwise)
				{
					double offset = angleConst*(1 - cos(offPhaseFactor*2*chain*PI))/2;
					quaternion q(rotationAxis, angleConst*(1 - cos(angle))/2 - offset);
					sinDipoleAngle = angleConst*(1 - cos(angle))/2 - offset; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
					prevDipole = externalDipole;
					externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
				}
				else
				{
					double offset = angleConst*(cos(offPhaseFactor*2*chain*PI) - 1)/2;
					quaternion q(rotationAxis, angleConst*(cos(angle) - 1)/2 - offset);
					sinDipoleAngle = angleConst*(cos(angle) - 1)/2 - offset; sinDipoleAngle = (roundf(sinDipoleAngle*10000) + 0.0)/10000;
					prevDipole = externalDipole;
					externalDipole = initialDipole; externalDipole.rotate(q); externalDipole.norm();
				}
			}
		}
	}

	inline void rotateOD()
	{
		double timeFrac = ((double)currStep)/maxTimeStep + offPhaseFactor;
		if (!sinDipole)
		{
			if (clockwise)
			{
				quaternion q(rotationAxis, (timeFrac - offPhaseFactor)*2*chain*PI*stepMult*offScale);
				prevOffDipole = offsetDipole;
				offsetDipole = initOffDipole; offsetDipole.rotate(q); offsetDipole.norm();
			}
			else
			{
				quaternion q(rotationAxis, -(timeFrac - offPhaseFactor)*2*chain*PI*stepMult*offScale);
				prevOffDipole = offsetDipole;
				offsetDipole = initOffDipole; offsetDipole.rotate(q); offsetDipole.norm();
			}
		}
		else
		{
			double angleConst = PI*sinDipoleFactor, angle = timeFrac*2*chain*PI;
			if (clockwise)
			{
				double offset = angleConst*offScale*(1 - cos(offPhaseFactor*2*chain*stepMult*PI))/2;
				quaternion q(rotationAxis, angleConst*(1 - cos(angle))/2*offScale - offset);
				prevOffDipole = offsetDipole;
				offsetDipole = initOffDipole; offsetDipole.rotate(q); offsetDipole.norm();
			}
			else
			{
				double offset = angleConst*offScale*(cos(offPhaseFactor*2*chain*stepMult*PI) - 1)/2;
				quaternion q(rotationAxis, angleConst*(cos(angle) - 1)/2*offScale - offset);
				//cout << angleConst*(cos(angle) - 1)/2*offScale - offset << '\n';
				prevOffDipole = offsetDipole;
				offsetDipole = initOffDipole; offsetDipole.rotate(q); offsetDipole.norm();
			}
		}
	}

	void findObservables()
	{
		observables.clear();
		intervalCount = 0;
		//printObservableList.clear();
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->obs != 0)
			{
				double temp;
				if (interactionList[i]->obs < 0) temp = interactionList[i]->obs + interactionList[i]->phase;
				else temp = interactionList[i]->obs - interactionList[i]->phase;
				temp *= 2*PI;
				observables.push_back(interactionList[i]);
				temp = (roundf(temp*10000) + 0.0)/10000;
				if (triggerObs.count(temp) == 0)
				{
					trobs v; v.vi.push_back(observables.size() - 1);
					triggerObs[temp] = v;
				}
				else
				{
					triggerObs[temp].vi.push_back(observables.size() - 1);
				}
			}
		}
		/*unordered_map<float, trobs>::iterator it;
		cout << "iterator: ";
		for (it = triggerObs.begin(); it != triggerObs.end(); ++it) cout << it->first << ' ' << it->second.vi.size() << '\n';
		cout << "obssize " << observables.size() << '\n';*/
		finalObservables.resize(observables.size());
		poltemp.clear(); poltemp.resize(observables.size());
	}
	
	inline void populateMatrix(int w, int step)
	{
		vector<float> vc; vc.resize(observables.size());
		double umbrella = 0;
		for (int i = 0; i < vc.size(); i++)
		{
			vc[i] = observables[i]->observable();
			umbrella += observables[i]->umbrella();
		}
		if (writePOL)
		{
			if (runNo > 1)
			{
				for (int i = 0; i < vc.size(); i++)
				{
					float k = roundf(vc[i]);
					map<float, int> &m = printObservableList[intervalCount][i];
					if (m.count(k) == 0) m[k] = 1;
					else m[k]++;
				}
				if (w == 0)
				{
					intervalCount++;
				}
			}
			else
			{	
				for (int i = 0; i < vc.size(); i++)
				{
					float k = roundf(vc[i]);
					if (poltemp[i].count(k) == 0) poltemp[i][k] = 1;
					else poltemp[i][k]++;
				}
				if (w == 0)
				{
					printObservableList.push_back(poltemp);
					poltemp.clear(); poltemp.resize(vc.size());
				}
			}
		}
		else
		{
			string s;
			for (int i = 0; i < vc.size(); i++)
			{
				char c = (unsigned char)((int)(roundf(vc[i]))<<24);
				s += c;
			}
			double LJumbExp = 1;
			if (LJumbrellaEnergy > 100) LJumbExp = 0;
			else LJumbExp =  exp(-LJumbrellaEnergy);
			if (observableMatrix.count(s) == 0)
			{
				observableMatrix[s].first = LJumbExp;
				observableMatrix[s].second = exp(-umbrella);
			}
			else observableMatrix[s].first += LJumbExp;
		}	
	}
	
	inline void recordFinalObservables()
	{
		int timeFrac;
		if (chain == 1) timeFrac = currStep;
		else
		{
			int chainLength = maxTimeStep/chain, remainder = currStep%chainLength;
			timeFrac = remainder;
		}
		if (triggerObs.count(sinDipoleAngle) == 1)
		{
			trobs &v = triggerObs[sinDipoleAngle];
			if (v.last == -1 || timeFrac - v.last > 0.05*maxTimeStep)
			{
				v.last = timeFrac;
				for (int i = 0; i < v.vi.size(); i++)
				{
					finalObservables[v.vi[i]].push_back(observables[v.vi[i]]->observable());
				}
			}
		}
	}

	inline void findFinalObservables()
	{
		unordered_map<float, trobs>::iterator it;
		for (it = triggerObs.begin(); it != triggerObs.end(); ++it)
		{
			it->second.last = -1;
		}
		finalObservableList.push_back(finalObservables);
		for (int i = 0; i < finalObservables.size(); i++) finalObservables[i].clear();
	}

	inline void catpInit()
	{
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->type == 5 && interactionList[i]->subtype == 0)
			{
				((constrainAxisTime*)interactionList[i])->init();
			}
		}
		externalDipole = initialDipole; prevDipole = initialDipole; expenditure = 0;
	}

	inline double constrainAxisTimePotential(bool calcForce, double adv, double &expend)
	{
		double total = 0.0; expend = 0.0;
		for (int i = 0; i < constrainAxisTimeList.size(); i++)
		{
			constrainAxisTimeList[i].fix(adv);
			expend += constrainAxisTimeList[i].expend;
			total += constrainAxisTimeList[i].potential();
			if (calcForce) constrainAxisTimeList[i].force();
		}
		return total;
	}

	inline double interactionPotential(bool calcForce, bool printStep, double t)
	{
		double total = 0.0;
		for (int i = 0; i < 6; i++) contribs[i] = 0;
		for (int i = 0; i < intListSize; i++)
		{
			//cout << i << ' ' << intListSize << '\n';
			if (isActive && interactionList[i]->active != 0 && active.count(interactionList[i]->active) == 0) {}
			else
			{	
				if (interactionList[i]->type == 5 && interactionList[i]->subtype == 0)
				{
					((constrainAxisTime*)interactionList[i])->fix(t);
					expenditure += ((constrainAxisTime*)interactionList[i])->expend;
				}
				else
				{
					interactionList[i]->fix(printStep, isMin, isEquil);
				}
				contribs[interactionList[i]->type] += interactionList[i]->potential();
				if (calcForce) interactionList[i]->force();
			}
		}
		for (int i = 0; i < 6; i++) total += contribs[i];
		return total;
	}
	#include "pairwise.h" //Lennard-Jones and electrostatic potentials + cell list
	#include "nativeGene.h"
	
	vector<vector<pair<int, int> > > adjList;
	vector<vector<int> > adjList2;

	inline void decompAll()
	{
		for (int i = 0; i < RBlistSize; i++)
		{
			RBlist[i].decompose();
			RBlist[i].force.assign(0, 0, 0); RBlist[i].torque.assign(0, 0, 0); RBlist[i].resetMoments();
		}
	}

	inline double customInteractions(){return 0.0;} //placeholder for custom interactions

	inline double totalPotential(double t)
	{
		double expend;
		return interactionPotential(0, 0, t) + LJpotential(0) + electrostatic(0) + customInteractions();
	}

	inline void printContributions(double t)
	{
		double expend;
		cout << "constraint: " << contribs[0] << '\n';
		cout << "constrainAxis: " << contribs[1] << '\n';
		cout << "constrainAxisTime: " << constrainAxisTimePotential(0, t, expend) << '\n';
		cout << "bond: " << contribs[2] << '\n';
		cout << "angle: " << contribs[3] << '\n';
		cout << "dihedral: " << contribs[4] << '\n';
		cout << "external dipole: " << contribs[5] << '\n';
		cout << "dispersion: " << LJcontrib << '\n';
		cout << "electrostatic: " << electrostatic(0) << '\n';
	}

	inline double totalForce(double t, bool printStep)
	{
		double expend = 0;
		double iP = interactionPotential(1, printStep, t);
		double eP = electrostatic(1);
		double ljP = LJpotential(1);
		for (int i = 0; i < RBlistSize; i++) RBlist[i].calcTorque();
		expenditure += expend;
		return iP + ljP + eP + customInteractions();
	}

	inline double SP() //Single point energy
	{
		initDecompAll();
		return totalPotential(0);
	}

	void printBonds()
	{
		cout << "bonds: " << bondList.size() << '\n';
		for (int i = 0; i < bondList.size(); i++) bondList[i].printtocout();
	}
	void initialise()
	{
		for (int i = 0; i < intListSize; i++)
		{
			if (interactionList[i]->type == 5 && interactionList[i]->subtype == 3)
			{
				//cout << i << '\n';
				staccatoDipole * sd = (staccatoDipole*)interactionList[i];
				sd->findAxes(sinDipoleFactor, chain, rotationAxis, initialDipole, clockwise);
			}
		}
	}
	//#include "thread.h"
	#include "langevin2.h" //langevin thermostat
	#include "top.h"
	#include "costf.h"
	
	void expectedKEfinder()
	{
		int nonDummies = RBlistSize;
		for (int i = 0; i < RBlistSize; i++) if (RBlist[i].isDummy) nonDummies--;
		expectedTKE = 1.5*nonDummies; expectedRKE = expectedTKE;
	}
	
	void KEfinder()
	{
		expectedKEfinder();
		sumTKE /= totalKEsample; sumRKE /= totalKEsample;
		cout << "Total TKE excess: " << sumTKE << ' ' << sumTKE - expectedTKE << '\n';
		cout << "Total RKE excess: " << sumRKE << ' ' << sumRKE - expectedRKE << '\n';
	}
	
	simulation()
	{
		printStep = 0; expenditure = 0; intListSize = 0; runNo = 1; maxRuns = 1;
		updateCellList = 1; chain = 1;
		assignLJparam();	
		interactionTime = 0, integratorTime = 0, liTime = 0, cellListTime = 0;
		interactionList = NULL; RBlist = NULL;
		initialDipole = defDipole; rotationAxis = defAxis; clockwise = defDirection; 
		initialDipole.norm(); rotationAxis.norm();
		externalDipole = initialDipole; prevDipole = initialDipole;
		offScale = offScaleGen; offPhaseFactor = offPhaseFactorGen;
		initOffDipole = defOffDipole;
		offsetDipole = initOffDipole; prevOffDipole = offsetDipole;
		isEquil = 0; isMin = 0; isMD = 0; isActive = 0; sinDipoleAngle = 0; circular = 0;
	}
	void overrideParameters()
	{
		printStep = 0; expenditure = 0; intListSize = 0; runNo = 1; maxRuns = 1;
		updateCellList = 1; chain = 1;
		assignLJparam();	
		interactionTime = 0, integratorTime = 0, liTime = 0, cellListTime = 0;
		interactionList = NULL; RBlist = NULL;
		initialDipole = defDipole; rotationAxis = defAxis; clockwise = defDirection;
		if (rotateDip != 0)
		{
			quaternion q(rotationAxis, rotateDip);
			initialDipole.rotate(q);
			
		}
		initialDipole.norm(); rotationAxis.norm();
		externalDipole = initialDipole; prevDipole = initialDipole;
		offScale = offScaleGen; offPhaseFactor = offPhaseFactorGen;
		initOffDipole = defOffDipole;
		offsetDipole = initOffDipole; prevOffDipole = offsetDipole;
		isEquil = 0; isMin = 0; isMD = 0; isActive = 0; sinDipoleAngle = 0; circular = Fcircular;
	}
	void deepCopy(simulation &target)
	{
		target = *this;
		target.RBlist = new rigidBody[RBlistSize];
		for (int i = 0; i < RBlistSize; i++)
		{
			target.RBlist[i] = RBlist[i]; target.RBlist[i].repointer(&target.minDim, &target.maxDim, &target.updateCellList);
		}
		target.interactionList = new interaction*[intListSize];
		string fname = "temp"; ofstream ofs(fname.c_str());
		for (int i = 0; i < intListSize; i++)
		{
			int type = stoi(to_string(interactionList[i]->type) + to_string(interactionList[i]->subtype));
			ofs << i << ' ' << type << ' ' << intListSize << '\n';
			switch (type)
			{
				case 20:
					target.interactionList[i] = new bond;
					break;
				case 30:
					target.interactionList[i] = new angle;
					break;
				case 31:
					target.interactionList[i] = new angleAxis;
					break;
				case 40:
					target.interactionList[i] = new dihedral;
					break;
				case 0:
					target.interactionList[i] = new constraint;
					*((constraint*)target.interactionList[i]) = *((constraint*)interactionList[i]);
					break;
				case 10:
					target.interactionList[i] = new constrainAxis;
					*((constrainAxis*)target.interactionList[i]) = *((constrainAxis*)interactionList[i]);
					target.interactionList[i]->redist();
					break;
				case 50:
					target.interactionList[i] = new constrainAxisTime;
					*((constrainAxisTime*)target.interactionList[i]) = *((constrainAxisTime*)interactionList[i]);
					target.interactionList[i]->redist();
					break;
				case 51:
					target.interactionList[i] = new extDipole;
					*((extDipole*)target.interactionList[i]) = *((extDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 52:
					target.interactionList[i] = new offDipole;
					*((offDipole*)target.interactionList[i]) = *((offDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.offsetDipole, &target.prevOffDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
				case 53:
					target.interactionList[i] = new staccatoDipole;
					*((staccatoDipole*)target.interactionList[i]) = *((staccatoDipole*)interactionList[i]);
					target.interactionList[i]->redist();
					target.interactionList[i]->repointer2(&target.externalDipole, &target.prevDipole, &target.expenditure, &target.sinDipoleAngle);
					break;
			}
			if (type > 10 && type < 50) *(target.interactionList[i]) = *(interactionList[i]);
			target.interactionList[i]->repointer(target.RBlist);
			interactionList[i]->print(ofs);
			target.interactionList[i]->print(ofs);
			//cout << "original: " << (this->interactionList[i])->ground << ' '  << (this->interactionList[i])->k << ' ' << (this->interactionList[i])->d << '\n';
			//cout << "copy: " << (target.interactionList[i])->ground << ' '  << (target.interactionList[i])->k << ' ' << (target.interactionList[i])->d << '\n';
		}
	}
	
	~simulation()
	{
		//cout << "deleting\n";
		/*cout << intListSize << ' ' << interactionList << '\n';
		for (int i = 0; i < intListSize; i++) 
		{
			cout << i << ' ' << interactionList[i] << '\n';
			delete interactionList[i];
		}*/
		//delete [] interactionList;
		//cout << "deleted intList\n";
		//delete [] RBlist;
		//cout << "deleted\n";
	}
	#include "nativeParse.h"
};
