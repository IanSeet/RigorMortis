#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <pthread.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
using namespace std;
bool writetocout = 1;
const double PI = 3.14159265359, C = 112.1483235, D = 1.0, CAdj = C/D, Big = 2e128, Smol = -2e128, sqrt3 = sqrt(3), intmax = 2147483647; //Pi, coulomb and dielectric constants
double h = 1e-3, h2, mu = 1.5, minh = 1e-4; //step size for integrator
const double scale = 0.1;
double LJcutoff = 4.0, verletRadius = 4.3, verletSkin = (verletRadius - LJcutoff)/2; //LJ potential cutoff distance
bool wca = 0, noise = 1, multi = 0, roteq = 0, dipoleOn = 1, sinDipole = 0, defDirection = 1, expSum = 0; //default direction is anticlockwise. expSum = 1 when expenditures are summed across each writeout
int minsteps = 20000, eqsteps = 20000, roteqsteps = 1e5, cntOffset = 0;
double massMult = 1, kT = 1.0, stepMult = 1; //massMult multiplies mass and stepMult multiplies the timestep by the specified constant
double wcaCutoff = 1.6, wcaRadius = 1.9, wcaSkin = (wcaRadius - wcaCutoff)/2;  //WCA potential constants
bool expFront = 0; double dipoleMass = 0; //these were for debugging the external dipole, do not modify
double sinDipoleFactor = 1, rotateDip = 0, takeObs = 0;
#include "quaternion.h" //Quaternion and 3d vector functions
vector3d defDipole, defOffDipole, defAxis;
double offScaleGen = 1.0, offPhaseFactorGen = 0.0; bool writePOL = 0;
double phaseObs = 0.0; bool inactivateDum = 1, Fcircular = 0;
bool isTest = 0, redirectCout = 0;
#include "Eigen/Eigenvalues" //Finds eigenvalues and eigenvectors of a hermitian 3x3 matrix
using namespace Eigen;
#include "rigidBody.h" //Rigid body functions
rigidBody * RBlistMain;
#include "interactionHarmonic.h" //Constraint, bond, angle and dihedral interactions
struct triple
{
	int a, b, c;
	inline triple(){}
	inline triple(int _a, int _b, int _c){a = _a; b = _b; c = _c;}
};

struct genAlgoParam
{
	double mutrate;
	int maxthreads;
	int steps, runs, eqsteps, roteqsteps;
	genAlgoParam();
	genAlgoParam(double _mutrate, int _maxthreads, int _steps, int _runs, int _eqsteps, int _roteqsteps):
	mutrate(_mutrate), maxthreads(_maxthreads), steps(_steps), runs(_runs), eqsteps(_eqsteps), roteqsteps(_roteqsteps){}
};

struct gene
{
	vector<pair<int, int> > vecParticle;
	vector<int> vecInt;
	int type;
	double shift;
	vector3d disp;
	gene(){shift = 1; disp.assign(0, 0, 0);}
};

#include "simulation.h"
vector<simulation> vecsim;
#include "gene.h"
#include "parser.h" //Parses input and writes to output
#include "multi.h"

int main(int argc, char *argv[])
{
	string s; int rng = 0, multi = 0; if (argc >= 2) s = argv[1];
	if (!multi && s != "unitTest")
	{
		int seed = time(NULL) + rng;
		srand(seed); //cout << seed << '\n';
	}
	if (s == "unitTest") {srand(0); isTest = 1;}
	if (isTest) redirectCout = 1;
	ofstream coutput("cout");
	streambuf *coutbuf = cout.rdbuf(); //save old buf
	if (redirectCout)
	{
		cout.rdbuf(coutput.rdbuf());
	}
	ifstream ifs(s.c_str());
	double t = clock();
	if (multi) parseInputMulti(ifs);
	else 
	parseInput(ifs);
	if (redirectCout)
	{
		cout.rdbuf(coutbuf);
	}
	
	for (int i = 0; i < vecsim.size(); i++)
	cout << i << " integrator time: " << vecsim[i].integratorTime << " interaction time: " << vecsim[i].interactionTime << " li time: " << vecsim[i].liTime << ' ' << "cell list time: " << vecsim[i].cellListTime << '\n';
	cout << "time to completion: " << clock() - t << '\n';
}
