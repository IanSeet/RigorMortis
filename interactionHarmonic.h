struct interaction 
{
	double ground, k, d; int active; bool valid, equil, min, umb; short type, subtype = 0; double obs, phase;
	vector<pair<int, int> > px;
	vector<int> rx;
	vector<rigidBody*> rbx;
	vector<vector3d*> vx;
	void resize(int i)
	{
		px.resize(i); rx.resize(i); rbx.resize(i); vx.resize(i);
	}
	void repointer(rigidBody *rb)
	{
		for (int i = 0; i < px.size(); i++)
		{
			rbx[i] = rb + px[i].first;
			if (px[i].second < 0) {vx[i] = &(rbx[i]->center); rx[i] = -1;}
			else {vx[i] = &(rbx[i]->specList[px[i].second].decomp); rx[i] = px[i].second;}
		}
	}
	virtual void redist(){}
	virtual void repointer2(vector3d *_extDip, vector3d *_prevDip, double *_expend, double *_sda){}
	void printPair(ofstream &ofs)
	{
		for (int i = 0; i < px.size(); i++) ofs << px[i].first << ' ' << px[i].second << ' ';
		ofs << '\n';
	}
	void printPairCout()
	{
		for (int i = 0; i < px.size(); i++) cout << px[i].first << ' ' << px[i].second << ' ';
		cout << '\n';
	}
	void parsePair(ifstream &ifs, int RBoffset)
	{
		for (int i = 0; i < px.size(); i++) ifs >> px[i].first >> px[i].second;
		for (int i = 0; i < px.size(); i++)  px[i].first += RBoffset;
	}
	inline virtual void fix(bool printStep, bool isMin, bool isEquil){}
	inline virtual double potential(){return 0;}
	inline virtual void force(){}
	inline virtual float observable(){return 0;}
	inline virtual double umbrella(){return 0;}
	inline virtual void print(ofstream &ofs){}
	inline virtual void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset){}
	inline virtual void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset, vector3d &center, quaternion &rotate){}
};

struct constraint: interaction //distance-based one-particle simple harmonic potential
{
	vector3d locus;
	inline constraint(){resize(1);}
	inline constraint(pair<int, int> &p, vector3d& _locus, double _ground, double _k, rigidBody * RBlist)
	{
		resize(1);
		ground = _ground; k = _k; locus = _locus; px[0] = p; type = 0; obs = 0;
		repointer(RBlist);
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; return;}
		if (equil && !isEquil) {d = 0; return;}
		d = dist(*vx[0], locus) - ground;
	}
	inline double potential()
	{
		return 0.5*k*d*d;
	}
	inline void force()
	{
		vector3d f1 = (locus - *vx[0]); 
		f1.norm(); f1 *= k*d; rbx[0]->force += f1;
		if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << ' ' << ground << ' ' << k << ' ';
		printPair(ofs);
		locus.print(ofs);
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k; parsePair(ifs, RBoffset);
		locus.read(ifs);
		ifs >> equil >> min >> active;
		type = 0; obs = 0;
		repointer(RBlist);
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset, vector3d &center, quaternion &rotate)
	{
		ifs >> ground >> k; parsePair(ifs, RBoffset);
		locus.read(ifs); locus += center;
		ifs >> equil >> min >> active;
		type = 0; obs = 0;
		repointer(RBlist);
	}
};

struct constrainAxis: interaction //constrains object to lie on a given axis
{
	vector3d axis, normv;
	double dist;
	void redist()
	{
		dist = (*vx[0] - rbx[0]->specList[0].decomp).mag();
	}
	inline constrainAxis(){resize(1);}
	inline constrainAxis(pair<int, int> &p, vector3d& _axis, double _ground, double _k, rigidBody * RBlist)
	{
		resize(1);
		ground = _ground; k = _k; axis = _axis; axis.norm(); px[0] = p; type = 1; obs = 0;
		repointer(RBlist); redist();
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 1; valid = 0; return;}
		if (equil && !isEquil) {d = 1; valid = 0; return;}
		normv = (*vx[0] - rbx[0]->specList[0].decomp); normv /= dist; //normv.printtocout();
		d = normv * axis;
		valid = 1;
	}
	inline double potential()
	{
		return k*(1 - d);
	}
	inline void force()
	{
		if (!valid) return;
		vector3d f1 = axis - (normv * d); 
		f1 *= (k / dist);
		rbx[0]->specList[rx[0]].moment += f1;
		rbx[0]->specList[0].moment -= f1;
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << ' ' << ground << ' ' << k << ' '; 
		printPair(ofs);
		axis.print(ofs);
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k; parsePair(ifs, RBoffset); px[0].first += RBoffset;
		axis.read(ifs);
		ifs >> equil >> min >> active;
		type = 1; obs = 0;
		repointer(RBlist); redist();
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset, vector3d &center, quaternion &rotate)
	{
		ifs >> ground >> k; parsePair(ifs, RBoffset);
		axis.read(ifs); axis.rotate(rotate);
		ifs >> equil >> min >> active;
		type = 1; obs = 0;
		repointer(RBlist); redist();

	}
};

struct constrainAxisTime: interaction //constrains object to lie on a given axis (time dependent)
{
	vector3d startAxis, endAxis, diff, currAxis, normv;
	double dist, expend, t;
	double duration, startTime, offset, on, off;
	double duration0, on0, off0, offset0;
	void redist()
	{
		dist = (*vx[0] - rbx[0]->specList[0].decomp).mag();
	}
	inline constrainAxisTime(){resize(1);}
	inline constrainAxisTime(pair<int, int> &p, vector3d& _startAxis, vector3d& _endAxis, double _ground, double _k, double _duration,
	 double _startTime, double _offset, double _on, double _off, rigidBody * RBlist)
	{
		resize(1);
		ground = _ground; k = _k; obs = 0;
		startAxis = _startAxis; startAxis.norm(); type = 5; subtype = 0;
		endAxis = _endAxis; endAxis.norm();
		diff = endAxis - startAxis; currAxis = startAxis;
		px[0] = p; duration = _duration; startTime = _startTime;
		repointer(RBlist); redist();
		offset = _offset; on = _on; off = _off; if (off < 0) off = 2147483647;
		expend = 0; t = 0;
	}
	inline void init()
	{
		t = startTime;
	}
	inline void adjust(int steps)
	{
		if (duration <= 1)
		{
			duration0 = duration*steps;
		}	
		if (offset <= 1)
		{
			offset0 = offset*steps;
		}		
		if (on <= 1)
		{
			on0 = on*steps;
		}
		if (off <= 1)
		{
			off0 = off*steps;
		}
	}
	inline void fix(double adv)
	{
		currAxis = startAxis + (diff * t); currAxis.norm();
		normv = (*vx[0] - rbx[0]->specList[0].decomp); normv /= dist;
		d = normv * currAxis;
		double P = k*(1 - d);
		t = (adv - offset0)/duration0; if (t < 0) t = 0; if (t > 1) t = 1;
		//cout << t << '\n';
		currAxis = startAxis + (diff * t); currAxis.norm();
		normv = (*vx[0] - rbx[0]->specList[0].decomp); normv /= dist;
		if (on0 <= adv && off0 >= adv)
		{
			//cout << on << ' ' << off << '\n';
			d = normv * currAxis;
			expend = k*(1 - d) - P;
		}
		else
		{
			d = 2; expend = 0;
		}
	}
	inline double potential()
	{
		if (d >= -1 && d <= 1) return k*(1 - d);
		else return 0;
	}
	inline void force()
	{
		if (d >= -1 && d <= 1)
		{
			vector3d f1 = currAxis - (normv * d); 
			f1 *= (k / dist);
			rbx[0]->specList[rx[0]].moment += f1;
			rbx[0]->specList[0].moment -= f1;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << subtype << ' ' << ground << ' ' << k << ' ';
		printPair(ofs);
		ofs << duration << ' ' << startTime << ' ' << offset << ' ' << on << ' ' << off << '\n';
		startAxis.print(ofs); endAxis.print(ofs);
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> duration >> startTime >> offset >> on >> off; px[0].first += RBoffset;
		startAxis.read(ifs); endAxis.read(ifs); type = 5; subtype = 0;
		ifs >> equil >> min >> active;
		diff = endAxis - startAxis; currAxis = startAxis;
		repointer(RBlist); redist();
		if (off < 0) off = 2147483647;
		expend = 0; t = 0;
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset, vector3d &center, quaternion &rotate)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> duration >> startTime >> offset >> on >> off; px[0].first += RBoffset;
		startAxis.read(ifs); endAxis.read(ifs); type = 5; subtype = 0;
		startAxis.rotate(rotate); endAxis.rotate(rotate);
		ifs >> equil >> min >> active;
		diff = endAxis - startAxis; currAxis = startAxis;
		repointer(RBlist); redist();
		if (off < 0) off = 2147483647;
		expend = 0; t = 0;
	}
};

struct bond : interaction //distance-based pairwise simple harmonic potential
{
	inline bond(){resize(2);}
	bool isChain = 0;
	inline bond(pair<int, int> &_pa0, pair<int, int> &_pa1, double _ground, double _k, rigidBody * RBlist)
	{
		resize(2);
		px[0] = _pa0; px[1] = _pa1; ground = _ground; k = _k; type = 2; obs = 0;
		repointer(RBlist);
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; valid = 0; return;}
		if (equil && !isEquil) {d = 0; valid = 0; return;}
		valid = 1;
		d = dist(*vx[0], *vx[1]) - ground; //Have to dereference because my function for distance takes references, not pointers. Efficiency should be unaffected on most compilers.
	}
	inline double potential()
	{
		return 0.5*k*d*d;
	}
	inline void force()
	{
		if (!valid) return;
		vector3d f1 = (*vx[1] - *vx[0]), f2; 
		f1.norm(); f1 *= k*d; rbx[0]->force += f1;
		if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
		f2 = (*vx[0] - *vx[1]);
		f2.norm(); f2 *= k*d; rbx[1]->force += f2;
		if (rx[1] != -1) rbx[1]->specList[rx[1]].moment += f2;
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << ' ' << ground << ' ' << k << ' ';
		printPair(ofs);
		ofs << equil << ' ' << min << ' ' << active << ' ' << isChain << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k;
		parsePair(ifs, RBoffset);
		ifs >> equil >> min >> active >> isChain;
		type = 2; obs = 0;
		repointer(RBlist);
	}
	inline void printtocout()
	{
		cout << type << ' ' << ground << ' ' << k << ' ';
		printPairCout();
		cout << equil << ' ' << min << ' ' << active << ' ' << isChain << '\n';
	}
};

struct angle : interaction //harmonic potential for angles
{
	vector3d a1, a2;
	double a1mag, a2mag, d2, dot;
	inline angle(){resize(3);}
	inline angle(pair<int, int> &_pa0, pair<int, int> &_pa1, pair<int, int> &_pa2, double _ground, double _k, rigidBody * RBlist)
	{
		resize(3);
		px[0] = _pa0; px[1] = _pa1; px[2] = _pa2; ground = _ground*PI/180; k = _k; type = 3; subtype = 0;
		repointer(RBlist);
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; valid = 0; return;}
		if (equil && !isEquil) {d = 0; valid = 0; return;}
		valid = 1;
		a1 = *vx[0] - *vx[1], a2 = *vx[2] - *vx[1];
		a1mag = a1.mag(); a2mag = a2.mag();
		if (a1mag == 0 || a2mag == 0) {valid = 0; return;}
		a1 /= a1mag; a2 /= a2mag;
		dot = a1 * a2;
		if (dot > 1) dot = 1; if (dot < -1) dot = -1;
		d = ground - acos(dot);
		if (obs == 1 && printStep && writetocout) 
		{
			cout << "Angle between:\n" <<
			px[0].first << ' ' << px[0].second << ' ' << px[1].first << ' ' << px[1].second << ' ' << px[2].first << ' ' << px[2].second << ": " << -d*180/PI << '\n';
		}
	}
	inline double potential()
	{
		if (!valid) return 0;
		return k*d*d/2;
	}
	inline float observable()
	{
		return (acos(dot)*180/PI + 0.5);
	}
	inline double umbrella()
	{
		double x = (double)(unsigned char)(acos(dot)*180/PI + 0.5);
		x *= PI/180;
		return k*(ground - x)*(ground - x)/2;
	}
	inline void force()
	{
		if (!valid) return;
		vector3d f1, f2, f3;
		double derivative = -k*d, sint = sqrt(1 - dot*dot);
		if (sint != 0)
		{
			derivative /= sint;
			f1 = a2 - a1 * dot; f1 *= (derivative / a1mag); rbx[0]->force += f1;
			if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
			f3 = a1 - a2 * dot; f3 *= (derivative / a2mag); rbx[2]->force += f3;
			if (rx[2] != -1) rbx[2]->specList[rx[2]].moment += f3;
			f2 = f1 + f3; f2 *= -1; rbx[1]->force += f2;
			if (rx[1] != -1) rbx[1]->specList[rx[1]].moment += f2;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << subtype << ' '  << ground << ' ' << k << ' ';
		printPair(ofs);
		ofs << obs << ' ' << umb << ' ' << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k;
		parsePair(ifs, RBoffset);
		ifs >> obs >> umb >> equil >> min >> active;
		type = 3; subtype = 0;
		repointer(RBlist);
	}
	inline void printtocout()
	{
		cout << type << subtype << ' '  << ground << ' ' << k << ' ';
		printPairCout();
		cout << obs << ' ' << umb << ' ' << equil << ' ' << min << ' ' << active << '\n';
	}
};

struct gaussAngle : interaction //harmonic potential for angles
{
	vector3d a1, a2;
	double a1mag, a2mag, d2, dot;
	inline gaussAngle(){resize(3);}
	inline gaussAngle(pair<int, int> &_pa0, pair<int, int> &_pa1, pair<int, int> &_pa2, double _ground, double _k, rigidBody * RBlist)
	{
		resize(3);
		px[0] = _pa0; px[1] = _pa1; px[2] = _pa2; rbx[0] = &RBlist[px[0].first]; rbx[1] = &RBlist[px[1].first]; rbx[2] = &RBlist[px[2].first]; ground = _ground*PI/180; k = _k; type = 3; subtype = 0;
		repointer(RBlist);
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; valid = 0; return;}
		if (equil && !isEquil) {d = 0; valid = 0; return;}
		valid = 1;
		a1 = *vx[0] - *vx[1], a2 = *vx[2] - *vx[1];
		a1mag = a1.mag(); a2mag = a2.mag();
		if (a1mag == 0 || a2mag == 0) {valid = 0; return;}
		a1 /= a1mag; a2 /= a2mag;
		dot = a1 * a2;
		if (dot > 1) dot = 1; if (dot < -1) dot = -1;
		d = ground - acos(dot);
		if (obs == 1 && printStep && writetocout) 
		{
			cout << "Angle between:\n" <<
			px[0].first << ' ' << px[0].second << ' ' << px[1].first << ' ' << px[1].second << ' ' << px[2].first << ' ' << px[2].second << ": " << -d*180/PI << '\n';
		}
	}
	inline double potential()
	{
		if (!valid) return 0;
		return k*exp(-d*d);
	}
	inline float observable()
	{
		return (acos(dot)*180/PI + 0.5);
	}
	inline double umbrella()
	{
		double x = (double)(unsigned char)(acos(dot)*180/PI + 0.5);
		x *= PI/180;
		return k*exp(-(ground - x)*(ground - x));
	}
	inline void force()
	{
		if (!valid) return;
		vector3d f1, f2, f3;
		double derivative = 2*d*k*exp(-d*d), sint = sqrt(1 - dot*dot);
		if (sint != 0)
		{
			derivative /= sint;
			f1 = a2 - a1 * dot; f1 *= (derivative / a1mag); rbx[0]->force += f1;
			if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
			f3 = a1 - a2 * dot; f3 *= (derivative / a2mag); rbx[2]->force += f3;
			if (rx[2] != -1) rbx[2]->specList[rx[2]].moment += f3;
			f2 = f1 + f3; f2 *= -1; rbx[1]->force += f2;
			if (rx[1] != -1) rbx[1]->specList[rx[1]].moment += f2;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << subtype << ' '  << ground;
		printPair(ofs);
		ofs << obs << ' ' << umb << ' ' << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k;
		parsePair(ifs, RBoffset);
		ifs >> obs >> umb >> equil >> min >> active;
		type = 3; subtype = 0;
		repointer(RBlist);
	}
	inline void printtocout()
	{
		cout << type << subtype << ' '  << ground;
		printPairCout();
		cout << obs << ' ' << umb << ' ' << equil << ' ' << min << ' ' << active << '\n';
	}
};

struct angleAxis : interaction
{

	vector3d a1, a2;
	double a1mag, a2mag, d2, dot;
	inline angleAxis(){resize(3);}
	inline angleAxis(pair<int, int> &_pa0, pair<int, int> &_pa1, vector3d &_axis, double _ground, double _k, rigidBody * RBlist)
	{
		resize(3);
		px[0] = _pa0; px[1] = _pa1; ground = _ground*PI/180; k = _k; a2 = _axis; a2mag = a2.mag(); type = 3; subtype = 1;
		repointer(RBlist);
		if (a2mag != 0) a2 /= a2mag;
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; valid = 0; return;}
		if (equil && !isEquil) {d = 0; valid = 0; return;}
		valid = 1;
		a1 = *vx[0] - *vx[1];
		a1mag = a1.mag();
		if (a1mag == 0 || a2mag == 0) {valid = 0; return;}
		a1 /= a1mag;
		dot = a1 * a2;
		if (dot > 1) dot = 1; if (dot < -1) dot = -1;
		d = ground - acos(dot);
		if (obs == 1 && printStep && writetocout) 
		{
			cout << "Angle between:\n" <<
			px[0].first << ' ' << px[0].second << ' ' << px[1].first << ' ' << px[1].second << " and " << a2.i << ' ' << a2.j << ' ' << a2.k << " : " << -d*180/PI << '\n';
		}
	}
	inline double potential()
	{
		if (!valid || obs) return 0;
		return k*d*d/2;
	}
	inline float observable()
	{
		return (acos(dot)*180/PI + 0.5);
	}
	inline double umbrella()
	{
		double x = (double)(unsigned char)(acos(dot)*180/PI + 0.5);
		x *= PI/180;
		return k*(ground - x)*(ground - x)/2;
	}
	inline void force()
	{
		if (!valid || obs) return;
		vector3d f1, f2;
		double derivative = -k*d, sint = sqrt(1 - dot*dot);
		if (sint != 0)
		{
			derivative /= sint;
			f1 = a2 - a1 * dot; f1 *= (derivative / a1mag); rbx[0]->force += f1;
			if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
			f2 = a1 - a2 * dot; f2 *= (derivative / a2mag); rbx[2]->force += f2;
			if (rx[1] != -1) rbx[1]->specList[rx[1]].moment += f2;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << subtype << ' '  << ground << ' ' << k << ' ' << px[0].first << ' ' << px[0].second << ' ' << px[1].first << ' ' << px[1].second << '\n';
		ofs << obs << ' ' << umb << ' ' << equil << ' ' << min << ' ' << active << ' ' << phase << '\n';
		a2.print(ofs);
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, int RBoffset)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> px[1].first >> px[1].second >> obs >> umb >> equil >> min >> active >> phase;
		//cout << "obs: " << obs << '\n';
		px[0].first += RBoffset; px[1].first += RBoffset;
		a2.read(ifs);
		a2mag = a2.mag(); type = 3; subtype = 1;
		repointer(RBlist);
		if (a2mag != 0) a2 /= a2mag;
	}
};

struct dihedral : interaction //periodic potential for dihedrals
{
	vector3d a1, a2, a3, n1, n2;
	double d2, n1mag, n2mag, dot;
	inline dihedral(){resize(4);}
	inline dihedral(pair<int, int> &_pa0, pair<int, int> &_pa1, pair<int, int> &_pa2, pair<int, int> &_pa3, double _ground, double _k, double _divisor, rigidBody * RBlist)
	{
		resize(4);
		px[0] = _pa0; px[1] = _pa1; px[2] = _pa2; px[3] = _pa3; type = 4; obs = 0;
		ground = _ground*PI/180; k = _k;
		repointer(RBlist);
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		if (min && !isMin) {d = 0; valid = 0; return;}
		if (equil && !isEquil) {d = 0; valid = 0; return;}
		valid = 1;
		a1 = *vx[0] - *vx[1], a2 = *vx[3] - *vx[2], a3 = *vx[2] - *vx[1]; a3.norm();
		n1 = a1 - a3 * (a1 * a3); n1mag = n1.mag(); if (n1mag == 0) {valid = 0; return;} n1 /= n1mag;
		n2 = a2 - a3 * (a2 * a3); n2mag = n2.mag(); if (n2mag == 0) {valid = 0; return;} n2 /= n2mag;
		dot = n1 * n2;
		if (dot > 1) dot = 1; if (dot < -1) dot = -1;
		d = ground - acos(dot); //a1.printtocout(); a2.printtocout(); a3.printtocout(); n1.printtocout(); n2.printtocout(); cout << ground << ' ' << acos(n1 * n2)*180/PI << '\n';
	}
	inline double potential()
	{
		if (!valid) return 0;
		return k*d*d/2;
	}
	inline void force()
	{
		if (!valid) return;
		vector3d f1, f2, f3, f4;
		double derivative = -k*d, sint = sqrt(1 - dot*dot);
		if (sint != 0)
		{
			derivative /= sint;
			f1 = n2 - n1 * dot; f1 *= (derivative / n1mag); rbx[0]->force += f1;
			if (rx[0] != -1) rbx[0]->specList[rx[0]].moment += f1;
			f2 = f1; f2 *= -1; rbx[1]->force += f2;
			if (rx[1] != -1) rbx[1]->specList[rx[1]].moment += f2;
			f4 = n1 - n2 * dot; f4 *= (derivative / n2mag); rbx[3]->force += f4;
			if (rx[3] != -1) rbx[3]->specList[rx[3]].moment += f4;
			f3 = f4; f3 *= -1; rbx[2]->force += f3;
			if (rx[2] != -1) rbx[2]->specList[rx[2]].moment += f3;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << ' ' << ground << ' ' << k << ' ';
		printPair(ofs);
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody *RBlist, int RBoffset)
	{
		ifs >> ground >> k; parsePair(ifs, RBoffset);
		ifs >> equil >> min >> active;
		type = 4; obs = 0;
		repointer(RBlist);
	}
	inline void printtocout()
	{
		cout << type << ' ' << ground << ' ' << k << ' ';
		printPairCout();
		cout << equil << ' ' << min << ' ' << active << '\n';
	}
};

struct extDipole : interaction
{
	vector3d normv, *externalDipole, *prevDipole;
	double dist, *expenditure, expend, deflection;
	quaternion q, qrev;
	void redist()
	{
		dist = (*vx[0] - rbx[0]->specList[0].decomp).mag();
	}
	inline extDipole(){resize(1);}
	inline extDipole(pair<int, int> &p, vector3d *_externalDipole, vector3d *_prevDipole, double *_expenditure, double _ground, double _k, rigidBody * RBlist, 
	double _deflection, vector3d &rotationAxis)
	{
		resize(1);
		ground = _ground; k = _k; px[0] = p; type = 5; obs = 0; externalDipole = _externalDipole; expenditure = _expenditure; prevDipole = _prevDipole;
		subtype = 1; deflection = _deflection;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
	}
	virtual void repointer2(vector3d *_extDip, vector3d *_prevDip, double *_expend, double *_sda)
	{
		externalDipole = _extDip; prevDipole = _prevDip; expenditure = _expend;
	}
	inline virtual void fix(bool printStep, bool isMin, bool isEquil)
	{
		normv = (*vx[0] - rbx[0]->specList[0].decomp);
		if (deflection != 0) normv.rotate(q);
		normv /= dist;
		d = normv * (*externalDipole);
	}
	inline double potential()
	{
		return k*(1 - d);
	}
	inline virtual void force()
	{
		expend = k*(normv * (*prevDipole) - d);
		*expenditure += expend;
		vector3d f1;
		if (!expFront) f1 = *externalDipole - (normv * d);
		else f1 = *prevDipole - (normv * d);
		//f1.printtocout(); normv.printtocout();
		if (deflection != 0) f1.rotate(qrev);
		f1 *= (k / dist);
		rbx[0]->specList[rx[0]].moment += f1;
		rbx[0]->specList[0].moment -= f1;
	}
	inline virtual void print(ofstream &ofs)
	{
		ofs << type << subtype << ' ' << ground << ' ' << k << ' ' << px[0].first << ' ' << px[0].second << ' ' << deflection << '\n';
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline virtual void parse(ifstream &ifs, rigidBody * RBlist, vector3d * _externalDipole, vector3d *_prevDipole, double *_expenditure, int RBoffset, vector3d &rotationAxis)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> deflection >> equil >> min >> active;
		px[0].first += RBoffset;
		type = 5; obs = 0; externalDipole = _externalDipole; expenditure = _expenditure; prevDipole = _prevDipole;
		subtype = 1;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
	}
};

struct offDipole : extDipole
{
	inline offDipole(){resize(1);}
	inline offDipole(pair<int, int> &p, vector3d *_offDipole, vector3d *_prevOffDipole, double *_expenditure, double _ground, double _k, rigidBody * RBlist,
	double _deflection, vector3d &rotationAxis)
	{
		resize(1);
		ground = _ground; k = _k; px[0] = p; type = 5; obs = 0; externalDipole = _offDipole; expenditure = _expenditure; prevDipole = _prevOffDipole;
		subtype = 2; deflection = _deflection;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, vector3d * _externalDipole, vector3d *_prevDipole, double *_expenditure, int RBoffset, vector3d &rotationAxis)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> deflection >> equil >> min >> active;
		px[0].first += RBoffset;
		type = 5; obs = 0; externalDipole = _externalDipole; expenditure = _expenditure; prevDipole = _prevDipole;
		subtype = 2;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
	}
};

struct staccatoDipole : extDipole
{
	vector3d normv, startAxis, haltAxis;
	double start, halt, silenceStart, silenceEnd;
	float dir; double *sinDipoleAngle;
	inline staccatoDipole(){resize(1);}
	inline staccatoDipole(pair<int, int> &p, vector3d *_externalDipole, vector3d *_prevDipole, double *_expenditure, double _ground, double _k, rigidBody * RBlist, 
	double _deflection, vector3d &rotationAxis, double _start, double _halt, double _silenceStart, double _silenceEnd, double * _sda)
	{
		resize(1);
		ground = _ground; k = _k; px[0] = p; type = 5; obs = 0; externalDipole = _externalDipole; expenditure = _expenditure; prevDipole = _prevDipole;
		subtype = 3; deflection = _deflection;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
		start = _start; halt = _halt; silenceStart = _silenceStart; silenceEnd = _silenceEnd; sinDipoleAngle = _sda;
		start *= 2*PI; halt *= 2*PI; silenceStart *= 2*PI; silenceEnd *= 2*PI;
	}
	void repointer2(vector3d *_extDip, vector3d *_prevDip, double *_expend, double *_sda)
	{
		externalDipole = _extDip; prevDipole = _prevDip; expenditure = _expend; sinDipoleAngle = _sda;
	}
	void findAxes(double sinDipoleFactor, double chain, vector3d &rotationAxis, vector3d &initialDipole, bool clockwise)
	{
		double angleConst = PI*sinDipoleFactor, angle = start*chain, angle2 = halt*chain;
		quaternion q(rotationAxis, angle);
		startAxis = initialDipole; startAxis.rotate(q); startAxis.norm();
		quaternion q2(rotationAxis, angle2);
		haltAxis = initialDipole; haltAxis.rotate(q2); haltAxis.norm();
		if (halt < start) dir = -1; else dir = 1;
	}
	inline void fix(bool printStep, bool isMin, bool isEquil)
	{
		normv = (*vx[0] - rbx[0]->specList[0].decomp);
		if (deflection != 0) normv.rotate(q);
		normv /= dist;
		if (silenceStart < *sinDipoleAngle && silenceEnd >= *sinDipoleAngle)
		{
			d = 1;
		}
		else
		{
			if (start*dir > *sinDipoleAngle*dir)
			{
				d = normv * startAxis;
			}
			else if (halt*dir < *sinDipoleAngle*dir)
			{
				d = normv * haltAxis;
			}
			else
			{
				d = normv * (*externalDipole);
			}
		}
	}
	inline void force()
	{
		expend = 0;
		vector3d f1;
		if (!(silenceStart < *sinDipoleAngle && silenceEnd >= *sinDipoleAngle))
		{
			if (start*dir > *sinDipoleAngle*dir)
			{
				f1 = startAxis - (normv * d);
			}
			else if (halt*dir < *sinDipoleAngle*dir)
			{
				f1 = haltAxis - (normv * d);
			}
			else
			{
				expend = k*(normv * (*prevDipole) - d);
				*expenditure += expend;
				if (!expFront) f1 = *externalDipole - (normv * d);
				else f1 = *prevDipole - (normv * d);
			}
			if (deflection != 0) f1.rotate(qrev);
			f1 *= (k / dist);
			rbx[0]->specList[rx[0]].moment += f1;
			rbx[0]->specList[0].moment -= f1;
		}
	}
	inline void print(ofstream &ofs)
	{
		ofs << type << subtype << ' ' << ground << ' ' << k << ' ' << px[0].first << ' ' << px[0].second << ' ' << deflection << ' ' << start << ' ' << halt << ' ' << silenceStart <<
		' ' << silenceEnd << '\n';
		ofs << equil << ' ' << min << ' ' << active << '\n';
	}
	inline void parse(ifstream &ifs, rigidBody * RBlist, vector3d * _externalDipole, vector3d *_prevDipole, double *_expenditure, int RBoffset, vector3d &rotationAxis, double *_sda)
	{
		ifs >> ground >> k >> px[0].first >> px[0].second >> deflection >> start >> halt >> silenceStart >> silenceEnd >> equil >> min >> active;
		px[0].first += RBoffset;
		type = 5; obs = 0; externalDipole = _externalDipole; expenditure = _expenditure; prevDipole = _prevDipole;
		rx[0] = px[0].second; subtype = 3;
		repointer(RBlist); redist();
		q = convert(rotationAxis, deflection);
		qrev = convert(rotationAxis, -deflection);
		sinDipoleAngle = _sda;
	}
};
