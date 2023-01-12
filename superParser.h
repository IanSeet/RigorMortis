void parseRigid(ifstream &ifs, string &s, vector3d &offset, int counter, quaternion &orient3, simulation &sim, int idx)
{
	rigidBody * RBlist = sim.RBlist;
	vector3d empty(0, 0, 0), omega;
	double a, b, c, d, m; string sp; int col = 0; bool isavrg = 0, overrideOmega = 0;
	getline(ifs, s); istringstream iss(s);
	quaternion orient;
	vector3d center, euler;
	iss >> s;
	if (s == "avr")
	{
		isavrg = 1;
		int a, b, c, d;
		iss >> a >> b >> c >> d; a += counter; c += counter;
		RBlist[a].decompose(); RBlist[c].decompose();
		vector3d v = RBlist[a].specList[b + 1].decomp + RBlist[c].specList[d + 1].decomp; v /= 2; v -= offset;
		center.assign(v);
		euler.assign(empty);
	}
	else
	{
		a = atof(s.c_str());
		iss >> b >> c; center.assign(a, b, c);
		iss >> a >> b >> c; euler.assign(a, b, c);
	}
	iss >> col;
	sim.colour.push_back(col);
	orient = eulerConv(euler);
	vector<vector3d> vspecial; vspecial.push_back(empty);
	vector<pair<double, vector3d> > charges;
	vector<cuboid> components;
	vector<char> ptype;
	while (true)
	{
		getline(ifs, s); if (s[0] == 'r' || s[0] == '!' || s[0] == 'p' || s[0] == 'd') break;
		if (s[0] != ';')
		{
			//cout << s << '\n';
			m = 1;
			istringstream iss(s);
			iss >> s;
			quaternion orient2;
			vector3d dim, center2, offset2(0, 0, 0); d = 1;
			if (s == "def")
			{
				ptype.push_back('A');
				iss >> a >> b >> c; dim.assign(a, b, c);
				iss >> a >> b >> c; center2.assign(a, b, c);
				iss >> sp; if (sp == "spacing") iss >> d;
				iss >> sp; if (sp == "mass") iss >> m;
				cuboid c(dim, center2, offset2, orient2, d, m);
				components.push_back(c);
			}
			else if (s == "sp")
			{
				iss >> a >> b >> c; dim.assign(a, b, c);
				vspecial.push_back(dim);
			}
			else if (s == "cg")
			{
				iss >> a >> b >> c >> d; dim.assign(a, b, c); pair<double, vector3d> p(d, dim); charges.push_back(p);
			}
			else if (s == "omega")
			{
				iss >> a >> b >> c; omega.assign(a, b, c); overrideOmega = 1;
			}
			else
			{
				ptype.push_back(s[0]); m = sim.paramTable[s[0]].mass;
				iss >> a >> b >> c; dim.assign(a, b, c);
				iss >> a >> b >> c; center2.assign(a, b, c);
				//iss >> a >> b >> c; vector3d euler2(a, b, c); orient2 = eulerConv(euler2);
				iss >> sp;
				if (sp == "orient")
				{
					iss >> a >> b >> c; vector3d euler2(a, b, c); orient2 = eulerConv(euler2);
				}
				else if (sp == "spacing") iss >> d;
				else if (sp == "mass") iss >> m;
				cuboid c(dim, center2, offset2, orient2, d, m);
				components.push_back(c);
			}
			//dim.printtocout();
		}
	}
	if (!isavrg)
	{
		center.rotate(orient3);
		orient = orient3 * orient;
	}
	center += offset;
	rigidBody rb(components, orient, center, ptype, &sim.minDim, &sim.maxDim, &sim.updateCellList);
	//cout << rb.massListSize << '\n';
	rb.specListSize = vspecial.size();
	//cout << idx << ' ' << rb.specListSize << '\n';
	rb.specList.resize(rb.specListSize);
	for (int i = 0; i < rb.specListSize; i++) rb.specList[i].initial = vspecial[i];
	rb.chargeListSize = charges.size();
	rb.chargeList.resize(rb.chargeListSize);
	for (int i = 0; i < rb.chargeListSize; i++)
	{
		rb.chargeList[i].initial = charges[i].second;
		rb.chargeList[i].charge = charges[i].first;
	}
	rb.decompose();
	if (overrideOmega)
	{
		rb.overrideAngmom = 1;
		rb.overAngmom = omega;
	}
	//cout << RBlist << '\n';
	RBlist[idx] = rb;
}

void parsePoint(ifstream &ifs, string &s, vector3d &offset, simulation &sim, int idx)
{
	vector3d empty(0, 0, 0);
	rigidBody * RBlist = sim.RBlist;
	double a, b, c, d, m; string sp;
	quaternion orient; vector3d center;
	istringstream iss(s); iss >> s;
	iss >> a >> b >> c; center.assign(a, b, c);
	vector<cuboid> components;
	vector<char> ptype; ptype.push_back('C');
	d = 1; m = 0.12;
	cuboid cb(empty, empty, empty, orient, d, m);
	components.push_back(cb);
	center += offset;
	rigidBody rb(components, orient, center, ptype, &sim.minDim, &sim.maxDim, &sim.updateCellList);
	rb.specList.resize(1); rb.specListSize = 1; rb.specList[0].initial.assign(0, 0, 0);
	rb.decompose();
	RBlist[idx] = rb;
}

void parseInteraction(ifstream &ifs, string &s, vector3d &offset, int counter, quaternion &rotate, simulation &sim, int &intIdx, int steps)
{
	rigidBody * RBlist = sim.RBlist;
	istringstream iss(s);
	//cout << s << '\n';
	string aux, aux2; float observe = 0; bool min = 0, equil = 0, umb = 0;
	double ground = 0, k = 0, divisor;
	int active = 0;
	vector<pair<int, int> > vp(4);
	iss >> s >> aux;
	if (aux == "obs") {observe = 1;}
	else if (aux == "obsv") {iss >> observe;}
	else
	{
		if (aux == "act")
		{
			iss >> active >> aux;
			sim.isActive = 1;
		}
		ground = atof(aux.c_str());
		iss >> k;
		if (!iss.eof())
		{
			iss >> aux2;
			//cout << aux2 << '\n';
			if (aux2 == "min") min = 1;
			else if (aux2 == "equil") equil = 1;
			else if (aux2 == "umb") {observe = 1; umb = 1;}
		}
	}
	if (s == "constraint")
	{
		string s2;
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> s2;
		if (s2 == "curr")
		{
			iss2 >> vp[0].first >> vp[0].second;
			vp[0].second++; vp[0].first += counter;
			locus = RBlist[vp[0].first].specList[vp[0].second].decomp;
		}
		else
		{
			vp[0].first = atof(s2.c_str()) + counter;
			iss2 >> vp[0].second >> locus.i >> locus.j >> locus.k;
			locus += offset;
			vp[0].second++;
		}
		constraint * c = new constraint(vp[0], locus, ground, k, RBlist); sim.interactionList[intIdx] = c; intIdx++;
	}
	if (s == "constrainAxis")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> vp[0].first >> vp[0].second >> locus.i >> locus.j >> locus.k;
		locus.rotate(rotate);
		vp[0].second++; vp[0].first += counter;
		constrainAxis * c = new constrainAxis(vp[0], locus, ground, k, RBlist); sim.interactionList[intIdx] = c; intIdx++;
	}
	if (s == "constrainAxisTime")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d startaxis, endaxis;
		double t = 0;
		double duration, offset = 0, on = -1, off = -1;
		iss2 >> vp[0].first >> vp[0].second >> startaxis.i >> startaxis.j >> startaxis.k >> endaxis.i >> endaxis.j >> endaxis.k >> duration >> offset >> on >> off;
		startaxis.rotate(rotate); endaxis.rotate(rotate);
		vp[0].second++; vp[0].first += counter;
		//cout << vp[0].first << ' ' << vp[0].second << '\n';
		constrainAxisTime * c = new constrainAxisTime(vp[0], startaxis, endaxis, ground, k, duration, t, offset, on, off, RBlist);
		c->adjust(steps);
		sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "bond" || s == "chain")
	{
		bool bisChain = 0;
		if (s == "chain") bisChain = 1;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second;
		for (int i = 0; i < 2; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		bond * b = new bond(vp[0], vp[1], ground, k, RBlist); sim.interactionList[intIdx] = b; intIdx++;
		b->isChain = bisChain;
		bond temp = *b;
		sim.bondList.push_back(temp);
	}
	else if (s == "angle")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second;
		for (int i = 0; i < 3; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		angle * a = new angle(vp[0], vp[1], vp[2], ground, k, RBlist); a->obs = observe; sim.interactionList[intIdx] = a; intIdx++;
	}
	/*else if (s == "quartAngle")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second;
		for (int i = 0; i < 3; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		quartAngle * a = new quartAngle(vp[0], vp[1], vp[2], ground, k, RBlist); a->obs = observe; sim.interactionList[intIdx] = a; intIdx++;
	}*/
	else if (s == "gaussAngle")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second;
		for (int i = 0; i < 3; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		gaussAngle * a = new gaussAngle(vp[0], vp[1], vp[2], ground, k, RBlist); a->obs = observe; sim.interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "angleAxis")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d axis;
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> axis.i >> axis.j >> axis.k;
		for (int i = 0; i < 2; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		angleAxis * a = new angleAxis(vp[0], vp[1], axis, ground, k, RBlist); a->obs = observe; sim.interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "dihedral")
	{
		iss >> divisor;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> vp[1].first >> vp[1].second >> vp[2].first >> vp[2].second >> vp[3].first >> vp[3].second;
		for (int i = 0; i < 4; i++)
		{
			vp[i].first += counter;
			vp[i].second++;
		}
		dihedral * d = new dihedral(vp[0], vp[1], vp[2], vp[3], ground, k, divisor, RBlist); sim.interactionList[intIdx] = d; intIdx++;
	}
	else if (s == "extDipole")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second;
		vp[0].second++; vp[0].first += counter;
		extDipole * c = new extDipole(vp[0], &sim.externalDipole, &sim.prevDipole, &sim.expenditure, ground, k, RBlist, 0, sim.rotationAxis); sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "extDipoleShift")
	{
		double angle = 0;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> angle; angle *= PI/180;
		vp[0].second++; vp[0].first += counter;
		extDipole * c = new extDipole(vp[0], &sim.externalDipole, &sim.prevDipole, &sim.expenditure, ground, k, RBlist, angle, sim.rotationAxis); sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "staccatoDipole")
	{
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second;
		vp[0].second++; vp[0].first += counter;
		double start = 0, halt = 1, silenceStart = 0, silenceEnd = 0;
		iss2 >> start >> halt >> silenceStart >> silenceEnd;
		staccatoDipole * c = new staccatoDipole(vp[0], &sim.externalDipole, &sim.prevDipole, &sim.expenditure, ground, k, RBlist, 0, sim.rotationAxis, start, halt, silenceStart,
		silenceEnd, &sim.sinDipoleAngle);
		sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "offDipole")
	{
		double offset = 0;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second;
		vp[0].second++; vp[0].first += counter;
		offDipole * c = new offDipole(vp[0], &sim.offsetDipole, &sim.prevOffDipole, &sim.expenditure, ground, k, RBlist, 0, sim.rotationAxis); sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "offDipoleShift")
	{
		double offset = 0, angle = 0;
		getline(ifs, s); istringstream iss2(s);
		iss2 >> vp[0].first >> vp[0].second >> angle; angle *= PI/180;
		vp[0].second++; vp[0].first += counter;
		offDipole * c = new offDipole(vp[0], &sim.offsetDipole, &sim.prevOffDipole, &sim.expenditure, ground, k, RBlist, angle, sim.rotationAxis); sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "LJstrength") sim.LJstrength[(int)ground + counter] = k;
	sim.interactionList[intIdx - 1]->min = min;
	sim.interactionList[intIdx - 1]->equil = equil;
	sim.interactionList[intIdx - 1]->umb = umb;
	sim.interactionList[intIdx - 1]->active = active;
}

void parseSuperInteraction(ifstream &ifs, string &s, vector<int> &RBlistCounter, simulation &sim, int &intIdx, int steps)
{
	rigidBody * RBlist = sim.RBlist;
	istringstream iss(s);
	double ground, k, divisor;
	vector<triple> vp(4); vector<pair<int, int> > p(4);
	bool observe = 0, min = 0, equil = 0, umb = 0, active = 0;
	iss >> s >> ground >> k;
	if (s == "constraint")
	{
		string s2;
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> s2;
		if (s2 == "curr")
		{
			iss2 >> vp[0].a >> vp[0].b >> vp[0].c;
			p[0].first = vp[0].b + RBlistCounter[vp[0].a]; p[0].second = vp[0].c + 1;
			locus = RBlist[p[0].first].specList[p[0].second].decomp;
		}
		else
		{
			vp[0].a = atoi(s2.c_str());
			iss2 >> vp[0].b >> vp[0].c >> locus.i >> locus.j >> locus.k;
			p[0].first = vp[0].b + RBlistCounter[vp[0].a]; p[0].second = vp[0].c + 1;
		}
		constraint * c = new constraint(p[0], locus, ground, k, RBlist); sim.interactionList[intIdx] = c; intIdx++;
	}
	if (s == "constrainAxis")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d locus;
		iss2 >> vp[0].a >> vp[0].b >> vp[0].c >> locus.i >> locus.j >> locus.k;
		p[0].first = vp[0].b + RBlistCounter[vp[0].a]; p[0].second = vp[0].c + 1;
		constrainAxis * c = new constrainAxis(p[0], locus, ground, k, RBlist); sim.interactionList[intIdx] = c; intIdx++;
	}
	if (s == "constrainAxisTime")
	{
		getline(ifs, s); istringstream iss2(s);
		vector3d startaxis, endaxis;
		double t = 0;
		double duration, offset = 0, on = -1, off = -1;
		iss2 >> vp[0].a >> vp[0].b >> vp[0].c >> startaxis.i >> startaxis.j >> startaxis.k >> endaxis.i >> endaxis.j >> endaxis.k >> duration >> offset >> on >> off;
		p[0].first = vp[0].b + RBlistCounter[vp[0].a]; p[0].second = vp[0].c + 1;
		constrainAxisTime * c = new constrainAxisTime(p[0], startaxis, endaxis, ground, k, duration, t, offset, on, off, RBlist);
		c->adjust(steps);
		sim.interactionList[intIdx] = c; intIdx++;
	}
	else if (s == "bond" || s == "chain")
	{
		bool bisChain = 0;
		if (s == "chain") bisChain = 1;
		getline(ifs, s); istringstream iss2(s);
		for (int i = 0; i < 2; i++)
		{
			iss2 >> vp[i].a >> vp[i].b >> vp[i].c;
			p[i].first = vp[i].b + RBlistCounter[vp[i].a]; p[i].second = vp[i].c + 1;
			//cout << "hi: " << p[i].first << '\n';
		}
		bond * b = new bond(p[0], p[1], ground, k, RBlist); sim.interactionList[intIdx] = b; intIdx++;
		b->isChain = bisChain;
		bond temp = *b;
		sim.bondList.push_back(temp);
	}
	else if (s == "angle")
	{
		getline(ifs, s); istringstream iss2(s);
		for (int i = 0; i < 3; i++)
		{
			iss2 >> vp[i].a >> vp[i].b >> vp[i].c;
			p[i].first = vp[i].b + RBlistCounter[vp[i].a]; p[i].second = vp[i].c + 1;
		}
		angle * a = new angle(p[0], p[1], p[2], ground, k, RBlist); sim.interactionList[intIdx] = a; intIdx++;
	}
	else if (s == "dihedral")
	{
		iss >> divisor;
		getline(ifs, s); istringstream iss2(s);
		for (int i = 0; i < 4; i++)
		{
			iss2 >> vp[i].a >> vp[i].b >> vp[i].c;
			p[i].first = vp[i].b + RBlistCounter[vp[i].a]; p[i].second = vp[i].c + 1;
		}
		dihedral * d = new dihedral(p[0], p[1], p[2], p[3], ground, k, divisor, RBlist); sim.interactionList[intIdx] = d; intIdx++;
	}
	sim.interactionList[intIdx - 1]->min = min;
	sim.interactionList[intIdx - 1]->equil = equil;
	sim.interactionList[intIdx - 1]->umb = umb;
	sim.interactionList[intIdx - 1]->active = active;
}

void probeInteraction(ifstream &ifs, string &s, int &intIdx)
{
	istringstream iss(s);
	string aux;
	iss >> s >> aux;
	if (s == "constraint")
	{
		intIdx++;
	}
	else if (s == "constrainAxis" || s == "extDipole" || s == "offDipole" || s == "extDipoleShift" || s == "offDipoleShift" || s == "staccatoDipole" || s == "constrainAxisTime")
	{
		intIdx++;
	}
	else if (s == "bond" || s == "chain")
	{
		intIdx++;
	}
	else if (s == "angle" || s == "angleAxis" || s == "gaussAngle")
	{
		intIdx++;
	}
	else if (s == "dihedral")
	{
		intIdx++;
	}
	getline(ifs, s);
}

void parseFrag(ifstream &ifs, vector3d &offset, int counter, quaternion &orient, simulation &sim, bool ind, triple &count, int steps)
{
	int i, j;
	string s;
	if (ind)
	{
		getline(ifs, s);
		sim.RBlistSize = 0;
		while (true) //parse list of rigid bodies
		{
			if (s[0] != ';')
			{
				//cout << s[0] << '\n';
				if (s[0] == '!') break;
				//cout << s << '\n';
				if (s[0] == 'r' || s[0] == 'p' || s[0] == 'd') sim.RBlistSize++;
			}
			getline(ifs, s);
		}
		getline(ifs, s);
		sim.intListSize = 0;
		if (s[0] != '!')
		{
			while (true) //parse list of interactions
			{
				if (s[0] != ';')
				{
					if (s[0] == '!') break;
					probeInteraction(ifs, s, sim.intListSize);
				}
				getline(ifs, s);
			}
		}
		ifs.clear();
		ifs.seekg(0, ios::beg);
		sim.RBlist = new rigidBody[sim.RBlistSize];
		sim.interactionList = new interaction*[sim.intListSize];
	}
	//cout << sim.RBlistSize << ' ' << sim.intListSize << "\n";
	rigidBody * RBlist = sim.RBlist; int startCount = count.a;
	getline(ifs, s);
	while (true) //parse list of rigid bodies
	{
		if (s[0] != ';')
		{
			//cout << s[0] << '\n';
			if (s[0] == '!') break;
			if (s[0] == 'r')
			{
				parseRigid(ifs, s, offset, counter, orient, sim, count.a);
				count.a++;
				//cout << count.a << '\n';
			}
			else if (s[0] == 'd')
			{
				parseRigid(ifs, s, offset, counter, orient, sim, count.a);
				RBlist[count.a].isDummy = 1;
				count.a++;
				//cout << count.a << '\n';
			}
			else if (s[0] == 'p')
			{
				parsePoint(ifs, s, offset, sim, count.a);
				getline(ifs, s);
				count.a++;
			}
		}
		else getline(ifs, s);
	}
	getline(ifs, s);
	//sim.decompAll();
	if (s[0] != '!')
	{
		while (true) //parse list of interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				parseInteraction(ifs, s, offset, counter, orient, sim, count.b, steps);
				//cout << count.b << '\n';
			}
			getline(ifs, s);
		}
	}
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true)
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				istringstream iss(s);
				iss >> s;
				if (s == "suppressLJ")
				{
					iss >> i >> j;
					i += counter; j += counter;
					//cout << i << ' ' << j << '\n';
					sim.suppressLJ.insert(sim.suppressf(i, j));
					sim.suppressLJ.insert(sim.suppressf(j, i));
					pair<int, int> p(i, j);
					sim.suppressed.push_back(p);
				}
				else if (s == "umbrellaLJ")
				{
					iss >> i >> j;
					i += counter; j += counter;
					//cout << i << ' ' << j << '\n';
					sim.umbrellaLJ.insert(sim.suppressf(i, j));
					sim.umbrellaLJ.insert(sim.suppressf(j, i));
					pair<int, int> p(i, j);
					sim.umbrellaed.push_back(p);
					//sim.suppressed.push_back(p);
				}
			}
			getline(ifs, s);
		}
	}
	for (int i = startCount; i < count.a; i++)
	{
		//cout << &RBlist[i] << '\n';
		//cout << "ID: " << i << '\n';
		RBlist[i].recenter();
		RBlist[i].MST();
		//RBlist[i].printAdjList();
	}
	if (ind)
	{
		getline(ifs, s);
		if (s[0] != '!')
		{
			while (true)
			{
				if (s[0] != ';')
				{
					if (s[0] == '!') break;
					istringstream iss(s);
					iss >> s;
					if (s == "spec")
					{
						iss >> i >> j;
						pair<int, int> p(i, j);
						sim.topSpecList.push_back(p);
					}
				}
				getline(ifs, s);
			}
		}
		sim.populateIntTable();
		sim.initDecompAll();
	}
}

void probeFrag(ifstream &ifs, int &RBlistSize, int &intListSize)
{
	string s;
	getline(ifs, s);
	while (true) //parse list of rigid bodies
	{
		if (s[0] != ';')
		{
			//cout << s[0] << '\n';
			if (s[0] == '!') break;
			//cout << s << '\n';
			if (s[0] == 'r' || s[0] == 'p' || s[0] == 'd') RBlistSize++;
		}
		getline(ifs, s);
	}
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true) //parse list of interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				probeInteraction(ifs, s, intListSize);
			}
			getline(ifs, s);
		}
	}
}

void parseSuper(ifstream &ifs, simulation &sim, int steps)
{
	vector<int> RBlistCounter; RBlistCounter.push_back(0);
	double a, b, c; int i = 0, j, k, l; sim.RBlistSize = 0; sim.intListSize = 0;
	string s;
	getline(ifs, s);
	while (true) //parse list of input files
	{
		//cout << "obj: " << i << '\n';
		if (s[0] != ';')
		{
			if (s[0] == '!') break;
			istringstream iss(s); iss >> s;
			ifstream ifsf(s.c_str());
			probeFrag(ifsf, sim.RBlistSize, sim.intListSize);
		}
		getline(ifs, s);
	}
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true) //probe super interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				probeInteraction(ifs, s, sim.intListSize);
			}
			getline(ifs, s);
		}
	}
	ifs.clear();
	ifs.seekg(0, ios::beg);
	sim.RBlist = new rigidBody[sim.RBlistSize];
	sim.interactionList = new interaction*[sim.intListSize];
	//cout << sim.RBlistSize << ' ' << sim.intListSize << '\n';
	triple indices(0, 0, 0);
	getline(ifs, s);
	while (true) //parse list of input files
	{
		//cout << "obj: " << i << '\n';
		if (s[0] != ';')
		{
			if (s[0] == '!') break;
			istringstream iss(s); iss >> s >> a >> b >> c;
			vector3d offset(a, b, c);
			iss >> a >> b >> c;
			vector3d euler(a, b, c);
			quaternion orient = eulerConv(euler);
			cout << "Parsing " << s << '\n';
			ifstream ifsf(s.c_str());
			parseFrag(ifsf, offset, RBlistCounter[i], orient, sim, 0, indices, steps);
			cout << "Parsed " << s << '\n';
			RBlistCounter.push_back(indices.a);
			i++;
		}
		getline(ifs, s);
	}
	sim.populateIntTable();
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true) //stitch together disparate structures with new interactions
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				parseSuperInteraction(ifs, s, RBlistCounter, sim, indices.b, steps);
			}
			getline(ifs, s);
		}
	}
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true)
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				istringstream iss(s);
				iss >> s;
				if (s == "suppressLJ")
				{
					iss >> i >> j >> k >> l;
					i = RBlistCounter[i] + j; j = RBlistCounter[k] + l;
					//cout << i << ' ' << j << '\n';
					sim.suppressLJ.insert(sim.suppressf(i, j));
					sim.suppressLJ.insert(sim.suppressf(j, i));
					pair<int, int> p(i, j);
					sim.suppressed.push_back(p);
				}
			}
			getline(ifs, s);
		}
	}
	getline(ifs, s);
	if (s[0] != '!')
	{
		while (true)
		{
			if (s[0] != ';')
			{
				if (s[0] == '!') break;
				istringstream iss(s);
				iss >> s;
				if (s == "spec")
				{
					iss >> i >> j >> k;
					j = RBlistCounter[i] + j;
					pair<int, int> p(j, k);
					sim.topSpecList.push_back(p);
				}
			}
			getline(ifs, s);
		}
	}
	sim.initDecompAll();
}
