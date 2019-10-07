//----------Model Path----------
	/*std::string model_path = "Samples/Advection/";

	//----------Add Nodes----------
	std::vector<Vector<double> > nodes;
	ImportNodesFromCSV(nodes, model_path + "Node.csv");

	//----------Add Elements----------
	std::vector<std::vector<int> > elements;
	ImportElementsFromCSV(elements, model_path + "Element.csv");

	//----------Add Field----------
	std::vector<int> field;
	int KDEGREE = 0;
	ImportFieldFromCSV(field, KDEGREE, nodes.size(), model_path + "Field.csv");

	//----------Add Dirichlet Condition----------
	std::vector<int> isufixed;
	std::vector<double> ufixed;
	ImportDirichletFromCSV(isufixed, ufixed, field, model_path + "Dirichlet.csv");

	//----------Add Initial Condition----------
	std::vector<double> u = std::vector<double>(nodes.size(), 0.0);
	u[12] = 100.0;

	//----------Time Step----------
	double dt = 0.01;
	double theta = 0.5;
	for (int k = 0; k < 1000; k++) {
		//----------Culculate Ke Ce and Assembling----------
		LILCSR<double> K = LILCSR<double>(KDEGREE, KDEGREE);
		LILCSR<double> C = LILCSR<double>(KDEGREE, KDEGREE);
		for (auto element : elements) {
			std::vector<std::vector<double> > Ke = AdvectionTri(nodes, element, 1.0, 1.0, 1.0);
			Assembling(K, Ke, element, field);
			std::vector<std::vector<double> > Ce = MassTri(nodes, element, 1.0);
			Assembling(C, Ce, element, field);
		}

		//----------Make Equation----------
		LILCSR<double> KC = (1.0 / dt)*C + theta * K;
		std::vector<double> F = ((1.0 / dt)*C - (1.0 - theta) * K) * u;
		std::cout << KC << std::endl;

		//----------Set Dirichlet Boundary Condition----------
		SetDirichlet(KC, F, isufixed, ufixed, 1.0e3);

		//----------Solve System Equation----------
		CSR<double> Kmod = CSR<double>(KC);
		CSR<double> M = ILU0(Kmod);
		std::vector<double> result = ILU0BiCGSTAB(Kmod, M, F, 10000, 1.0e-10);

		//----------Post Process----------
		u.clear();
		FieldResultToNodeValue(result, u, field);

		//----------Save file----------
		std::ofstream fout(model_path + "result" + std::to_string(k) + ".vtk");
		MakeHeadderToVTK(fout);
		AddPointsToVTK(nodes, fout);
		AddElementToVTK(elements, fout);
		std::vector<int> et = std::vector<int>(32, 5);
		AddElementTypes(et, fout);
		AddPointScalers(u, "u", fout);
		fout.close();
	}*/