#include "RadiationFlatUniform.h"

template<class Coord>
RadiationFlatUniform<Coord>::RadiationFlatUniform(SimulationData<Coord>& simData_):
simData(simData_), grid(simData_.metric.grid), nDir(simData_.nDir), stencil(simData_.uniStencil)
{
	E        = new double[grid.n12]();
	Fx       = new double[grid.n12]();
	Fy       = new double[grid.n12]();
	Pxx      = new double[grid.n12]();
	Pxy      = new double[grid.n12]();
	Pyy      = new double[grid.n12]();
	Fnorm	 = new double[grid.n12]();
	kappa0   = new double[grid.n12]();
	kappa1   = new double[grid.n12]();
	kappaA   = new double[grid.n12]();
	eta      = new double[grid.n12]();
	I        = new double[grid.n12 * nDir]();
	temp     = new double[grid.n12 * nDir]();

	LoadInitialData();
	ComputeMoments();
	NormalizeInitialData();
	LoadInitialData(); // loads normalized data into rotated stencil.
}



template<class Coord>
RadiationFlatUniform<Coord>::~RadiationFlatUniform()
{
	delete[] E;
	delete[] Fx;
	delete[] Fy;
	delete[] Pxx;
	delete[] Pxy;
	delete[] Pyy;
	delete[] kappa0;
	delete[] kappa1;
	delete[] kappaA;
	delete[] eta;
	delete[] I;
	delete[] temp;
}



template<class Coord>
void RadiationFlatUniform<Coord>::LoadInitialData()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		kappa0[ij] = simData.initialKappa0[ij];
		eta[ij] = simData.initialEta[ij];
		double intensity = simData.initialE[ij];
		if(intensity > 0)
		{
			for(int d=0; d<nDir; d++)
				I[ij + d*grid.n12] = intensity;
		}
	}
}



template<class Coord>
void RadiationFlatUniform<Coord>::NormalizeInitialData()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		if (E[ij] > 0)
			simData.initialE[ij] *= simData.initialE[ij] / E[ij];
}



template<class Coord>
void RadiationFlatUniform<Coord>::InitVakuumMicrophysics()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{	
		kappa0[ij] = 0;
		kappa1[ij] = 0;
		kappaA[ij] = 0;
		eta   [ij] = 0;
	}
}



template<class Coord>
void RadiationFlatUniform<Coord>::ComputeMoments()
{
	constexpr double twoPiInv = 1.0 / (2.0 * M_PI);

	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		E [ij]  = 0.0;
		Fx[ij]  = 0.0;
		Fy[ij]  = 0.0;
		Pxx[ij] = 0.0;
		Pxy[ij] = 0.0;
		Pyy[ij] = 0.0;
		for(int d = 0; d < nDir; d++)
		{
			E [ij]  += stencil.W(d) * I[ij + d*grid.n12];
			Fx[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d);
			Fy[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d);
			Pxx[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d) * stencil.Cx(d);
			Pxy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d) * stencil.Cy(d);
			Pyy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d) * stencil.Cy(d);
		}
        E [ij]  *= twoPiInv;
        Fx[ij]  *= twoPiInv;
        Fy[ij]  *= twoPiInv;
        Pxx[ij] *= twoPiInv;
        Pxy[ij] *= twoPiInv;
        Pyy[ij] *= twoPiInv;
		Fnorm[ij] = Tensor2<xy,IF>(Fx[ij],Fy[ij]).EuklNorm();
	}
}



template<class Coord>
void RadiationFlatUniform<Coord>::Stream()
{
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);
		int ijd = grid.Index(i,j,d);
		
		// Flat Streaming:
		Tensor2<xy,IF> vTemp = stencil.Cxy(d);
		Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
		xTemp[1] -= vTemp[1] * grid.dt;
		xTemp[2] -= vTemp[2] * grid.dt;

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= BilinearInterpolation(iTemp-i0,jTemp-j0,
		 I[grid.Index(i0,j0) + d*grid.n12],I[grid.Index(i0,j1) + d*grid.n12],
		 I[grid.Index(i1,j0) + d*grid.n12],I[grid.Index(i1,j1) + d*grid.n12]);
	}
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}



template<class Coord>
void RadiationFlatUniform<Coord>::Collide()
{
	#pragma omp parallel for
	for(int j = 0; j < grid.n2; j++)
	for(int i = 0; i < grid.n1; i++)
	{
		int ij = grid.Index(i,j);
		if(kappa0[ij] > 0)
		{
			for(int d=0; d<nDir; d++)
				I[grid.Index(i,j,d)] = 0;
		}
		if(eta[ij] > 0)
		{
			for(int d=0; d<nDir; d++)
				I[grid.Index(i,j,d)] = simData.initialE[i,j];
		}
	}
//	// TODO:
//	// Steife DGL?
//	
//	#pragma omp parallel for
//	for(int j = 0; j < metric.grid.n2; j++)
//	for(int i = 0; i < metric.grid.n1; i++)
//	{
//		Tensor3 x = metric.grid.xyzCoord(i,j,1);
//		if(metric.ToCloseToBH(x))
//			continue;
//
//		int ij = i+j*metric.grid.n1;
//
//		// Select correct stencil:
//		Stencil& stencil = rotation[ij] < 0 ? uniStencil : dirStencil;
//		
//		double Gamma = eta[ij];
//		for(int d = 0; d < velSt.nDir; d++)
//		{
//			int ijd = ij + d*metric.grid.n12;
//			Gamma += -kappaA[ij] * I[ijd] + kappa0[ij]*(E[ij]-I[ijd]) + 3*(velSt.vx[d]*Fx[ij] + velSt.vy[d]*Fy[ij]);
//			I[ijd] += metric.GetAlpha(x) * Gamma * metric.grid.dt;
//		}
//	}
//	
}



template<class Coord>
void RadiationFlatUniform<Coord>::RunSimulation(int writeFrequency, bool updateFourierHarmonics, bool keepSourceNodesActive, bool writeData, bool printToTerminal)
{
	// -------------------- Initialization --------------------
	Timer<10>::Reset("MainLoop");
	Timer<20>::Reset("ComputeMoments");
	Timer<30>::Reset("WriteDataToJson");
	Timer<40>::Reset("Collide");
	Timer<60>::Reset("Stream");
	Timer<70>::Reset("KeepSourceNodesActive");

	// Initial data output:
	if (printToTerminal)
	{
		// std::cout << " sigma=" << simData.sigma << std::endl;
		// std::cout << " n1=" << grid.n1 << std::endl;
		// std::cout << " n2=" << grid.n2 << std::endl;
		// std::cout << " nDir=" << nDir << std::endl;
		// std::cout << " simTime=" << simData.simTime << std::endl;
		// std::cout << " dt=" << grid.dt << std::endl;
		// std::cout << " timeSteps=" << simData.timeSteps << std::endl;
		// std::cout << " filesToWrite=" << simData.timeSteps/writeFrequency << std::endl;
	}
	// --------------------------------------------------------



	// ----------------- Main simulation Loop -----------------
	Timer<10>::Start();
	for(int n=0; n<simData.timeSteps; n++)
	{
		std::cout << n << "," << std::flush;

		Timer<20>::Start();
		ComputeMoments();
		Timer<20>::End();

		Timer<40>::Start();
		Collide();
		Timer<40>::End();

		if (writeData && (n % writeFrequency) == 0)
		{
			Timer<30>::Start();
			grid.WriteFrametoJson(n*grid.dt,E,Fx,Fy,Fnorm,n,simData.directoryPath + "/Moments");
			Timer<30>::End();
		}

		Timer<60>::Start();
		Stream();
		Timer<60>::End();

		if (keepSourceNodesActive)
		{
			Timer<70>::Start();
			LoadInitialData();
			Timer<70>::End();
		}
	}
	std::cout << std::endl;
	if (writeData)
	{
		Timer<30>::Start();
		grid.WriteFrametoJson(simData.timeSteps*grid.dt,E,Fx,Fy,Fnorm,simData.timeSteps,simData.directoryPath + "/Moments");
		Timer<30>::End();
	}
	Timer<10>::End();
	
	if (printToTerminal)
	{
		Timer<10>::Print();
		Timer<20>::Print();
		Timer<30>::Print();
		Timer<40>::Print();
		Timer<60>::Print();
		Timer<70>::Print();
	}
	// --------------------------------------------------------



	// ---------------------- Termination ---------------------
	simData.AddTimeMeasurement(Timer<10>::name, Timer<10>::time);
	simData.AddTimeMeasurement(Timer<20>::name, Timer<20>::time);
	simData.AddTimeMeasurement(Timer<30>::name, Timer<30>::time);
	simData.AddTimeMeasurement(Timer<40>::name, Timer<40>::time);
	simData.AddTimeMeasurement(Timer<60>::name, Timer<60>::time);
	simData.AddTimeMeasurement(Timer<70>::name, Timer<70>::time);
	simData.LogSimulationParameters();
	// --------------------------------------------------------
}



template class RadiationFlatUniform<xy>;
template class RadiationFlatUniform<rph>;