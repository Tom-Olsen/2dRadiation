#include "RadiationFlat.h"

template<class Coord>
RadiationFlat<Coord>::RadiationFlat(SimulationData<Coord>& simData_):
simData(simData_), grid(simData_.metric.grid), nDir(simData_.nDir),
uniStencil(simData_.uniStencil), dirStencil(simData_.dirStencil), momentStencil(simData_.momentStencil)
{
	rotation = new double[grid.n12]();
	E        = new double[grid.n12]();
	Fx       = new double[grid.n12]();
	Fy       = new double[grid.n12]();
	Pxx      = new double[grid.n12]();
	Pxy      = new double[grid.n12]();
	Pyy      = new double[grid.n12]();
	Fnorm    = new double[grid.n12]();
	kappa0   = new double[grid.n12]();
	kappa1   = new double[grid.n12]();
	kappaA   = new double[grid.n12]();
	eta      = new double[grid.n12]();
	I        = new double[grid.n12 * nDir]();
	temp     = new double[grid.n12 * nDir]();

	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		rotation[ij] = -1;

	LoadInitialData();
	ComputeMoments();
	NormalizeInitialData();
	LoadInitialData(); // loads normalized data into rotated stencil.
}


template<class Coord>
RadiationFlat<Coord>::~RadiationFlat()
{
	delete[] rotation;
	delete[] E;
	delete[] Fx;
	delete[] Fy;
	delete[] Pxx;
	delete[] Pxy;
	delete[] Pyy;
	delete[] Fnorm;
	delete[] kappa0;
	delete[] kappa1;
	delete[] kappaA;
	delete[] eta;
	delete[] I;
	delete[] temp;
}



template<class Coord>
void RadiationFlat<Coord>::LoadInitialData()
{
	auto IntensityDistribution = [](double phi, double sigma)
	{
		return exp(-0.5 * pow((phi) / sigma, 2))
		     + exp(-0.5 * pow((phi - 2.0*M_PI) / sigma, 2))
		     + exp(-0.5 * pow((phi + 2.0*M_PI) / sigma, 2));
	};

	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		kappa0[ij] = simData.initialKappa0[ij];
		eta[ij] = simData.initialEta[ij];
		double intensity = simData.initialE[ij];
		if(intensity > 0)
		{
			rotation[ij] = simData.initialAngle[ij];
			if(rotation[ij] > 0)
			{
				Stencil& stencil = dirStencil;
				// Set intensities such that peak is at I[0] and I[nDir-1].
				// This is done by NOT rotating the stencil yet and leaving it to be pointing right.
				for(int d=0; d<nDir; d++)
					I[ij + d*grid.n12] = intensity * IntensityDistribution(stencil.Phi(d), simData.sigma);
			}
			else
			{
				for(int d=0; d<nDir; d++)
					I[ij + d*grid.n12] = intensity;
			}
		}
	}
}



template<class Coord>
void RadiationFlat<Coord>::NormalizeInitialData()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		if (E[ij] > 0)
			simData.initialE[ij] *= simData.initialE[ij] / E[ij];
}



template<class Coord>
void RadiationFlat<Coord>::InitVakuumMicrophysics()
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
void RadiationFlat<Coord>::ComputeMoments()
{
	constexpr double twoPiInv = 1.0 / (2.0 * M_PI);
	
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		Stencil& stencil = rotation[ij] < 0 ? uniStencil : dirStencil;
		
		E [ij]  = 0.0;
		Fx[ij]  = 0.0;
		Fy[ij]  = 0.0;
		Pxx[ij] = 0.0;
		Pxy[ij] = 0.0;
		Pyy[ij] = 0.0;
		for(int d=0; d<momentStencil.nDir; d++)
		{
			double phi = momentStencil.Phi(d);			// angle of current direction in moment stencil
			double k = stencil.Index(phi,rotation[ij]);	// k value of above angle in original stencil
	
			int kFloor = floor(k);
			int kM1 = (kFloor - 1 + nDir) % nDir;
			int kP0 = (kFloor + 0 + nDir) % nDir;
			int kP1 = (kFloor + 1 + nDir) % nDir;
			int kP2 = (kFloor + 2 + nDir) % nDir;
	
			double Id = CubicInterpolation (k - kFloor,I[ij + kM1*grid.n12],I[ij + kP0*grid.n12],I[ij + kP1*grid.n12],I[ij + kP2*grid.n12]);
			// double Id = LinearInterpolation(k - kFloor,I[ij + kP0*grid.n12],I[ij + kP1*grid.n12]);

			E [ij]  += momentStencil.W(d) * Id;
			Fx[ij]  += momentStencil.W(d) * Id * momentStencil.Cx(d);
			Fy[ij]  += momentStencil.W(d) * Id * momentStencil.Cy(d);
			Pxx[ij] += momentStencil.W(d) * Id * momentStencil.Cx(d) * momentStencil.Cx(d);
			Pxy[ij] += momentStencil.W(d) * Id * momentStencil.Cx(d) * momentStencil.Cy(d);
			Pyy[ij] += momentStencil.W(d) * Id * momentStencil.Cy(d) * momentStencil.Cy(d);
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
void RadiationFlat<Coord>::RotateStencil()
{
	#pragma omp parallel for
	for(int i=1; i<grid.n1-1; i++)
	for(int j=1; j<grid.n2-1; j++)
	{
		int ij = grid.Index(i,j);
		double rotationOld = rotation[ij];

		// Find new rotation:
		{
			// Average F on 9 point cluster:
			Tensor2<Coord,IF> averageF(0.0);
			for(int a=-1; a<=1; a++)
				for(int b=-1; b<=1; b++)
				{
					int index = grid.Index(i+a,j+b);
					averageF[1] += Fx[index];
					averageF[2] += Fy[index];
				}
			averageF[1] /= 9.0;
			averageF[2] /= 9.0;

			constexpr double rotationThreshold = 1e-3;
			if(averageF.EuklNorm() > rotationThreshold)
				rotation[ij] = averageF.Angle();
			else
				rotation[ij] = -1;
		}
		double rotationNew = rotation[ij];
		
		// Set stencils
		Stencil& stencilOld = rotationOld < 0 ? uniStencil : dirStencil;
		Stencil& stencilNew = rotationNew < 0 ? uniStencil : dirStencil;

		// Interpolate intensities to new stencil:
		if(rotationNew != rotationOld)
		{
			double Irotated[nDir];
			for(int d=0; d<nDir; d++)
			{
				double phi = stencilNew.Phi(d,rotationNew);		// angle of current direction in new stencil
				double k = stencilOld.Index(phi,rotationOld);	// k value of above angle in old stencil

				int kFloor = floor(k);
				int kM1 = (kFloor - 1 + nDir) % nDir;
				int kP0 = (kFloor + 0 + nDir) % nDir;
				int kP1 = (kFloor + 1 + nDir) % nDir;
				int kP2 = (kFloor + 2 + nDir) % nDir;

				Irotated[d] = CubicInterpolation (k - kFloor,I[ij + kM1*grid.n12],I[ij + kP0*grid.n12],I[ij + kP1*grid.n12],I[ij + kP2*grid.n12]);
				// Irotated[d] = LinearInterpolation(k - kFloor,I[ij + kP0*grid.n12],I[ij + kP1*grid.n12]);
			}
			for(int d=0; d<nDir; d++)
				I[ij + d * grid.n12] = Irotated[d];
		}
	}
}



template<class Coord>
inline __attribute__((always_inline)) double RadiationFlat<Coord>::InterpolateIntensity(int ij, double phi)
{
	double k = rotation[ij] < 0 ? uniStencil.Index(phi,rotation[ij]) : dirStencil.Index(phi,rotation[ij]);

	int kFloor = floor(k);
	int kM1 = (kFloor - 1 + nDir) % nDir;
	int kP0 = (kFloor + 0 + nDir) % nDir;
	int kP1 = (kFloor + 1 + nDir) % nDir;
	int kP2 = (kFloor + 2 + nDir) % nDir;

	return CubicInterpolation (k - kFloor,I[ij + kM1*grid.n12],I[ij + kP0*grid.n12],I[ij + kP1*grid.n12],I[ij + kP2*grid.n12]);
	// return LinearInterpolation(k - kFloor,I[ij + kP0*grid.n12],I[ij + kP1*grid.n12]);
}



template<class Coord>
void RadiationFlat<Coord>::Stream()
{
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Select correct stencil:
		Stencil& stencil = rotation[ij] < 0 ? uniStencil : dirStencil;

		// Flat Streaming:
		Tensor2<xy,IF> vTemp = stencil.Cxy(d,rotation[ij]);
		Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
		xTemp[1] -= vTemp[1] * grid.dt;
		xTemp[2] -= vTemp[2] * grid.dt;

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Intensity interpolation:
		double angle = stencil.Phi(d,rotation[ij]);
		double intensityAt_i0j0 = InterpolateIntensity(grid.Index(i0,j0),angle);
		double intensityAt_i0j1 = InterpolateIntensity(grid.Index(i0,j1),angle);
		double intensityAt_i1j0 = InterpolateIntensity(grid.Index(i1,j0),angle);
		double intensityAt_i1j1 = InterpolateIntensity(grid.Index(i1,j1),angle);

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= BilinearInterpolation(iTemp-i0,jTemp-j0,
		 intensityAt_i0j0,intensityAt_i0j1,
		 intensityAt_i1j0,intensityAt_i1j1);
	}
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}



template<class Coord>
void RadiationFlat<Coord>::Collide()
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
void RadiationFlat<Coord>::RunSimulation(int writeFrequency, bool updateFourierHarmonics, bool keepSourceNodesActive, bool writeData, bool printToTerminal)
{
	// -------------------- Initialization --------------------
	Timer<10>::Reset("MainLoop");
	Timer<20>::Reset("ComputeMoments");
	Timer<30>::Reset("WriteDataToJson");
	Timer<40>::Reset("Collide");
	Timer<60>::Reset("Stream");
	Timer<70>::Reset("KeepSourceNodesActive");
	Timer<90>::Reset("Rotate Stencil");

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

		Timer<90>::Start();
		RotateStencil();
		Timer<90>::End();

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
		Timer<90>::Print();
		std::cout << std::endl;
	}
	// --------------------------------------------------------



	// ---------------------- Termination ---------------------
	simData.AddTimeMeasurement(Timer<10>::name, Timer<10>::time);
	simData.AddTimeMeasurement(Timer<20>::name, Timer<20>::time);
	simData.AddTimeMeasurement(Timer<30>::name, Timer<30>::time);
	simData.AddTimeMeasurement(Timer<40>::name, Timer<40>::time);
	simData.AddTimeMeasurement(Timer<60>::name, Timer<60>::time);
	simData.AddTimeMeasurement(Timer<70>::name, Timer<70>::time);
	simData.AddTimeMeasurement(Timer<90>::name, Timer<90>::time);
	simData.LogSimulationParameters();
	// --------------------------------------------------------
}



template class RadiationFlat<xy>;
template class RadiationFlat<rph>;