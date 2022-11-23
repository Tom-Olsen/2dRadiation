#include "RadiationUniform.h"

template<class Coord>
RadiationUniform<Coord>::RadiationUniform(SimulationData<Coord>& simData_):
simData(simData_), grid(simData_.metric.grid), metric(simData_.metric), nDir(simData_.nDir),
stencil(simData_.uniStencil), fourierStencil(simData_.fourierStencil)
{
	E        = new double[grid.n12]();
	Fx       = new double[grid.n12]();
	Fy       = new double[grid.n12]();
	Pxx      = new double[grid.n12]();
	Pxy      = new double[grid.n12]();
	Pyy      = new double[grid.n12]();
	E_LF     = new double[grid.n12]();
	Fx_LF    = new double[grid.n12]();
	Fy_LF    = new double[grid.n12]();
	Pxx_LF   = new double[grid.n12]();
	Pxy_LF   = new double[grid.n12]();
	Pyy_LF   = new double[grid.n12]();
	Fnorm_LF = new double[grid.n12]();
	kappa0   = new double[grid.n12]();
	kappa1   = new double[grid.n12]();
	kappaA   = new double[grid.n12]();
	eta      = new double[grid.n12]();
	I        = new double[grid.n12 * nDir]();
	temp     = new double[grid.n12 * nDir]();
	coefficientsS   = new double[grid.n12 * fourierStencil.nDir];
	coefficientsX1  = new double[grid.n12 * fourierStencil.nDir];
	coefficientsX2  = new double[grid.n12 * fourierStencil.nDir];
	coefficientsCx1 = new double[grid.n12 * fourierStencil.nDir];
	coefficientsCx2 = new double[grid.n12 * fourierStencil.nDir];

	LoadInitialData();
	ComputeMomentsIF();
	ComputeMomentsLF();
	NormalizeInitialData();
	LoadInitialData(); // loads normalized data into rotated stencil.
	UpdateFourierCoefficients();
}



template<class Coord>
RadiationUniform<Coord>::~RadiationUniform()
{
	delete[] E;
	delete[] Fx;
	delete[] Fy;
	delete[] Pxx;
	delete[] Pxy;
	delete[] Pyy;
	delete[] E_LF;
	delete[] Fx_LF;
	delete[] Fy_LF;
	delete[] Pxx_LF;
	delete[] Pxy_LF;
	delete[] Pyy_LF;
	delete[] Fnorm_LF;
	delete[] kappa0;
	delete[] kappa1;
	delete[] kappaA;
	delete[] eta;
	delete[] I;
	delete[] temp;
	delete[] coefficientsS;
	delete[] coefficientsX1;
	delete[] coefficientsX2;
	delete[] coefficientsCx1;
	delete[] coefficientsCx2;
}



template<class Coord>
void RadiationUniform<Coord>::LoadInitialData()
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
			double rot = simData.initialAngle[ij];
			if(rot > 0)
			{
				// Set intensities such that peak is at I[0] and I[nDir-1].
				// This is done by NOT rotating the stencil yet and leaving it to be pointing right.
				for(int d=0; d<nDir; d++)
					I[ij + d*grid.n12] = intensity * IntensityDistribution(fmod(stencil.Phi(d) + rot,2.0*M_PI), simData.sigma);
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
void RadiationUniform<Coord>::NormalizeInitialData()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		if (E_LF[ij] > 0)
			simData.initialE[ij] *= simData.initialE[ij] / E_LF[ij];
}



template<class Coord>
void RadiationUniform<Coord>::UpdateFourierCoefficients()
{
	#pragma omp parallel for
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);
		double dataS  [fourierStencil.nDir];
		double dataX1 [fourierStencil.nDir];
		double dataX2 [fourierStencil.nDir];
    	double dataCx1[fourierStencil.nDir];
    	double dataCx2[fourierStencil.nDir];
		Coordinate2<Coord> xStart = grid.x12Coord(i,j);
		Coordinate2<xy> xyCoord = grid.xyCoord(i,j);
		double alpha = metric.GetAlpha(ij);
		for(int d=0; d<fourierStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coordinate2<Coord> x = xStart;
            Tensor2<xy,IF> cxy = fourierStencil.Cxy(d);
            Tensor2<Coord,IF> c = cxy.template Transform<Coord>(xyCoord);
            Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
            Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);

			// Solve geodesic equation backwards:
			if(!metric.InsideBH(x))
			{
    	    	// s *= Euler_GeodesicEquation<-1>(grid.dt,x,v,metric);
    	    	s *= RK45_GeodesicEquation<-1>(grid.dt,x,v,metric);
			}
			else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
				v = Tensor2<Coord,LF>(0.0);

			// Final data points for fourier expansion:
			dataS[d]   = 1.0/s;
			dataX1[d]  = x[1];
			dataX2[d]  = x[2];
			dataCx1[d] = v[1];
			dataCx2[d] = v[2];
    	}
		std::vector<double> cS  = Fourier::Expansion::GetCoefficients(fourierStencil,dataS);
		std::vector<double> cX1 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX1);
		std::vector<double> cX2 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX2);
		std::vector<double> cCx1= Fourier::Expansion::GetCoefficients(fourierStencil,dataCx1);
		std::vector<double> cCx2= Fourier::Expansion::GetCoefficients(fourierStencil,dataCx2);
		for(int f=0; f<fourierStencil.nDir; f++)
		{
			int fij = f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1;
			coefficientsS  [fij] = cS  [f];
			coefficientsX1 [fij] = cX1 [f];
			coefficientsX2 [fij] = cX2 [f];
			coefficientsCx1[fij] = cCx1[f];
			coefficientsCx2[fij] = cCx2[f];
		}
	}
	// Save data to file for debugging:
	// SaveStreamingFourierExtrapolation();
}
template<class Coord>
void RadiationUniform<Coord>::SaveStreamingFourierExtrapolation()
{
	std::ofstream fileX("output/FourierHarmonicsExpansionCoord.csv");
	std::ofstream fileV("output/FourierHarmonicsExpansionVeloc.csv");
	fileX << "x,y,z,s" << std::endl;
	fileV << "x,y,z,s" << std::endl;
	for(int j=0; j<grid.n2; j++)
	for(int i=0; i<grid.n1; i++)
	{
		int ij = grid.Index(i,j);
		Coordinate2<xy> x = grid.xyCoord(i,j);
		fileX << x[1] << "," << x[2] << "," << 0 << "," << 1 << std::endl;
		fileV << x[1] << "," << x[2] << "," << 0 << "," << 1 << std::endl;
		for(int d=0; d<nDir; d++)
		{
			double phi = stencil.Phi(d);
			double s = GetFrequencyShift(i,j,phi);
			Coordinate2<Coord> xTemp = GetTempCoordinate(i,j,phi);
			Coordinate2<xy> xTempxy = xTemp.template Transform<xy>();
			Tensor2<Coord,LF> vTemp = GetTemp2Velocity(i,j,phi);
			Tensor2<xy,LF> vTempxy = vTemp.template Transform<xy>(xTemp);
			fileX << xTemp[1] << "," << xTemp[2] << "," << 0 << "," << s << std::endl;
			fileV << x[1] + vTemp[1]/20 << "," << x[2] + vTemp[2]/20 << "," << 0 << "," << s << std::endl;
		}
	}
	fileX.close();
	fileV.close();
	exit_on_error("Test data has been saved.");
}




template<class Coord>
void RadiationUniform<Coord>::InitVakuumMicrophysics()
{
	#pragma omp parallel for
	for(int ij=0; ij<metric.grid.n12; ij++)
	{	
		kappa0[ij] = 0;
		kappa1[ij] = 0;
		kappaA[ij] = 0;
		eta   [ij] = 0;
	}
}



template<class Coord>
void RadiationUniform<Coord>::ComputeMomentsIF()
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
	}
}

template<class Coord>
void RadiationUniform<Coord>::ComputeMomentsLF()
{
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		Tensor3x3<Coord,IF> EnergyMomentumTensorIF(E[ij],Fx[ij],Fy[ij], Fx[ij],Pxx[ij],Pxy[ij], Fy[ij],Pxy[ij],Pyy[ij]);
		Tensor3x3<Coord,LF> EnergyMomentumTensorLF = EnergyMomentumTensorIF.template Transform<LF>(metric.GetTetrad(ij));
		E_LF[ij] = EnergyMomentumTensorLF[{0,0}];
		Fx_LF[ij] = EnergyMomentumTensorLF[{0,1}];
		Fy_LF[ij] = EnergyMomentumTensorLF[{0,2}];
		Pxx_LF[ij] = EnergyMomentumTensorLF[{1,1}];
		Pxy_LF[ij] = EnergyMomentumTensorLF[{1,2}];
		Pyy_LF[ij] = EnergyMomentumTensorLF[{2,2}];
		Fnorm_LF[ij] = Tensor2<Coord,LF>(Fx_LF[ij],Fy_LF[ij]).Norm(metric.GetGamma_ll(ij));
	}
}



template<class Coord>
Coordinate2<Coord> RadiationUniform<Coord>::GetTempCoordinate(int i, int j, double phi)
{
	std::vector<double> cX1(fourierStencil.nDir);
	std::vector<double> cX2(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
	{
		cX1[f] = coefficientsX1[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
		cX2[f] = coefficientsX2[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	}
	return Coordinate2<Coord>(Fourier::Expansion::GetValue(phi,cX1), Fourier::Expansion::GetValue(phi,cX2));
}
template<class Coord>
Tensor2<Coord,LF> RadiationUniform<Coord>::GetTemp2Velocity(int i, int j, double phi)
{
	std::vector<double> cCx1(fourierStencil.nDir);
	std::vector<double> cCx2(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
	{
		cCx1[f] = coefficientsCx1[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
		cCx2[f] = coefficientsCx2[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	}
	return Tensor2<Coord,LF>(Fourier::Expansion::GetValue(phi,cCx1), Fourier::Expansion::GetValue(phi,cCx2));
}
template<class Coord>
double RadiationUniform<Coord>::GetFrequencyShift(int i, int j, double phi)
{
	std::vector<double> cS(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
		cS[f] = coefficientsS[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	return Fourier::Expansion::GetValue(phi,cS);
}
template<class Coord>
double RadiationUniform<Coord>::InterpolateIntensity(int ij, double phi)
{
	double k = stencil.Index(phi);

	int kFloor = floor(k);
	int kM1 = (kFloor - 1 + nDir) % nDir;
	int kP0 = (kFloor + 0 + nDir) % nDir;
	int kP1 = (kFloor + 1 + nDir) % nDir;
	int kP2 = (kFloor + 2 + nDir) % nDir;

	return CubicInterpolation (k - kFloor,I[ij + kM1*grid.n12],I[ij + kP0*grid.n12],I[ij + kP1*grid.n12],I[ij + kP2*grid.n12]);
	// return LinearInterpolation(k - kFloor,I[ij + kP0*grid.n12],I[ij + kP1*grid.n12]);
}



template<class Coord>
void RadiationUniform<Coord>::Stream()
{
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Skip LPs which are inside BH:
		if(metric.InsideBH(i,j))
		{
			temp[ijd] = 0;
			continue;
		}

		// Fourier Streaming:
		double phi = stencil.Phi(d);
		double s = GetFrequencyShift(i,j,phi);
		Coordinate2<Coord> xTemp = GetTempCoordinate(i,j,phi);
		Tensor2<Coord,LF> vTempLF = GetTemp2Velocity(i,j,phi);
		Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
		Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);

		// Skip temporary Grid Points inside BH:
		if(metric.InsideBH(xTemp))
		{
			temp[ijd] = 0;
			continue;
		}

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Intensity interpolation:
		double angle = vTempIFxy.Angle();
		double alpha = metric.GetAlpha(ij);
		double intensityAt_i0j0 = pow(alpha / metric.GetAlpha(grid.Index(i0,j0)),4) * InterpolateIntensity(grid.Index(i0,j0),angle);
		double intensityAt_i0j1 = pow(alpha / metric.GetAlpha(grid.Index(i0,j1)),4) * InterpolateIntensity(grid.Index(i0,j1),angle);
		double intensityAt_i1j0 = pow(alpha / metric.GetAlpha(grid.Index(i1,j0)),4) * InterpolateIntensity(grid.Index(i1,j0),angle);
		double intensityAt_i1j1 = pow(alpha / metric.GetAlpha(grid.Index(i1,j1)),4) * InterpolateIntensity(grid.Index(i1,j1),angle);

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
		 intensityAt_i0j0,intensityAt_i0j1,
		 intensityAt_i1j0,intensityAt_i1j1);
	}
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}



template<class Coord>
void RadiationUniform<Coord>::Collide()
{
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
void RadiationUniform<Coord>::RunSimulation(int writeFrequency, bool updateFourierHarmonics, bool keepSourceNodesActive, bool writeData, bool printToTerminal)
{
	// -------------------- Initialization --------------------
	Timer<10>::Reset("MainLoop");
	Timer<20>::Reset("ComputeMomentsIF");
	Timer<30>::Reset("WriteDataToJson");
	Timer<40>::Reset("Collide");
	Timer<50>::Reset("UpdateFourierCoefficients");
	Timer<60>::Reset("Stream");
	Timer<70>::Reset("KeepSourceNodesActive");
	Timer<80>::Reset("Transform Moments");

	// Initial data output:
	if (printToTerminal)
	{
		std::cout << " sigma=" << simData.sigma << std::endl;
		std::cout << " n1=" << grid.n1 << std::endl;
		std::cout << " n2=" << grid.n2 << std::endl;
		std::cout << " nDir=" << nDir << std::endl;
		std::cout << " nFourier=" << simData.nFourier << std::endl;
		std::cout << " simTime=" << simData.simTime << std::endl;
		std::cout << " dt=" << grid.dt << std::endl;
		std::cout << " timeSteps=" << simData.timeSteps << std::endl;
		std::cout << " filesToWrite=" << simData.timeSteps/writeFrequency << std::endl;
	}
	// --------------------------------------------------------



	// ----------------- Main simulation Loop -----------------
	Timer<10>::Start();
	for(int n=0; n<simData.timeSteps; n++)
	{
		Timer<20>::Start();
		ComputeMomentsIF();
		Timer<20>::End();
		Timer<40>::Start();
		Collide();
		Timer<40>::End();

		if (writeData && (n % writeFrequency) == 0)
		{
			Timer<80>::Start();
			ComputeMomentsLF();
			Timer<80>::End();
			Timer<30>::Start();
			grid.WriteFrametoJson(n*grid.dt,E_LF,Fx_LF,Fy_LF,Fnorm_LF,n,simData.directoryPath + "/Moments");
			Timer<30>::End();
		}


		if (updateFourierHarmonics)
		{
			Timer<50>::Start();
			UpdateFourierCoefficients();
			Timer<50>::End();
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
	if (writeData)
	{
		Timer<30>::Start();
		grid.WriteFrametoJson(simData.timeSteps*grid.dt,E_LF,Fx_LF,Fy_LF,Fnorm_LF,simData.timeSteps,simData.directoryPath + "/Moments");
		Timer<30>::End();
	}
	Timer<10>::End();

	if (printToTerminal)
	{
		Timer<10>::Print();
		Timer<30>::Print();
	}
	// --------------------------------------------------------



	// ---------------------- Termination ---------------------
	simData.AddTimeMeasurement(Timer<10>::name, Timer<10>::time);
	simData.AddTimeMeasurement(Timer<20>::name, Timer<20>::time);
	simData.AddTimeMeasurement(Timer<30>::name, Timer<30>::time);
	simData.AddTimeMeasurement(Timer<40>::name, Timer<40>::time);
	simData.AddTimeMeasurement(Timer<50>::name, Timer<50>::time);
	simData.AddTimeMeasurement(Timer<60>::name, Timer<60>::time);
	simData.AddTimeMeasurement(Timer<70>::name, Timer<70>::time);
	simData.AddTimeMeasurement(Timer<80>::name, Timer<80>::time);
	simData.LogSimulationParameters();
	// --------------------------------------------------------
}



template class RadiationUniform<xy>;
template class RadiationUniform<rph>;