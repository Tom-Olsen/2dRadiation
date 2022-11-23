#include "../src/Stencil.hh"
#include "../src/Interpolation.h"
#include "../src/Utility.h"

using namespace std;



void Test_DtoXmapping()
{
    constexpr int nDir = 10;
    // UniformStencil<nDir>& stencil = UniformStencil<nDir>::GetInstance();
    DirectedStencil<nDir>& stencil = DirectedStencil<nDir>::GetInstance();
    double rot = 1;

    cout << "Output should be between 0 and 1 (inclusive)." << endl;
    for(double d=-0.5; d<=nDir-0.5; d+=0.5)
        cout << "X(" + Format(d,1) + ") =" << Format(stencil.X(d),3) << endl;

    cout << endl;
    cout << "d at start and end should be the same:" << endl;
    for(double d=0; d<nDir; d++)
        cout << "d=" << Format(d,1) << ", " << "phi=" << Format(stencil.Phi(d,rot),3) << ", d=" << Format(stencil.Index(stencil.Phi(d,rot),rot),1) << endl;
}



void Test_Weights()
{
    constexpr int nDir = 20;
    UniformStencil<nDir>& uniStencil = UniformStencil<nDir>::GetInstance();
    DirectedStencil<nDir>& dirStencil = DirectedStencil<nDir>::GetInstance();

    uniStencil.Print();
    dirStencil.Print();

    int formatAcc = 16;
    {
        cout << "uniform stencil:" << endl;
        int N = nDir/2;
        int M = uniStencil.accuracy + 4;
        double I[M];
        I[0] = 1.0;
        for(int d=1; d<M; d++)
            I[d] = (pow(-1,d) + 1) * tgamma((d+1.0)/2.0) / (sqrt(M_PI)*d*tgamma(d/2.0));
    
        double w[N];
        for(int d=0; d<N; d++)
            w[d] = uniStencil.W(d);
    
        double A[N*M];
        for(int i=0; i<M; i++)
            for(int j=0; j<N; j++)
                A[i*N+j] = pow(cos(uniStencil.Phi(j)),i);
    
        double Itest[M];
        for(int i=0; i<M; i++)
        {
            Itest[i] = 0;
            for(int j=0; j<N; j++)
                Itest[i] += A[i*N+j] * w[j];
        }
    
        for(int d=0; d<M; d++)
        {
            if(d<=uniStencil.accuracy)
                cout << Format(d,0) << ": I, Itest = " << "(" << Format(I[d],formatAcc) << "," << Format(Itest[d],formatAcc) << ")\n";
            else 
                cout << Format(d,0) << ": I, Itest = " << "(" << Format(I[d],formatAcc) << "," << Format(Itest[d],formatAcc) << ") <- should be wrong!\n";
        }
    }
 cout << endl;
    {
        cout << "directed stencil:" << endl;
        int N = nDir/2;
        int M = dirStencil.accuracy + 4;
        double I[M];
        I[0] = 1.0;
        for(int d=1; d<M; d++)
            I[d] = (pow(-1,d) + 1) * tgamma((d+1.0)/2.0) / (sqrt(M_PI)*d*tgamma(d/2.0));
    
        double w[N];
        for(int d=0; d<N; d++)
            w[d] = dirStencil.W(d);
    
        double A[N*M];
        for(int i=0; i<M; i++)
            for(int j=0; j<N; j++)
                A[i*N+j] = pow(cos(dirStencil.Phi(j)),i);
    
        double Itest[M];
        for(int i=0; i<M; i++)
        {
            Itest[i] = 0;
            for(int j=0; j<N; j++)
                Itest[i] += A[i*N+j] * w[j];
        }
    
        for(int d=0; d<M; d++)
        {
            if(d<=dirStencil.accuracy)
                cout << Format(d,0) << ": I, Itest = " << "(" << Format(I[d],formatAcc) << "," << Format(Itest[d],formatAcc) << ")\n";
            else 
                cout << Format(d,0) << ": I, Itest = " << "(" << Format(I[d],formatAcc) << "," << Format(Itest[d],formatAcc) << ") <- should be wrong!\n";
        }
    }

    Double2(dirStencil.Phi(5),dirStencil.Phi(5,1)).Print("phi,phiRot");
}



void Test_Quadrature()
{
    constexpr int nDir = 200;
    UniformStencil<nDir>& uniStencil = UniformStencil<nDir>::GetInstance();
    DirectedStencil<nDir>& dirStencil = DirectedStencil<nDir>::GetInstance();

    auto f = [](double phi, int nFreq)
    {
        double value = 0;
        for(int k=0; k<nFreq; k++)
            value += k/10.0 * cos(k * phi);
        return value;
    };

    // Fill values:
    int highestFrequency = 20;
    double uniI[nDir];
	double dirI[nDir];
    for(int d=0; d<nDir; d++)
    {
        uniI[d] = f(uniStencil.Phi(d),highestFrequency);
        dirI[d] = f(dirStencil.Phi(d),highestFrequency);
    }

    // Extract frequencies:
    for(int k=0; k<highestFrequency; k++)
    {
        double uniFreq = 0;
        double dirFreq = 0;
        for(int d=0; d<nDir; d++)
        {
            uniFreq += uniI[d] * cos(k*uniStencil.Phi(d)) * uniStencil.W(d);
            dirFreq += dirI[d] * cos(k*dirStencil.Phi(d)) * dirStencil.W(d);
        }
        Double2(uniFreq,dirFreq).Print("Frequency " + to_string(k));
    }
}



void Test_InterpolationBetweenUniAndDirStencil()
{
    // allocation:
    constexpr int nDir = 100;
    UniformStencil<nDir>& uniStencil = UniformStencil<nDir>::GetInstance();
    DirectedStencil<nDir>& dirStencil = DirectedStencil<nDir>::GetInstance();
    double uniI[nDir];
    double dirI[nDir];
    double uniI2[nDir];

    // fill uniform stencil intensities:
    int dir0 = nDir/4;    // direction of peak
    double rotation = uniStencil.Phi(dir0);
    for(int k=-dir0; k<nDir - dir0; k++)
        uniI[(k+dir0)%nDir] = exp(-pow(k/2.0,2));
        
    int dir1 = 7*nDir/8;    // direction of peak
    for(int k=-dir1; k<nDir - dir1; k++)
        uniI[(k+dir1)%nDir] += exp(-pow(k/2.0,2));

    cout << "uniform distribution:" << endl;
    for(int k=0; k<nDir; k++)
        cout << uniI[k] << endl;
    cout << endl;
    
    // interpolate from uniform to directed stencil:
    {
        Stencil& stencilOld = uniStencil;
        Stencil& stencilNew = dirStencil;
        double rotationOld = -1;
        double rotationNew = rotation;
        for(int k=0; k<nDir; k++)
        {
            double phi = stencilNew.Phi(k,rotationNew);	// angle of current direction in new stencil
            double d   = stencilOld.Index(phi,rotationOld); // d value of above angle in old stencil

	        int dFloor = floor(d);
	        int m1 = (dFloor - 1 + nDir) % nDir;
	        int p0 = (dFloor + 0 + nDir) % nDir;
	        int p1 = (dFloor + 1 + nDir) % nDir;
	        int p2 = (dFloor + 2 + nDir) % nDir;

	        // dirI[k] = CubicInterpolation(d - dFloor, uniI[m1], uniI[p0], uniI[p1], uniI[p2]);
            dirI[k] = LinearInterpolation(d - dFloor, uniI[p0], uniI[p1]);
        }
    }
    
    cout << "directed distribution:" << endl;
    for(int k=0; k<nDir; k++)
        cout << dirI[k] << endl;
    cout << endl;



    // interpolate from directed to uniform stencil:
    {
        Stencil& stencilOld = dirStencil;
        Stencil& stencilNew = uniStencil;
        double rotationOld = rotation;
        double rotationNew = -1;
        for(int k=0; k<nDir; k++)
        {
            double phi = stencilNew.Phi(k,rotationNew); // angle of current direction in new stencil
            double d   = stencilOld.Index(phi,rotationOld); // d value of above angle in old stencil

	        int dFloor = floor(d);
	        int m1 = (dFloor - 1 + nDir) % nDir;
	        int p0 = (dFloor + 0 + nDir) % nDir;
	        int p1 = (dFloor + 1 + nDir) % nDir;
	        int p2 = (dFloor + 2 + nDir) % nDir;

	        // uniI2[k] = CubicInterpolation(d - dFloor, dirI[m1], dirI[p0], dirI[p1], dirI[p2]);
            uniI2[k] = LinearInterpolation(d - dFloor, dirI[p0], dirI[p1]);
        }
    }
    


    cout << "uniform2 distribution:" << endl;
    for(int k=0; k<nDir; k++)
        cout << uniI2[k] << endl;
    cout << endl;

    // calculate integrals
    double uniE = 0;
    double dirE = 0;
    double uniE2 = 0;
    double uniFx = 0;
    double dirFx = 0;
    double uniFx2 = 0;
    double uniFy = 0;
    double dirFy = 0;
    double uniFy2 = 0;
    for(int k=0; k<nDir; k++)
    {
        uniE   += uniStencil.W(k) * uniI[k];
        dirE   += dirStencil.W(k) * dirI[k];
        uniE2  += uniStencil.W(k) * uniI2[k];
        uniFx  += uniStencil.W(k) * uniI[k]  * uniStencil.Cx(k,rotation);
        dirFx  += dirStencil.W(k) * dirI[k]  * dirStencil.Cx(k,rotation);
        uniFx2 += uniStencil.W(k) * uniI2[k] * uniStencil.Cx(k,rotation);
        uniFy  += uniStencil.W(k) * uniI[k]  * uniStencil.Cy(k,rotation);
        dirFy  += dirStencil.W(k) * dirI[k]  * dirStencil.Cy(k,rotation);
        uniFy2 += uniStencil.W(k) * uniI2[k] * uniStencil.Cy(k,rotation);
    }

    Double3(uniE,dirE,uniE2).Print("Integrals");
    Double3(uniFx,dirFx,uniFx2).Print("Integrals");
    Double3(uniFy,dirFy,uniFy2).Print("Integrals");

    // for(int d=0; d<nDir; d++)
        // cout << "Phi: " << dirStencil.Phi(d) << ", " << dirStencil.Cxy(d).Angle() << endl;
}



void Test_WeightRotation()
{
    constexpr int nDir = 10;
    double rotation0 = M_PI/4;
    double rotation1 = M_PI/2;
    DirectedStencil<nDir>& stencil = DirectedStencil<nDir>::GetInstance();

    for(int d=0; d<nDir; d++)
        Double4(stencil.W(d), stencil.Cx(d,rotation0), stencil.Cy(d,rotation0), stencil.Phi(d,rotation0)).Print("w,cx,cy,phi");

    cout << endl << endl;

    for(int d=0; d<nDir; d++)
        Double4(stencil.W(d), stencil.Cx(d,rotation1), stencil.Cy(d,rotation1), stencil.Phi(d,rotation1)).Print("w,cx,cy,phi");
}



void Test_StencilInterpolation()
{
    constexpr int nDir = 50;
    DirectedStencil<nDir>& stencil = DirectedStencil<nDir>::GetInstance();

    // Grid arangement:
    // 6 7 8  - | /
    // 3 4 5  - / -
    // 0 1 2  / | |
    // angles:
    // -: 2.0 * M_PI
    // /: M_PI / 4.0
    // |: M_PI / 2.0

    double rotation[9];
    rotation[0] = M_PI / 4.0; // /
    rotation[1] = M_PI / 2.0; // |
    rotation[2] = M_PI / 2.0; // |
    rotation[3] = 2.0 * M_PI; // -
    rotation[4] = M_PI / 4.0; // /
    rotation[5] = 2.0 * M_PI; // -
    rotation[6] = 2.0 * M_PI; // -
    rotation[7] = M_PI / 2.0; // |
    rotation[8] = M_PI / 4.0; // /

	auto GaussPeak = [](double phi, double direction)
	{
        double sigma = 0.2;
        return exp(-0.5 * pow((phi - direction) / sigma, 2))
             + exp(-0.5 * pow((phi - direction + 2.0 * M_PI) / sigma, 2))
             + exp(-0.5 * pow((phi - direction - 2.0 * M_PI) / sigma, 2));
	};
    auto DoubleGaussPeak = [](double phi, double direction, double seperation)
    {
        double sigma = 0.2;
        return exp(-0.5 * pow((phi + seperation - direction) / sigma, 2))
             + exp(-0.5 * pow((phi + seperation - direction + 2.0 * M_PI) / sigma, 2))
             + exp(-0.5 * pow((phi + seperation - direction - 2.0 * M_PI) / sigma, 2))
             + exp(-0.5 * pow((phi - seperation - direction) / sigma, 2))
             + exp(-0.5 * pow((phi - seperation - direction + 2.0 * M_PI) / sigma, 2))
             + exp(-0.5 * pow((phi - seperation - direction - 2.0 * M_PI) / sigma, 2));
    };
    double I[9][nDir];
    for(int i=0; i<9; i++)
    {
        for(int d=0; d<nDir; d++)
        {
            double phi = stencil.Phi(d,rotation[i]);
            if(i==0 || i==4)
                I[i][d] = DoubleGaussPeak(phi,rotation[i],rotation[i]);
            else
                I[i][d] = GaussPeak(phi,rotation[i]);
        }
    }

    double E[9];
    double Fx[9];
    double Fy[9];
	constexpr double twoPiInv = 1.0 / (2.0 * M_PI);
    for(int i=0; i<9; i++)
    {
        E [i]  = 0.0;
	    Fx[i]  = 0.0;
	    Fy[i]  = 0.0;
	    for(int d = 0; d < nDir; d++)
	    {
	    	E [i]  += stencil.W(d) * I[i][d];
	    	Fx[i]  += stencil.W(d) * I[i][d] * stencil.Cx(d,rotation[i]);
	    	Fy[i]  += stencil.W(d) * I[i][d] * stencil.Cy(d,rotation[i]);
	    }
        E [i]  *= twoPiInv;
        Fx[i]  *= twoPiInv;
        Fy[i]  *= twoPiInv;
    }

    cout << "Momenta of the 9 original stencils:" << endl;
    Double3(E[6],E[7],E[8]).Print(" E(6,7,8)");
    Double3(E[3],E[4],E[5]).Print(" E(3,4,5)");
    Double3(E[0],E[1],E[2]).Print(" E(0,1,2)",true);

    Double3(Fx[6],Fx[7],Fx[8]).Print("Fx(6,7,8)");
    Double3(Fx[3],Fx[4],Fx[5]).Print("Fx(3,4,5)");
    Double3(Fx[0],Fx[1],Fx[2]).Print("Fx(0,1,2)",true);

    Double3(Fy[6],Fy[7],Fy[8]).Print("Fy(6,7,8)");
    Double3(Fy[3],Fy[4],Fy[5]).Print("Fy(3,4,5)");
    Double3(Fy[0],Fy[1],Fy[2]).Print("Fy(0,1,2)",true);


	auto InterpolatedIntensity = [](double phi, int i, double* rotation, double I[9][50])
    {
        constexpr int nDir = 50;
        DirectedStencil<nDir>& stencil = DirectedStencil<nDir>::GetInstance();
        float d = stencil.Index(phi,rotation[i]);
        Double2(phi,d).Print(to_string(i));

        int dFloor = (int)floor(d);
        int p0 = (dFloor + 0 + nDir) % nDir;
        int p1 = (dFloor + 1 + nDir) % nDir;

        return LinearInterpolation(d - dFloor, I[i][p0], I[i][p1]);
    };
    double rotationNew = M_PI/4.0;
    double Inew[nDir];
    for(int d=0; d<nDir; d++)
    {
        float phi = stencil.Phi(d,rotationNew);
        Coordinate2<xy> xTemp(-cos(phi), -sin(phi));

        float i = xTemp[1] + 1.0;
        float j = xTemp[2] + 1.0;
        int i0 = (int)floor(i);
        int j0 = (int)floor(j);
        int i1 = i0 + 1;
        int j1 = j0 + 1;

        cout << "(" << i << "," << j << ")\n"; 
        double Itemp[4];
        Itemp[0] = InterpolatedIntensity(phi, i0 + 3 * j0,rotation,I);
        Itemp[1] = InterpolatedIntensity(phi, i0 + 3 * j1,rotation,I);
        Itemp[2] = InterpolatedIntensity(phi, i1 + 3 * j0,rotation,I);
        Itemp[3] = InterpolatedIntensity(phi, i1 + 3 * j1,rotation,I);
        cout << endl;
        //Double4(Itemp[0],Itemp[1],Itemp[2],Itemp[3]).Print(to_string(d));
        Inew[d] = BilinearInterpolation(i - i0, j - j0, Itemp[0], Itemp[1], Itemp[2], Itemp[3]);
    }
    double Enew = 0;
    double Fxnew = 0;
    double Fynew = 0;
	for(int d = 0; d < nDir; d++)
	{
		Enew  += stencil.W(d) * Inew[d];
		Fxnew += stencil.W(d) * Inew[d] * stencil.Cx(d,rotationNew);
		Fynew += stencil.W(d) * Inew[d] * stencil.Cy(d,rotationNew);
	}
    Enew  *= twoPiInv;
    Fxnew *= twoPiInv;
    Fynew *= twoPiInv;

    cout << "Momenta of new stencils:" << endl;
    Double3(Enew,Fxnew,Fynew).Print("(E,Fx,Fy)new",true);

    cout << "Intensity comparicon to Unity:" << endl;
    for(int d=0; d<nDir; d++)
    {
        // PrintDouble(stencil.W(d),to_string(d));
        // PrintDouble(stencil.Phi(d,rotationNew),to_string(d));
        // PrintDouble(Inew[d],to_string(d));
        // PrintDouble(I[3][d],to_string(d));
    }
}



int main()
{
    // cout << endl << "--------------------------------------------------" << endl << endl;
    // Test_DtoXmapping();
    // cout << endl << "--------------------------------------------------" << endl << endl;
    // Test_Weights();
    // cout << endl << "--------------------------------------------------" << endl << endl;
    // Test_Quadrature();
    // cout << endl << "--------------------------------------------------" << endl << endl;
    // Test_InterpolationBetweenUniAndDirStencil();
    // cout << endl << "--------------------------------------------------" << endl << endl;
    // Test_WeightRotation();
    cout << endl << "--------------------------------------------------" << endl << endl;
    Test_StencilInterpolation();
    cout << endl << "--------------------------------------------------" << endl << endl;
}