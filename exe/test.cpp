#include <iostream>
#include "../src/Radiation.h"
using namespace std;

void SphereWaveShadowRegion()
{
    Coord p0 = Coord(0.5732233, 0.9267766);
    Coord p1 = Coord(0.9267766, 0.5732233);
    double rOut = 1.80;
    double rIn = p0.EuklNorm();
    double phi0 = M_PI / 2.0;
    double phi1 = p0.Phi();
    double phi2 = p1.Phi();
    double phi3 = 0;

    int n = 20;
    cout << "#x, y, z" << endl;
    for (int i = 0; i < 20; i++)
    {
        double t = i / (n - 1.0);
        double phi = (phi1 - phi0) * t + phi0;
        double x = rOut * MyCos(phi);
        double y = rOut * MySin(phi);
        cout << x << ", " << y << ", 0" << endl;
    }
    cout << p0[1] << ", " << p0[2] << ", 0" << endl;
    cout << p1[1] << ", " << p1[2] << ", 0" << endl;

    for (int i = 0; i < 20; i++)
    {
        double t = i / (n - 1.0);
        double phi = (phi3 - phi2) * t + phi2;
        double x = rOut * MyCos(phi);
        double y = rOut * MySin(phi);
        cout << x << ", " << y << ", 0" << endl;
    }
    cout << endl;
}

void MomentsTransformation()
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 101;
    Coord start(2, 2);
    Coord end(3, 3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Coord xy(2.5, 2.5);
    Tensor3x3 tetrad = metric.GetTetrad(xy);
    Tensor3x3 g_ll = metric.GetMetric_ll(xy);
    Tensor3x3 eta_ll = metric.GetMinkowskiMetric_ll(xy);
    Tensor3x3 eta_ll_(0.0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int I = 0; I < 3; I++)
                for (int J = 0; J < 3; J++)
                    eta_ll_[{I, J}] += tetrad[{i, I}] * tetrad[{j, J}] * g_ll[{i, j}];

    tetrad.Print("tetrad", true);
    g_ll.Print("  g_ll", true);
    eta_ll.Print("eta_ll", true);
    eta_ll_.Print("eta_ll");

    double E_ = 10;
    double Fx_ = 8;
    double Fy_ = 2;
    double Pxx_ = 1;
    double Pxy_ = 2;
    double Pyy_ = 3;
    Tensor3x3 T_(E_, Fx_, Fy_,
                 Fx_, Pxx_, Pxy_,
                 Fy_, Pxy_, Pyy_);

    Tensor3x3 T(0.0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int I = 0; I < 3; I++)
                for (int J = 0; J < 3; J++)
                    T[{i, j}] += tetrad[{i, I}] * tetrad[{j, J}] * T_[{I, J}];

    double E = IntegerPow<2>(tetrad[{0, 0}]) * E_;
    double Fx = tetrad[{1, 0}] * tetrad[{0, 0}] * E_ + tetrad[{0, 0}] * (tetrad[{1, 1}] * Fx_ + tetrad[{1, 2}] * Fy_);
    double Fy = tetrad[{2, 0}] * tetrad[{0, 0}] * E_ + tetrad[{0, 0}] * (tetrad[{2, 1}] * Fx_ + tetrad[{2, 2}] * Fy_);

    cout << "E:  " << T[{0, 0}] << ", " << E << endl;
    cout << "Fx: " << T[{0, 1}] << ", " << Fx << endl;
    cout << "Fy: " << T[{0, 2}] << ", " << Fy << endl;

    double E_new = E / IntegerPow<2>(tetrad[{0, 0}]);
    double A = Fx / tetrad[{0, 0}] - tetrad[{1, 0}] * E_new;
    double B = Fy / tetrad[{0, 0}] - tetrad[{2, 0}] * E_new;
    double a = tetrad[{1, 1}];
    double b = tetrad[{1, 2}];
    double c = tetrad[{2, 1}];
    double d = tetrad[{2, 2}];
    double Fx_new = (A * d - B * b) / (a * d - c * b);
    double Fy_new = (A * c - B * a) / (b * c - d * a);

    cout << "E_:  " << E_ << ", " << E_new << endl;
    cout << "Fx_: " << Fx_ << ", " << Fx_new << endl;
    cout << "Fy_: " << Fy_ << ", " << Fy_new << endl;
}

void DistributionNormalization()
{
    Stencil stencil(1000);
    double E0 = 20;
    double Fx0 = 0.55 * E0 * MyCos(M_PI / 4);
    double Fy0 = 0.55 * E0 * MySin(M_PI / 4);

    double E = 0;
    double Fx = 0;
    double Fy = 0;

    double normalization;
    double sigma = -1;

    double currentF = 0;
    double targetF = Tensor2(Fx0, Fy0).EuklNorm();
    int refinement = 0;
    while (abs(currentF - targetF) / targetF > 0.001) // while difference bigger 0.1%
    {
        if (currentF < targetF)
            sigma += pow(10, -refinement);
        if (currentF > targetF)
        {
            sigma -= pow(10, -refinement);
            refinement++;
            currentF = 0;
            continue;
        }

        normalization = 0;
        for (int d = 0; d < stencil.nDir; d++)
        {
            double c = stencil.W(d) * exp(sigma * (MyCos(stencil.Phi(d)) - 1.0));
            normalization += c;
        }
        normalization = 1.0 / normalization;

        double currentFx = 0;
        double currentFy = 0;
        for (int d = 0; d < stencil.nDir; d++)
        {
            Tensor2 dir = stencil.C(d);
            double vonMisesAngle = Tensor2(Fx0, Fy0).Phi();

            double angle = fmod(stencil.Phi(d) - vonMisesAngle + 2.0 * M_PI, 2.0 * M_PI);
            double Id = E0 * exp(sigma * (MyCos(angle) - 1.0)) * normalization;
            double c = stencil.W(d) * Id;
            currentFx += c * dir[1];
            currentFy += c * dir[2];
        }
        currentF = Tensor2(currentFx, currentFy).EuklNorm();
        cout << "Test: " << Format(currentF, 4) << ", " << Format(Tensor2(Fx0, Fy0).EuklNorm(), 4) << ", " << sigma << endl;
    }

    for (int d = 0; d < stencil.nDir; d++)
    {
        Tensor2 dir = stencil.C(d);
        double vonMisesAngle = Tensor2(Fx0, Fy0).Phi();

        double angle = fmod(stencil.Phi(d) - vonMisesAngle + 2.0 * M_PI, 2.0 * M_PI);
        double Id = E0 * exp(sigma * (MyCos(angle) - 1.0)) * normalization;
        double c = stencil.W(d) * Id;
        E += c;
        Fx += c * dir[1];
        Fy += c * dir[2];
    }

    int acc = 6;
    cout << "sigma: " << sigma << endl;
    cout << "E:     " << Format(E0, acc) << ", " << Format(E, acc) << endl;
    cout << "Fx:    " << Format(Fx0, acc) << ", " << Format(Fx, acc) << endl;
    cout << "Fy:    " << Format(Fy0, acc) << ", " << Format(Fy, acc) << endl;
    cout << "F:     " << Format(Tensor2(Fx0, Fy0).EuklNorm(), acc) << ", " << Format(Tensor2(Fx, Fy).EuklNorm(), acc) << endl;
    cout << "angle: " << Format(Tensor2(Fx0, Fy0).Phi(), acc) << ", " << Format(Tensor2(Fx, Fy).Phi(), acc) << endl;
}

void GrammSchmidt()
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 101;
    Coord start(2, 2);
    Coord end(3, 3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Coord xy(2.5, 2.5);
    metric.GetMetric_ll(xy).Print("  g_ll", true);

    Profiler::Session &session = Profiler::Session::Get();
    session.Start("config.name", "output/testProfileResults.json");
    {
        PROFILE_SCOPE("Gramm Schmidt");
        for (int i = 0; i < 10000000; i++)
        {
            Tensor3x3 g_ll = metric.GetMetric_ll(xy);

            Tensor3 v0(1, 0, 0);
            Tensor3 u0 = v0;
            double u0u0 = Dot(u0, u0, g_ll);
            Tensor3 e0 = u0 / sqrt(abs(u0u0));

            Tensor3 v1(0, 0, 1);
            Tensor3 u1 = v1 - u0 * Dot(v1, u0, g_ll) / u0u0;
            double u1u1 = Dot(u1, u1, g_ll);
            Tensor3 e1 = u1 / sqrt(abs(u1u1));

            Tensor3 v2(0, 1, 0);
            Tensor3 u2 = v2 - u0 * Dot(v2, u0, g_ll) / u0u0 - u1 * Dot(v2, u1, g_ll) / u1u1;
            double u2u2 = Dot(u2, u2, g_ll);
            Tensor3 e2 = u2 / sqrt(abs(u2u2));

            Tensor3x3 tetradGS(e0[0], e1[0], e2[0],
                               e0[1], e1[1], e2[1],
                               e0[2], e1[2], e2[2]);
        }

        // Printing:
        // Tensor3x3 etaGS(0.0);
        // for (int i = 0; i < 3; i++)
        //    for (int j = 0; j < 3; j++)
        //        for (int I = 0; I < 3; I++)
        //            for (int J = 0; J < 3; J++)
        //                etaGS[{I, J}] += tetradGS[{i, I}] * tetradGS[{j, J}] * g_ll[{i, j}];
        // tetradGS.Print("tetradGS", true);
        // etaGS.Print("etaGS", true);
    }
    {
        PROFILE_SCOPE("My Method");
        for (int i = 0; i < 10000000; i++)
        {
            Tensor3 n1 = metric.uEulObs(xy);
            Tensor2x2 g1 = metric.GetGamma_ll(xy);
            Tensor2x2 g1Inv = metric.GetGamma_uu(xy);
            double a = 1.0 / sqrt(g1Inv[{1, 1}]);
            double b = a * a * g1Inv[{1, 2}];
            Tensor2x2 m1(1 / a, 0, b / a, 1);
            double g2 = m1[{2, 2}] * g1[{2, 2}] * m1[{2, 2}];

            Tensor3x3 tetradME(n1[0], 0, 0,
                               n1[1], m1[{1, 1}], 0,
                               n1[2], m1[{2, 1}], 1 / sqrt(g2));
        }

        // Printing:
        // Tensor3x3 g_ll = metric.GetMetric_ll(xy);
        // tetradME.Print("tetrad ME", true);
        // Tensor3x3 etaME(0.0);
        // for (int i = 0; i < 3; i++)
        //    for (int j = 0; j < 3; j++)
        //        for (int I = 0; I < 3; I++)
        //            for (int J = 0; J < 3; J++)
        //                etaME[{I, J}] += tetradME[{i, I}] * tetradME[{j, J}] * g_ll[{i, j}];
        // etaME.Print("etaME");
    }
    session.End();
    std::vector<std::string> names = session.GetAllFunctionNames();
    for (int i = 0; i < names.size(); i++)
        session.PrintFunctionDuration(names[i]);
}

double Dy(double y)
{
    return -15 * y;
}
void LambdaIteration()
{
    double dt = 0.25;
    double y0 = 1;

    std::vector<double> solution;
    solution.push_back(y0);
    double t = 0;
    double y = y0;
    while (t < 1)
    {
        double dy = Dy(y);
        y += dy * dt;
        solution.push_back(y);
        t += dt;
    }

    PrintList(solution, "solution");
}

void UnitConversion()
{
    double centralMass = 1.0;     // mass of central object in solar masses.
    double mSolar = 1.4776 * 1e5; // solar mass in geometrized cgs units.

    double massCgsToCode = 1.0 / mSolar * centralMass;
    double massCodeToCgs = 1.0 / massCgsToCode;
    double lengthCgsToCode = 1.0 / mSolar * centralMass;
    double lengthCodeToCgs = 1.0 / lengthCgsToCode;
    double timeCgsToCode = 1.0 / mSolar * centralMass;
    double timeCodeToCgs = 1.0 / timeCgsToCode;

    double energyDensityCgsToCode = massCgsToCode / (lengthCgsToCode * timeCgsToCode * timeCgsToCode);
    double energyDensityCodeToCgs = 1.0 / energyDensityCgsToCode;
    double fluxDensityCgsToCode = massCgsToCode / (timeCgsToCode * lengthCgsToCode * lengthCgsToCode);
    double fluxDensityCodeToCgs = 1.0 / fluxDensityCgsToCode;

    cout << "1 cgs energy density in code units: " << energyDensityCgsToCode << endl;
    cout << "1 cgs flux density in code units: " << fluxDensityCgsToCode << endl;
}

void LorentzBoostToFluidFrame()
{
    Tensor2 uFluid(0.5, 0.5);

    double betaSq = Tensor2::Dot(uFluid, uFluid);
    double gamma = 1.0 / sqrt(1.0 - betaSq);
    Tensor3x3 boost(gamma, -gamma * uFluid[1], -gamma * uFluid[2],
                    -gamma * uFluid[1], 1 + (gamma - 1) * uFluid[1] * uFluid[1] / betaSq, (gamma - 1) * uFluid[1] * uFluid[2] / betaSq,
                    -gamma * uFluid[2], (gamma - 1) * uFluid[1] * uFluid[2] / betaSq, 1 + (gamma - 1) * uFluid[2] * uFluid[2] / betaSq);
    Tensor3x3 boostInv = boost.Invert();

    double E = 1;
    double Fx = 0.5;
    double Fy = 0.5;
    double Pxx = 0.25;
    double Pxy = 0.125;
    double Pyy = 0.25;
    Tensor3x3 T_LF(E, Fx, Fy,
                   Fx, Pxx, Pxy,
                   Fy, Pxy, Pyy);
    Tensor3x3 T_FF(0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int I = 0; I < 3; I++)
                for (int J = 0; J < 3; J++)
                    T_FF[{I, J}] += T_LF[{i, j}] * boost[{I, i}] * boost[{J, j}];

    T_LF.Print("T_LF");
    T_FF.Print("T_FF");
}

void NullVelocityTransformations()
{
    size_t nx, ny;
    nx = ny = 101;
    Coord start(2, 2);
    Coord end(3, 3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Coord xy(2.4, 2.6);

    Tensor2 vIF = Tensor2(0.1, 0.2).EuklNormalized();
    Tensor3 uIF = metric.GetAlpha(xy) * Tensor3(1, vIF[1], vIF[2]);
    Tensor3x3 tetrad = metric.GetTetrad(xy);
    Tensor3 uLF = TransformIFtoLF(uIF, tetrad);

    vIF.Print("vIF");
    PrintDouble(vIF.EuklNorm(), "|vIF|");
    uIF.Print("uIF");
    PrintDouble(Norm2(uIF,metric.GetMinkowskiMetric_ll(xy)), "|uIF|");

    uLF.Print("uLF");
    PrintDouble(Norm2(uLF, metric.GetMetric_ll(xy)), "|uLF|");

    vIF = Vec2ObservedByEulObs<LF, IF>(uLF, xy, metric);
    vIF.Print("vIF");
    PrintDouble(vIF.EuklNorm(), "|vIF|", true);

    Coord xy2(2.41, 2.61);
    Tensor3 uLF2 = NullNormalize(uLF, metric.GetMetric_ll(xy2));
    Tensor2 vIF2 = Vec2ObservedByEulObs<LF, IF>(uLF2, xy2, metric);
    
    uLF2.Print("uLF2");
    PrintDouble(Norm2(uLF2, metric.GetMetric_ll(xy2)), "|uLF2|");
    vIF2.Print("vIF2");
    PrintDouble(vIF2.EuklNorm(), "|vIF2|", true);
}

void UndefinedTest()
{
    // Partial test, can be deleted
    size_t nx, ny;
    nx = ny = 101;
    Coord start(2, 2);
    Coord end(3, 3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Coord xy0(2.4, 2.6);
    Coord xy1(2.42, 2.62);

    Tensor3x3 tetrad0 = metric.GetTetrad(xy0);
    Tensor3x3 tetrad1 = metric.GetTetrad(xy1);
    Tensor3x3 invTetrad0 = tetrad0.Invert();
    Tensor3x3 invTetrad1 = tetrad1.Invert();

    Tensor2 vIF0 = Tensor2(0.1, 0.2).EuklNormalized();
    Tensor3 uIF0 = metric.GetAlpha(xy0) * Tensor3(1, vIF0[1], vIF0[2]);
    vIF0.Print("vIF0");
    PrintDouble(vIF0.EuklNorm(), "|vIF0|");
    uIF0.Print("uIF0");
    PrintDouble(Norm2(uIF0,metric.GetMinkowskiMetric_ll(xy0)), "|uIF0|");

    Tensor3 uLF = TransformIFtoLF(uIF0, tetrad0);

    Tensor3 uIF1 = TransformLFtoIF(uLF, tetrad1);
    Tensor2 vIF1 = Tensor2(uIF1[1], uIF1[2]) / metric.GetAlpha(xy1);
    vIF1.Print("vIF1");
    PrintDouble(vIF1.EuklNorm(), "|vIF1|");
    uIF1.Print("uIF1");
    PrintDouble(Norm2(uIF1,metric.GetMinkowskiMetric_ll(xy1)), "|uIF1|");

    Tensor3 n0 = Tensor3(1, -metric.GetBeta_u(xy0)[1], -metric.GetBeta_u(xy0)[2]) / metric.GetAlpha(xy0);
    Tensor3 n1 = Tensor3(1, -metric.GetBeta_u(xy1)[1], -metric.GetBeta_u(xy1)[2]) / metric.GetAlpha(xy1);
    n0.Print("n0");
    n1.Print("n1");
}

void BasisChange()
{
    size_t nx, ny;
    nx = ny = 101;
    Coord start(2, 2);
    Coord end(3, 3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.5);

    Coord xy0(2.4, 2.6);
    Coord xy1(2.42, 2.62);

    Tensor3x3 tetrad0 = metric.GetTetrad(xy0);
    Tensor3x3 tetrad1 = metric.GetTetrad(xy1);

    Tensor2x2 subTetrad0(tetrad0[{1,1}], tetrad0[{1,2}], tetrad0[{2,1}], tetrad0[{2,2}]);
    Tensor2x2 subTetrad1(tetrad1[{1,1}], tetrad1[{1,2}], tetrad1[{2,1}], tetrad1[{2,2}]);

    tetrad0.Print("tetrad0");
    tetrad1.Print("tetrad1");
    subTetrad0.Print("subTetrad0");
    subTetrad1.Print("subTetrad1");

    Tensor2x2 T = subTetrad1.Invert() * subTetrad0;
    T.Print("T");

    Tensor2 v0 = Tensor2(0.3,0.5).EuklNormalized();
    Tensor2 v1 = T * v0;
    v0.Print("v0");
    v1.Print("v1");
    PrintDouble(v1.EuklNorm(), "|v1|");
    v1.EuklNormalized().Print("v1");

    Tensor2x2 A(3,  1, 1, 3);
    Tensor2x2 B(3, -1, 0, 3);
    Tensor2x2 C = B.Invert() * A;
    Tensor2 vA(1,1);
    Tensor2 vB = C * vA;
    A.Print("A");
    B.Print("B");
    C.Print("C");
    vA.Print("vA");
    vB.Print("vB");
}

void Rank2TensorRotation()
{
    Tensor2x2 T(1,0, 0,2);
    double angle = M_PI / 4.0;
    double s = MySin(angle);
    double c = MyCos(angle);
    Tensor2x2 R(c, -s, s, c);

    Tensor2x2 Trot;
    for(int i=1; i<3; i++)
    for(int j=1; j<3; j++)
    {
        Trot[{i,j}] = 0;
        for(int k=1; k<3; k++)
        for(int l=1; l<3; l++)
            Trot[{i,j}] += R[{i,k}] * R[{j,l}] * T[{k,l}];
    }

    T.Print("T");
    R.Print("R");
    Trot.Print("Trot");
}

int main()
{
    // SphereWaveShadowRegion();
    // MomentsTransformation();
    // DistributionNormalization();
    // GrammSchmidt();
    // LambdaIteration();
    // UnitConversion();
    // LorentzBoostToFluidFrame();
    // NullVelocityTransformations();
    // UndefinedTest();
    // BasisChange();
    Rank2TensorRotation();
    

    // for(int i = 0; i < 16; i++)
    // {
        // double phi = i * 2.0 * M_PI / 16;
        // double x = MyCos(phi);
        // double y = MySin(phi);
        // 
        // cout << phi << ", " << fmod(MyAtan2(y, x) + 2.0 * M_PI, 2.0 * M_PI) << endl;
    // }
}