#include "Stencil.h"



// --------------- InterpolationGrid ---------------
InterpolationGrid::InterpolationGrid(size_t nGrid, const Stencil& stencil) : nGrid(nGrid)
{
    neighbourIndexes.resize(nGrid);
    for(int i=0; i<nGrid; i++)
    {
        double phi = 2.0 * M_PI * i / nGrid;
        neighbourIndexes[i] = SurroundingNeighbourIndexes(phi, stencil);
    }

    neighbourWeights.resize(nGrid);
    for(int i=0; i<nGrid; i++)
    {
        double phi = 2.0 * M_PI * i / nGrid;
        neighbourWeights[i] = LagrangeWeights(phi, neighbourIndexes[i], stencil);
    }
}



double InterpolationGrid::d(double phi) const
{ return nGrid * phi / (2.0 * M_PI); }



std::array<size_t,4> InterpolationGrid::SurroundingNeighbourIndexes(double phi, const Stencil& stencil)
{
    // Find nearest grid index:
    int index = -1;
    double minDist = 10;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double delta = std::abs(phi-stencil.Phi(d));
        double dist = std::min(delta, 2.0 * M_PI - delta);
        if(dist < minDist)
        {
            minDist = dist;
            index = d;
        }
    }

    // Check whether angle phiMin is to the left or right of phi.
    // To resolve jumping angles close to phi=0, phiMin is rotated such that phi->M_PI.
    double phiMin = std::fmod(stencil.Phi(index) + M_PI - phi, 2.0 * M_PI);
    std::array<size_t,4> indices;
    if(phiMin < M_PI)
    {
        indices[0] = (index - 1 + stencil.nDir) % stencil.nDir;
        indices[1] = (index + 0 + stencil.nDir) % stencil.nDir;
        indices[2] = (index + 1 + stencil.nDir) % stencil.nDir;
        indices[3] = (index + 2 + stencil.nDir) % stencil.nDir;
    }
    else
    {
        indices[0] = (index - 2 + stencil.nDir) % stencil.nDir;
        indices[1] = (index - 1 + stencil.nDir) % stencil.nDir;
        indices[2] = (index + 0 + stencil.nDir) % stencil.nDir;
        indices[3] = (index + 1 + stencil.nDir) % stencil.nDir;
    }
    return indices;
}



std::array<double,4> InterpolationGrid::LagrangeWeights(double phi, std::array<size_t,4> neighbours, const Stencil& stencil)
{
    // To resolve jumping angles close to phi=0, angles are rotated such that phi->M_PI.
    std::array<double,4> weights;
    for(int i=0; i<4; i++)
    {
        double weight = 1;
        double phi_i = std::fmod(stencil.Phi(neighbours[i]) + M_PI - phi + 2.0 * M_PI, 2.0 * M_PI);
        for(int j=0; j<4; j++)
        {
            double phi_j = std::fmod(stencil.Phi(neighbours[j]) + M_PI - phi + 2.0 * M_PI, 2.0 * M_PI);
            if(i!=j)
                weight *= (M_PI - phi_j) / (phi_i - phi_j);
        }
        weights[i] = weight;
    }
    return weights;
}



void InterpolationGrid::Print() const
{
    std::cout << "        d,\tneighbours,\tweights,\n";
    for(size_t d=0; d<nGrid; d++)
    {
        std::cout << Format(d) << ",\t(" << neighbourIndexes[d][0] << "," << neighbourIndexes[d][1] << "," << neighbourIndexes[d][2] << "," << neighbourIndexes[d][3] << ")";
        std::cout << ",\t(" << Format(neighbourWeights[d][0]) << "," << Format(neighbourWeights[d][1]) << "," << Format(neighbourWeights[d][2]) << "," << Format(neighbourWeights[d][3]) << "),\n";
    }
    std::cout << std::endl;
}
// -------------------------------------------------



// -------------------- Stencil --------------------
Stencil::Stencil(size_t nOrder, int nGhost)
{
    name = "Stencil" + std::to_string(nOrder) + "." + std::to_string(nGhost);
    this->nDir = nOrder + nGhost;
    this->nGhost = nGhost;
    this->nOrder = nOrder;
    this->nCoefficients = nOrder;
    
    w.resize(nDir);
    phi.resize(nDir);
    cx.resize(nDir);
    cy.resize(nDir);
    
    // Add fourier directions:
    for(int d=0; d<nDir-nGhost; d++)
    {
        w[d] = 2.0 / (nDir - nGhost);
        phi[d] = 2.0 * M_PI * (d + 0.5) / (nDir - nGhost);
        cx[d] = MyCos(phi[d]);
        cy[d] = MySin(phi[d]);
    }

    if (nGhost > 0)
    {
        if(nOrder % 2 != 0)
            ExitOnError("Stencil odd nOrder not allowed when using nGhost>0.");

        // Find Buckets:
        double deltaPhi = M_PI / 8.0 + 1e-8;
        std::vector<double> bucketStart;
        std::vector<double> bucketEnd;
        std::vector<int> bucketCount;
        for(int d=0; d<nDir-nGhost; d++)
        {
            double phi0 = phi[d];
            double phi1 = phi[std::fmod(d + 1, nDir - nGhost)];
            double delta0 = std::abs(phi0 - M_PI);
            double delta1 = std::abs(phi1 - M_PI);
            double dist0 = std::min(delta0, 2.0 * M_PI - delta0);
            double dist1 = std::min(delta1, 2.0 * M_PI - delta1);
            if (dist0 <= deltaPhi && dist1 <= deltaPhi)
            {
                bucketStart.push_back(phi0);
                bucketEnd.push_back(phi1);
                bucketCount.push_back(0);
            }
        }

        // Start of first and end of last bucket:
        double minAngle = bucketStart[0];
        double maxAngle = bucketEnd[bucketEnd.size() - 1];

        // Fill Buckets:
        auto Distribution = [](double x)
        { return std::acos(1.0 - 2.0 * x) / M_PI; };
        for(int i=0; i<nGhost; i++)
        {
            double x = (i+0.5) / nGhost;
            double D = Distribution(x);
            double angle = minAngle + (maxAngle - minAngle) * D;

            for(int j=0; j<bucketCount.size(); j++)
            {
                double bucketCenter = (bucketStart[j] + bucketEnd[j]) / 2.0;
                double bucketSize = bucketEnd[j] - bucketStart[j];

                double delta = std::abs(angle - bucketCenter);
                double dist = std::min(delta, 2.0 * M_PI - delta);

                if(dist < bucketSize / 2.0)
                    bucketCount[j]++;
            }
        }

        // Add angles uniformly in buckets:
        double phiGhost[nGhost];
        int index = 0;
        for(int i=0; i<bucketCount.size(); i++)
        for(int j=0; j<bucketCount[i]; j++)
        {
            phiGhost[index] = bucketStart[i] + (bucketEnd[i] - bucketStart[i]) * (j + 0.5) / bucketCount[i] + M_PI;
            index++;
        }

        // Add ghost directions to stencil:
        index = nDir - nGhost;
        for(int d=0; d<nGhost; d++)
        {
            w[index] = 0;
            phi[index] = fmod(phiGhost[d] + 2.0 * M_PI, 2.0 * M_PI);
            cx[index] = MyCos(phi[index]);
            cy[index] = MySin(phi[index]);
            index++;
        }
    }
    
    SortDirections();
    interpolationGrid = InterpolationGrid(10 * nOrder, *this);
}



void Stencil::SortDirections()
{
    size_t index[nDir];
    for(size_t i = 0; i < nDir; i++)
        index[i] = i;

    // Sort index array based on custom compare function:
    std::sort(index, index + nDir, [this](size_t i, size_t j)
    { return phi[i] < phi[j]; });

    // Create temporary arrays to hold sorted values:
    double sorted_w[nDir];
    double sorted_phi[nDir];
    double sorted_cx[nDir];
    double sorted_cy[nDir];

    // Copy values from original arrays to temporary arrays:
    for(size_t i = 0; i < nDir; i++)
    {
        sorted_w[i] = w[index[i]];
        sorted_phi[i] = phi[index[i]];
        sorted_cx[i] = cx[index[i]];
        sorted_cy[i] = cy[index[i]];
    }

    // Overwrite original arrays with sorted values:
    for(size_t i = 0; i < nDir; i++)
    {
        w[i] = sorted_w[i];
        phi[i] = sorted_phi[i];
        cx[i] = sorted_cx[i];
        cy[i] = sorted_cy[i];
    }
}



double Stencil::W(size_t d) const
{ return w[d]; }
double Stencil::Phi(size_t d) const
{ return phi[d]; }
double Stencil::Cx(size_t d) const
{ return cx[d]; }
double Stencil::Cy(size_t d) const
{ return cy[d]; }
Tensor2 Stencil::C(size_t d) const
{ return Tensor2(cx[d], cy[d]); }



void Stencil::Print() const
{
    std::cout << "        d,\t        w,\t      phi,\t       cx,\t       cy,\n";
    for(size_t d=0; d<nDir; d++)
        std::cout << Format(d) << ",\t" << Format(W(d)) << ",\t" << Format(Phi(d)) << ",\t" << Format(Cx(d)) << ",\t" << Format(Cy(d)) << ",\n";
    std::cout << std::endl;
}
void Stencil::WriteToCsv() const
{
    std::ofstream fileOut(OUTPUTDIR + name + ".csv");
    fileOut << "#d,w,phi,cx,cy,cz\n";
    for(size_t d=0; d<nDir; d++)
        fileOut << Format(d) << "," << Format(W(d)) << "," << Format(Phi(d)) << "," << Format(Cx(d)) << "," << Format(Cy(d)) << ",0\n";
}
// -------------------------------------------------