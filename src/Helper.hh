#ifndef __INCLUDE_GUARD_Helper_hh__
#define __INCLUDE_GUARD_Helper_hh__
#include <vector> // vector

// Return vector of Factors for n
inline std::vector<int> Factors(int n)
{
    std::vector<int> factors;
    for (int i = 1; i <= n; i++)
        if (n % i == 0)
            factors.push_back(i);
    return factors;
}

// determinate subgrid size and maximum usable nodes
// overwrites: bestN, bestX, bestY
inline void FindSubgridSize(int nx, int ny, int nodes, int &bestX, int &bestY, int &bestN)
{
    std::vector<int> widthFactors = Factors(nx);
    std::vector<int> heightFactors = Factors(ny);
    int totalArea = nx * ny;

    float bestU = 0; // circumference 2x+2y

    for (int j = 0; j < heightFactors.size(); j++)
    {
        for (int i = 0; i < widthFactors.size(); i++)
        {
            int area = widthFactors[i] * heightFactors[j];
            if (totalArea % area == 0 && totalArea / area <= nodes)
            {
                if (totalArea / area > bestN)
                {
                    bestX = widthFactors[i];
                    bestY = heightFactors[j];
                    bestN = totalArea / area;
                    bestU = bestX * 2 + bestY * 2;
                }
                else if (totalArea / area == bestN && widthFactors[i] * 2 + heightFactors[j] * 2 < bestU)
                {
                    bestX = widthFactors[i];
                    bestY = heightFactors[j];
                    bestN = totalArea / area;
                    bestU = bestX * 2 + bestY * 2;
                }
            }
        }
    }
}

#endif //__INCLUDE_GUARD_Helper_hh__