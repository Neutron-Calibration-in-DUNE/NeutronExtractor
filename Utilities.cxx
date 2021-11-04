#include "Utilities.h"

namespace neutron
{
    double EuclideanDistance(std::vector<double> p1, std::vector<double> p2)
    {
        double dist = 0;
        for (size_t i = 0; i < p1.size(); i++) {
            dist += (p1[i] - p2[i])*(p1[i] - p2[i]);
        }
        return sqrt(dist);
    }

    double EuclideanDistance(double x1, double y1, double z1, 
                             double x2, double y2, double z2)
    {
        return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
    }
}