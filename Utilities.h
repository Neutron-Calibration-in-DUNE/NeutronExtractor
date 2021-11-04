/**
 * @file    Utilities.h
 * @brief   A class for various Utility functions
 * @ingroup Utilities
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/
#pragma once
#include <vector>
#include <string>
#include <cmath>

namespace neutron
{
    double EuclideanDistance(std::vector<double> p1, std::vector<double> p2);
    double EuclideanDistance(double x1, double y1, double z1, 
                             double x2, double y2, double z2);
}