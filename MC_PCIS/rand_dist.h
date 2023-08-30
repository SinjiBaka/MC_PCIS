#pragma once
#ifndef RAND_DIST_H
#define RAND_DIST_H

#include <random>
#include <iostream>

#include "constants.h"


class Rand_dist
{
	std::mt19937 mt{ std::random_device{}() };
	std::uniform_real_distribution<> dis{ 0., 1. };

public:

	Rand_dist() {}
	~Rand_dist() {}

	double explorenz_zheng(double x, double a, double u0 = 3.);
	double explorenz_seon(double x, double a);
	double normal(double mean = 0., double sigma = std::sqrt(0.5));
	void isotropic_dir(double& nx, double& ny, double& nz);
	void dipole_dir(double& nx, double& ny, double& nz);
	double uniform01() { return dis(mt); }
	double uniform(double a, double b);
	

private:
	double p(double beta);
	double Q_star(double beta, double x);
	double Q(double beta, double x, double a);
};

#endif
