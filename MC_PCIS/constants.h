#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

const double pi = 3.1415926535;
const double c = 299792458.;			// speed of light [m/s]
const double k = 1.380649e-23;			// постоянная Больцмана [Дж/К]
const double mp = 1.67262192369e-27;	// масса протона [кг]

namespace lya
{
	const double wl = 121.567e-9;		// wave length [m]
	const double sigmaD = 5.9e-18;		// доплеровское сечение поглощения для функции Фойгта [m^2]
	const double sigmaN = 2.6e-23;		// естественное сечение поглощения для функции Фойгта [m^2]
}

namespace medium
{
	const double T = 2e4;				// температура среды [K]
	const double n = 2e22;				// плотность [m^{-3}]
}

#endif
