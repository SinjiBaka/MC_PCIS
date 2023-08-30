#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

const double pi = 3.1415926535;
const double c = 299792458.;			// speed of light [m/s]
const double k = 1.380649e-23;			// ���������� ��������� [��/�]
const double mp = 1.67262192369e-27;	// ����� ������� [��]

namespace lya
{
	const double wl = 121.567e-9;		// wave length [m]
	const double sigmaD = 5.9e-18;		// ������������ ������� ���������� ��� ������� ������ [m^2]
	const double sigmaN = 2.6e-23;		// ������������ ������� ���������� ��� ������� ������ [m^2]
}

namespace medium
{
	const double T = 2e4;				// ����������� ����� [K]
	const double n = 2e22;				// ��������� [m^{-3}]
}

#endif
