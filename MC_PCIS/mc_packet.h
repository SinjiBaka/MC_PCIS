#pragma once
#ifndef MC_PACKET_H
#define NC_PACKET_H

#include <random>
#include <iostream>

#include "constants.h"
#include "rand_dist.h"

class MC_packet
{
	Rand_dist dist_class;

public:
	double x, y, z;		// ���������� ������
	double nx, ny, nz;	// ����������� ���������������
	double v;			// ������������� �������� ����� [�/�], 
						// �� ������� ���������� �������, 
						//������������ ����������� ������� ������ ����� (v=0 ������ freq = f_{Lya})
	bool is_active;

	// ����� ����� ��������� �������������� ���������
	size_t scat_num; // ����� ������� ���������
	double distance; // ���������� ����������


	MC_packet();
	MC_packet(double velocity);


};

#endif

