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
	double x, y, z;		// координаты пакета
	double nx, ny, nz;	// направление распространения
	double v;			// допплеровская скорость атома [м/с], 
						// на котором происходит событие, 
						//относительно центральной частоты Лайман альфа (v=0 значит freq = f_{Lya})
	bool is_active;

	// может нужны небольшие статистические параметры
	size_t scat_num; // число событий рассеяния
	double distance; // пройденное расстояние


	MC_packet();
	MC_packet(double velocity);


};

#endif

