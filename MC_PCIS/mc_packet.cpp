#include "mc_packet.h"

MC_packet::MC_packet()
{
	// начальное положение
	x = 0.;
	y = 0.;
	z = 0.;

	// направление распространения
	dist_class.isotropic_dir(nx, ny, nz);

	// частота по умолчанию Lya или скорость 0
	v = 0.;

	is_active = true;
	scat_num = 0;
	distance = 0.;
}

MC_packet::MC_packet(double velocity)
{
	// начальное положение
	x = 0.;
	y = 0.;
	z = 0.;

	// направление распространения
	dist_class.isotropic_dir(nx, ny, nz);

	// частота по умолчанию Lya
	v = velocity;

	is_active = true;
	scat_num = 0;
	distance = 0.;
}