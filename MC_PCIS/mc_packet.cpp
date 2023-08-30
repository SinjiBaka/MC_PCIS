#include "mc_packet.h"

MC_packet::MC_packet()
{
	// ��������� ���������
	x = 0.;
	y = 0.;
	z = 0.;

	// ����������� ���������������
	dist_class.isotropic_dir(nx, ny, nz);

	// ������� �� ��������� Lya ��� �������� 0
	v = 0.;

	is_active = true;
	scat_num = 0;
	distance = 0.;
}

MC_packet::MC_packet(double velocity)
{
	// ��������� ���������
	x = 0.;
	y = 0.;
	z = 0.;

	// ����������� ���������������
	dist_class.isotropic_dir(nx, ny, nz);

	// ������� �� ��������� Lya
	v = velocity;

	is_active = true;
	scat_num = 0;
	distance = 0.;
}