#pragma once
#ifndef MCRT_H
#define MCRT_H


#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <filesystem>

#include "constants.h"
#include "mc_packet.h"



class Mcrt
{
	double R;	// радиус области расчета
	std::vector<MC_packet> packs;
	Rand_dist dist_class;

public:

	Mcrt() { R = 0; }
	Mcrt(double R, size_t Npacks);

	void execute_mcrt(int save_interval);
	std::vector<MC_packet> get_packs() { return packs; }
	
private:

	void propagate();
	double calc_distance2event(double tau, MC_packet p);
	size_t active_packs();
	void PIC_scatter(MC_packet& p, double tempr);
	void PIC_scatter_seon(MC_packet& p, double tempr);


	double voigt(double x2, double Temp, double vel2);
	double voigt_z(double x2);
	double voigt_q(double x2);
	double voigt_a(double T) { return 4.70574e-2 / sqrt(T); } // параметр функции Фойгта в зависимости от температуры в [К]

	void save_results();
	void save_avgs(std::string path);
	void save_spec(std::string path);
	void temp_save();
	
};

#endif