#pragma warning(disable : 4996)
#include "mcrt.h"


Mcrt::Mcrt(double R, size_t Npacks)
{
	this->R = R;
	for (size_t i = 0; i < Npacks; i++)
		packs.emplace_back();
}

void Mcrt::execute_mcrt(int save_interval)
{
	size_t Nactive = active_packs();

	auto begin = std::chrono::steady_clock::now();
	int t = 1;

	while(Nactive > 0)
	{
		propagate();
		Nactive = active_packs();
		

		if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - begin).count() >= t * 60 * save_interval)
		{
			temp_save();
			std::cout << "\tMinutes passed:\t" << save_interval * t << '\n';
			std::cout << "\tPhoton remain:\t" << Nactive << "\n";
			t++;
		}

	}

	auto end = std::chrono::steady_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "Elapsed time:\t" << elapsed_time.count() << "\tms\n";

	save_results();
}

void Mcrt::propagate()
{

	for (size_t i = 0; i < packs.size(); i++)
	{
		if (packs[i].is_active)
		{
			// расстояние, которое должен пролететь
			double tau = -std::log(dist_class.uniform01());
			double l = calc_distance2event(tau, packs[i]);
			
			// новые координаты
			packs[i].x = packs[i].x + l * packs[i].nx;
			packs[i].y = packs[i].y + l * packs[i].ny;
			packs[i].z = packs[i].z + l * packs[i].nz;
			packs[i].distance += l;

			double r = std::sqrt(packs[i].x * packs[i].x + packs[i].y * packs[i].y + packs[i].z * packs[i].z);
			if (r > R)
			{
				packs[i].distance -= l - R;
				packs[i].is_active = false;
			}
			else
			{
				// если не вылетел, то рассеивается
				PIC_scatter(packs[i], medium::T);
				packs[i].scat_num++;
			}
		}	
	}
}

void Mcrt::PIC_scatter(MC_packet& p, double tempr)
{
	// новое направление
	double old_nx = p.nx;
	double old_ny = p.ny;
	double old_nz = p.nz;

	dist_class.isotropic_dir(p.nx, p.ny, p.nz);
	

	double mu = old_nx * p.nx + old_ny * p.ny + old_nz * p.nz;

	// новая доплеровская скорость (частота)
	double vth = std::sqrt(2. * k * tempr / mp);

	double u = dist_class.explorenz_seon(p.v / vth, voigt_a(tempr));
	double w = dist_class.normal();

	double new_x = p.v / vth - u + u * mu + w * std::sqrt(1. - mu * mu);
	p.v = vth * new_x;
}

void Mcrt::PIC_scatter_seon(MC_packet& p, double tempr)
{
	/* немного другая обработка рессеяния, результаты совпадают с PIC_scatter*/

	// новое направление
	double old_nx = p.nx;
	double old_ny = p.ny;
	double old_nz = p.nz;

	dist_class.dipole_dir(p.nx, p.ny, p.nz);
	double mu = old_nx * p.nx + old_ny * p.ny + old_nz * p.nz;

	// новая доплеровская скорость (частота)
	double vth = std::sqrt(2. * k * tempr / mp);
	double xi1 = dist_class.uniform01();
	double xi2 = dist_class.uniform01();

	double ux = dist_class.explorenz_seon(p.v / vth, voigt_a(tempr));
	double uy = std::sqrt(-std::log(xi1)) * std::cos(2. * pi * xi2);
	double uz = std::sqrt(-std::log(xi1)) * std::sin(2. * pi * xi2);

	double new_x = (p.v / vth) - ux + p.nx * ux + p.ny * uy + p.nz * uz;

	// эффект отдачи (если число столкновений меньше 1е10 он мал)
	bool recoil = false;
	if (recoil)
		new_x -= (1. - mu) * (2.536e-4 / std::sqrt(tempr / 1e4));

	p.v = vth * new_x;
}

size_t Mcrt::active_packs()
{
	size_t Nact = 0;
	for (size_t i = 0; i < packs.size(); i++)
	{
		if (packs[i].is_active)
			Nact++;
	}

	return Nact;
}

double Mcrt::calc_distance2event(double tau, MC_packet p)
{
	double x2 = mp * pow(p.v, 2) / 2. / k / medium::T;

	double sigma = voigt(x2, medium::T, pow(p.v, 2));

	return tau / medium::n / sigma;
}

double Mcrt::voigt(double x2, double Temp, double vel2)
{
	double res = lya::sigmaD * 1e2 / std::sqrt(Temp) * std::exp(-x2);
	double q = voigt_q(x2);

	if(q != 0)
		res += lya::sigmaN * 1e10 / vel2 * voigt_q(x2);

	return res;
}

double Mcrt::voigt_z(double x2)
{
	return (x2 - 0.855) / (x2 + 3.42);
}

double Mcrt::voigt_q(double x2)
{
	double z = voigt_z(x2);
	if (z <= 0)
		return 0.;
	return (21. + x2) / (1. + x2) * z * (0.1117 + z * (4.421 + z * (5.674 * z - 9.207)));
}

void Mcrt::save_results()
{
	std::string path1, path2;

	std::string pref = "slns/spec_";
	std::string suf = ".dat";
	std::string root = std::to_string(R);
	size_t dot = root.find(".");
	root.replace(dot, 1, "_");

	path1 = pref + root + suf;
	path2 = "slns/avgs.dat";

	save_spec(path1);
	save_avgs(path2);

	if(std::filesystem::exists("slns/temporary_avgs.dat"))
		std::remove("slns/temporary_avgs.dat");
}

void Mcrt::temp_save()
{
	std::string path1, path2;

	std::string pref = "slns/spec_";
	std::string suf = ".dat";
	std::string root = std::to_string(R);
	size_t dot = root.find(".");
	root.replace(dot, 1, "_");

	path1 = pref + root + suf;
	path2 = "slns/temporary_avgs.dat";

	save_spec(path1);
	save_avgs(path2);
}

void Mcrt::save_avgs(std::string path)
{

	double tau0, lt = 0, Nscat = 0;
	int in_act = 0;
	tau0 = voigt(0., medium::T, 0.) * R * medium::n;

	for (size_t i = 0; i < packs.size(); i++)
	{
		if (!packs[i].is_active)
		{
			lt += packs[i].distance / c;
			Nscat += packs[i].scat_num;
			in_act++;
		}
	}
	lt = lt / in_act;
	Nscat = Nscat / in_act;


	std::ofstream fout;
	fout.open(path, std::ofstream::app);
	if (!fout.is_open())
	{
		std::cout << "Cannot open file to write average information!" << std::endl;
	}
	else
	{
		fout << "temperature\tdensity\tradius\ttau0\tlifetime\tscater_number\tpackets_num\n";
		fout << medium::T << '\t' << medium::n << '\t' << R << '\t' << tau0 << '\t' << lt << '\t' << Nscat << '\t' << in_act << '\n';
	}
	fout.close();
}

void Mcrt::save_spec(std::string path)
{
	std::ofstream fout;
	fout.open(path, std::ofstream::trunc);
	if (!fout.is_open())
	{
		std::cout << "Cannot open file to write spectral results!" << std::endl;
	}
	else
	{
		fout << "Dopler_velocity\tScattering_num\tDistance\tLifetime\n";
		for (size_t i = 0; i < packs.size(); i++)
		{
			if (!packs[i].is_active)
			{
				fout << packs[i].v << '\t' << packs[i].scat_num << '\t' << packs[i].distance << '\t' << packs[i].distance / c << '\n';
			}
		}
	}

	fout.close();
}


