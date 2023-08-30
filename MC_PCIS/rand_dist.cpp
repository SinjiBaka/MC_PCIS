#include "rand_dist.h"

double Rand_dist::explorenz_zheng(double x, double a, double u0)
{
	double ux = std::abs(x);
	int sx = (x >= 0 ? 1 : -1);
	x = ux;

	double theta0 = std::atan2(u0 - x, a);
	double p = (theta0 + pi / 2.) / (theta0 * (1. - std::exp(-u0 * u0)) + pi / 2. * (1. + std::exp(-u0 * u0)));

	double R1, R2, theta, u;
	while (true)
	{
		R1 = dis(mt);
		R2 = dis(mt);

		if (R1 <= p)
		{
			theta = dis(mt) * (theta0 + pi / 2.) - pi / 2.;
			u = a * std::tan(theta) + x;
			if (R2 < std::exp(-u * u))
				return sx * u;
		}
		else
		{
			theta = dis(mt) * (pi / 2. - theta0) + theta0;
			u = a * std::tan(theta) + x;
			if (R2 < std::exp(-u * u) / std::exp(-u0 * u0))
				return sx * u;
		}
	}
}

double Rand_dist::explorenz_seon(double x, double a)
{
   
    if (x < 1.)
    {
        return explorenz_zheng(x, a, 0.);
    }

    double ux = std::abs(x);
    int sx = (x >= 0 ? 1 : -1);
    x = ux;

    double xc = 2.41421; // 1 + sqrt(2)
    double b0 = std::exp(-x * x / 2.);
    double b1 = b0 + std::sqrt(2. * a / pi * (1. - b0) * b0 * x);

    double h0 = b0 / 2. / a;
    double h1 = Q_star(b1, x);
    double h2 = 0.3861 / (x * x - 1.373);


    double b, C;
    while (true)
    {
        if (h2 >= 2. * h0)
        {
            b = dis(mt);
            C = h2;
        }
        else if (h2 >= h0 && h2 < 2. * h0)
        {
            double S0 = b0 * h0;
            double S1 = (1. - b0) * h2;

            if (dis(mt) <= (S0 / (S0 + S1)))
            {
                b = b0 * std::sqrt(dis(mt));
            }
            else
            {
                b = b0 + (1. - b0) * dis(mt);
            }
            C = (b <= b0) ? b / a : h2;
        }
        else
        {
            double S0 = b0 * h0;
            double S1 = (b1 - b0) * h0;
            double S2 = (x < xc) ? (1. - b1) * h1 : (1. - b0) * std::max(h1, h2);
            double St = S0 + S1 + S2;

            double xi1 = dis(mt);
            if (xi1 <= (S0 / St))
            {
                b = b0 * std::sqrt(dis(mt));
            }
            else if (xi1 > (S0 / St) && xi1 <= (1. - S2 / St))
            {
                b = b0 + (b1 - b0) * dis(mt);
            }
            else
            {
                b = b1 + (1. - b1) * dis(mt);
            }

            if (b <= b0)
                C = b / a;
            else if (b > b0 && b <= b1)
                C = h0;
            else if (b > b1 && x <= xc)
                C = h1;
            else
                C = std::max(h1, h2);
        }

        if (dis(mt) <= Q(b, x, a) / C)
        {
            return sx * (x + a * std::tan(a * pi / b * Q(b, x, a) * dis(mt) - std::atan2(p(b) + x, a)));
        }
    }
}

double Rand_dist::normal(double mean, double sigma)
{
	std::normal_distribution<double> norm_dist{ mean, sigma };
	return norm_dist(mt);
}

void Rand_dist::isotropic_dir(double& nx, double& ny, double& nz)
{
	double phi = 2. * pi * dis(mt);
	double theta = std::acos(2. * dis(mt) - 1.);

	nx = std::sin(theta) * std::cos(phi);
	ny = std::sin(theta) * std::sin(phi);
	nz = std::cos(theta);
}

void Rand_dist::dipole_dir(double& nx, double& ny, double& nz)
{
	double mu = dis(mt);
	mu = std::cbrt((4. * mu - 2.) + std::sqrt((4. * mu - 2.) * (4. * mu - 2.) + 1.)) - 1./std::cbrt((4. * mu - 2.) + std::sqrt((4. * mu - 2.) * (4. * mu - 2.) + 1.));

	double phi = 2. * pi * dis(mt);
	double theta = std::acos(mu);

	nx = std::sin(theta) * std::cos(phi);
	ny = std::sin(theta) * std::sin(phi);
	nz = std::cos(theta);
}

double Rand_dist::uniform(double a, double b)
{
	return b - dis(mt) * (b - a);
}

double Rand_dist::p(double beta)
{
	return std::sqrt(-2. * std::log(beta));
}

double Rand_dist::Q_star(double beta, double x)
{
	return 2. * beta / pi * p(beta) / (x * x - std::pow(p(beta), 2.));
}

double Rand_dist::Q(double beta, double x, double a)
{
	return beta / pi / a * (std::atan2(p(beta) - x, a) + std::atan2(p(beta) + x, a));
}