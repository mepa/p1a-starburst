#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

using namespace std;

// assume that V_rot(R)=constant

// R^2*Omega = J = R*V_rot

double R_J(const double J, const double V) // assume V=const
{
  return J / V;
}

double dJ_dR(const double V)
{
  return V;
}

double r_J(const double J, const double V, const double theta)
{
  return R_J(J, V) / sin(theta);
}

// 1 / abs(diff(R V(R), R)) R=R_J  = 1 / V

double Dirac_factor(const double V)
{
  return 1.0 / V;
}

double rho(const double r, const double rho_0, const double R_e, const double n)
{
  const double p = 1.0 - 0.6097 / n + 0.05563 / (n * n);  // semi-empirical
  const double b = 2.0 * n - 1.0 / 3.0 + 0.009876 / n;
  const double x = r / R_e;
  return rho_0 * pow(x, - p) * exp(- b * pow(x, 1.0 / n));
}

double dN_dm(double m)
{
  const double Msun = 2.0e33;                   // g
  const double percentage_above_Msun = 0.36;    // from Kroupa, MNRAS (2001)
  const double slope = - 2.3;                   // "
  const double A = percentage_above_Msun * (- slope) * pow(Msun, - slope - 1.0);
  return A * pow(m, slope);
}

double mAverage()
{
  const double Msun = 2.0e33;                   // g
  return 0.2 * Msun;
}

double m_of_t(double t)
{
  const double Msun = 2.0e33;                   // g
  const double power = 3.0;
  const double year = 3.0e7;                    // sec
  const double t_sun = 1.0e10 * year;           // Sun's lifetime is 10 bln yrs (current age is 5 bln yrs)
  return Msun * pow(t / t_sun, - 1.0 / power);  // t_age ~ m^-3 (t_age = L/M)
}

double dm_dt(double t)
{
  const double Msun = 2.0e33;                  
  const double power = 3.0;
  const double year = 3.0e7;                    
  const double t_sun = 1.0e10 * year;
  return Msun * (1.0 / power) * pow(t / t_sun, - 1.0 / power - 1.0) / t_sun;


double delta_m(double m)
{
  const double Msun = 2.0e33;
  const double m_white_dwarf = Msun + (m - 6.0 * Msun) * 0.1;   // ???
  return m - m_white_dwarf;
}

double integrand(const double theta, const double J, const double V, const double t,
		 const double rho_0, const double R_e, const double n)
{
  const double _r_J = r_J(J, V, theta);
  const double m = m_of_t(t);
  const double _rho = rho(_r_J, rho_0, R_e, n);

  return 2.0 * M_PI * _r_J * _r_J * (_rho / mAverage()) * dN_dm(m) * dm_dt(t) * delta_m(m) * Dirac_factor(V);
}

double calc_rho_0(const double R_e, const double n, const double M_galaxy)
{
  const double kpc = 3.0e21;                    // cm
  const double r_min = 0.01 * kpc;
  const double r_max = 10.0 * kpc;
  const size_t N = 1000;
  vector<double> radii(N);

  for(size_t i = 0; i < N; i ++)
    {
      radii[i] = r_min * exp(double(i) * log(r_max / r_min) / double(N - 1));
    } 

  double sum = 0.0;

  for(size_t i = 0; i < (N - 1); i ++)
    {
      sum += 0.5 * (radii[i + 1] - radii[i]) *
	(4.0 * M_PI * radii[i] * radii[i] * rho(radii[i], 1.0, R_e, n) +
	 4.0 * M_PI * radii[i + 1] * radii[i + 1] * rho(radii[i + 1], 1.0, R_e, n));
    } // Integral

  return M_galaxy / sum;  // M_galaxy = rho_0 * Int[4*Pi*r^2*f(r)*dr]
}

double dsqM_dJ_dt(const double M_galaxy, const double rho_0, const double R_e, const double n, const double J, const double V, const double t)
{
  const size_t N = 1000;
  const double delta_theta = M_PI / double(N + 1);
  double sum = 0.0;

  for(size_t i = 0; i < N; i ++)
    {
      const double theta = double(i + 1) * delta_theta;

      const double value = integrand(theta, J, V, t, rho_0, R_e, n);
      sum += delta_theta * value;
    }

  return sum;
}

double MassEnclosed(const double rho_0,
		    const double n,
		    const double R_e,
		    const double r)
{
  static const size_t c = 1000;
  static double * radii = new double [c];
  static double * dens_r_sq = new double [c];
  
  const double r_min = 1.0e-3 * r;
  const double factor = - log(1.0e-3) / double(c - 1);

  for(size_t i = 0; i < c; i ++)
    {
      radii[i] = r_min * exp(double(i) * factor);

      dens_r_sq[i] = rho(radii[i], rho_0, R_e, n) * radii[i] * radii[i];
    }

  double sum = 0.0;

  for(size_t i = 0; i < (c - 1); i ++)
    {
      sum += 0.5 * (radii[i + 1] - radii[i]) * (dens_r_sq[i] + dens_r_sq[i + 1]);
    }

  return 4.0 * M_PI * sum;
}

double V_circular(const double rho_0,
		  const double n,
		  const double R_e,
		  const double r)
{
  const double G = 6.67e-8;

  return sqrt(G * MassEnclosed(rho_0, n, R_e, r) / r);
}


double R_circularization(const double J, 
			 const double rho_0,
			 const double n,
			 const double R_e)
{
  const double kpc = 3.0e21;
  double r_min = 1.0e-6 * kpc;
  double r_max = 1.0e3 * kpc;
  const double r_tiny = 1.0e-6 * kpc;

  while(r_max - r_min > r_tiny)
    {
      const double r_middle = 0.5 * (r_min + r_max);

      const double J_middle = r_middle * V_circular(rho_0, n, R_e, r_middle);

      if(J < J_middle) r_max = r_middle;
      else if(J > r_middle) r_min = r_middle;
      else return r_middle;
    }
  
  return 0.5 * (r_min + r_max);
}
			 
double dJ_dR_circ(const double R_circ, 
		  const double rho_0,
		  const double n,
		  const double R_e)
{
  const double V_circ = V_circular(rho_0, n, R_e, R_circ);
  const double M_encl = MassEnclosed(rho_0, n, R_e, R_circ);
  const double dM_encl_dR = 4.0 * M_PI * R_circ * R_circ * rho(R_circ, rho_0, R_e, n);
  const double dV_circ_dR = -0.5 * V_circ / R_circ + 0.5 * (V_circ / M_encl) * dM_encl_dR;

  return V_circ + R_circ * dV_circ_dR;
}

int main()
{
  const double M_sun = 2.0e33;                      // g
  const double M_galaxy = 1.0e11 * M_sun;
  const double n = 4.0;
  const double kpc = 3.0e21;                        // cm
  const double R_e = 2.0 * kpc;                     // why 2 kpc?

  const double rho_0 = calc_rho_0(R_e, n, M_galaxy);

  const size_t N = 100;

  vector<double> J(N);
  const double J_min = 0.01 * kpc * 5.0e6;
  const double J_max = 100.0 * kpc * 5.0e6;

  for(size_t i = 0; i < N; i ++)
    {
      J[i] = J_min * exp(double(i) * log(J_max / J_min) / double(N - 1));
    }

  const double V = 5.0e6;                           // cm
  const double year = 3.0e7;                        // sec
  const double t = 1.0e8 * year;                    // Approximate age of typical AGB star
       
  ofstream f1_out("d2MdRdt_R.txt");
  for(size_t i = 0; i < N; i ++)
    {
      const double R_circ = R_circularization(J[i], rho_0, n, R_e);

      f1_out << R_circ / kpc << " "
	     << kpc * dJ_dR_circ(R_circ, rho_0, n, R_e) * t * dsqM_dJ_dt(M_galaxy, rho_0, R_e, n, J[i], V, t) / M_sun 
	     << endl;
    }
  
//  const double kpc = 3.0e21;                      // cm
  const double r_min = 0.01 * kpc;
  const double r_max = 10.0 * kpc;
//  const size_t N = 1000;
  vector<double> radii(N);
  
  ofstream f2_out("rho_r.txt");
  for(size_t i = 0; i < N; i ++)
  {
    radii[i] = r_min * exp(double(i) * log(r_max / r_min) / double(N - 1));
    f2_out << radii[i] << " " << rho(radii[i], rho_0, R_e, n) << endl;
  }  

  cout << rho_0 << endl;
  
  //  system("PAUSE");
  return EXIT_SUCCESS;
  //return 0;
}

