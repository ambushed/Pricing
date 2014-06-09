#include <cmath>
#include <complex>
#include <iostream>
#include "BoostRandom.h" 

using namespace std;

double S_t = 110.0;
double L_t = 100.0;
double K = 100.0;
double T = 1.0;
double SigmaS = 0.1;
double SigmaL = 0.15;
double r = 0.05;
double repoS = 0.03;
double repoL = 0.02;
double rho = 0.0;

double SolveForCallValueWithBlackScholes()
{
	double d1 = (log(S_t/K)+(r+0.5*SigmaS*SigmaS)*T)/(SigmaS*sqrt(T));
	double d2 = d1 - SigmaS*sqrt(T);

	BoostRandom random(666);

	return S_t*random.NormalCdf(d1)-random.NormalCdf(d2)*K*exp(-r*T);
}


class COS
{

	public:

		complex<double> char_func;
		vector<double> psi_k;
		vector<double> chi_k;
		vector<double> omega_k;
		double a;
		double b;

	public:
		COS(size_t N)
		{
			psi_k.assign(N,0.0);
			chi_k.assign(N,0.0);
			omega_k.assign(N,0.0);
			char_func.resize(N);
			double PI = acos(-1);
			for(size_t i = 0; i<N;i++)
			{
				omega_k.at(i) = i*PI/(b-a);
			}
		}

		double CalculateCharacteristicFunction()
		{
			std::complex<double> i(0,1);
			double mu = (r - 0.5*SigmaS*SigmaS)*T;
			for(size_t k = 0; k<N; k++)
			{
				char_func.at(k) = exp(i*omega_k.at(k)*T*mu-0.5*SigmaS*SigmaS*omega_k.at(k)*omega_k.at(k)*T);
			}

		}

		double CalculateV(double c, double d)
		{

			for (size_t k = 0;k<numTerms;k++)
			{
			
				chi_k.at(k) = 1.0/(1+pow(omega_k.at(k),2))*(cos(omega_k.at(k)*(d - a))*exp(x_2)-cos(omega_k.at(k)*(c - a)*exp(x_1))+
					omega_k.at(k)*sin(omega_k.at(k)*(d - a))*exp(x_2)-omega_k.at(k)*sin(omega_k.at(k)*(c - a))*exp(x_1));
			
				psi_k.at(k) = (sin(omega_k.at(k)*(d-a))-sin(omega_k.at(k)*(c-a)))*(b-a)/(k*PI);
			}
			
			psi_k.at(0) = d - c;


		}

		double SolveForCallValueWithCos()
		{
			double a = -5.0;
			double b = 5.0;
			double sum = 0.0;

			CalculateCharacteristicFunction();

			// c and d for Call option are 0 and b
			// ----------- Put option are a and 0
			CalculateV(0.0, b);

			for (size_t k = 0;k<N;k++)
			{
				complex<double> secondMultiplier = exp(i*omega_k.at(k)*(-a));

				double A_k = 2.0/(b-a)*real(char_func.at(k)*secondMultiplier);
				double V_k = 2.0/(b-a)*K*(chi_k.at(k)-psi.at(k));
				sum += A_k*V_k;

				if (k == 0)
					sum *= 0.5;
			}

			return exp(-r*T) * (b-a)/2.0 * sum;
		}
}

int main()
{
	COS cosMethod(10);

	cout<<endl<<"Call Value with BS    : "<<SolveForCallValueWithBlackScholes()<<endl<<endl;
	cout<<endl<<"Call Value with COS   : "<<cosMethod.SolveForCallValue()<<endl<<endl;
	
	return 1;
}
