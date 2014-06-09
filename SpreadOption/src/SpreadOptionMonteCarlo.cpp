#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "BoostRandom.h"
#include "inverseCumulatives.h" 
#include "sobol_gold.h"
#include "sobol_primitives.h" 
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

using namespace std;
using namespace boost::accumulators;

double S_t = 110.0;
double L_t = 100.0;
double K = 0.0;
double T = 1.0;
double SigmaS = 0.1;
double SigmaL = 0.15;
double r = 0.05;
double repoS = 0.03;
double repoL = 0.02;
double rho = 0.0;

void CalculateMeanAndVariance(vector<double>::const_iterator& itBegin,vector<double>::const_iterator itEnd, size_t count, double& mean, double& variance)
{
	double sum = 0.0;
	double sumSquared = 0.0;
	for (vector<double>::const_iterator it = itBegin;it!=itEnd;it++)
	{
		sum+=*it;
		sumSquared+=*it**it;
	}
	
	mean = sum/count;
	variance = sumSquared/count;
}

void GenerateGaussianVariates(vector<double>& variates, size_t count, bool is_sobol, size_t dimensionCount = 1)
{
	double mean = 0.0;
	double variance = 0.0;

	variates.resize(count*dimensionCount);		
	if(is_sobol)
	{	
		//cout<<"Generating "<<count<<" Sobol Normal Variates for "<<dimensionCount<<" dimensions!"<<endl;
		vector<unsigned int> directions_vec(dimensionCount);
		initSobolDirectionVectors(dimensionCount,&directions_vec[0]);
		sobolCPU(count,  dimensionCount,  &directions_vec[0], &variates[0]);
				
		inverseCumulative<float> inv;
		transform(variates.begin(), variates.end(),variates.begin(),inv);
		
#ifdef PRINT_MEAN_AND_VARIANCE
		CalculateMeanAndVariance(variates.begin()+1,variates.begin()+count,count-1,mean,variance);
		cout<<"Dim 1 Sobol Mean     : "<<mean<<endl;
		cout<<"Dim 1 Sobol Variance : "<<variance<<endl;

		if (dimensionCount==1)
			return;
		CalculateMeanAndVariance(variates.begin()+count+1,variates.end(),count-1,mean,variance);
		cout<<"Dim 2 Sobol Mean     : "<<mean<<endl;
		cout<<"Dim 2 Sobol Variance : "<<variance<<endl;
#endif
	}
	else
	{
		//cout<<"Generating "<<count<<" Pseudo Random Normal Variates!"<<endl;
		BoostRandom random(666);
		for(size_t i = 0;i<count*dimensionCount;i++)
		{
			variates[i]=random.NextGaussian();
		}
	
#ifdef PRINT_MEAN_AND_VARIANCE
		CalculateMeanAndVariance(variates.begin()+count+1,variates.end(),count-1,mean,variance);
		cout<<"Mean     : "<<mean<<endl;
		cout<<"Variance : "<<variance<<endl;
#endif
	}
}

double SolveForSpreadOptionCallValue(bool is_sobol, size_t numberOfSimulations)
{
	double payoffSum = 0.0;
	double payoffSumSquared = 0.0;

	vector<double> gaussianVariates;
	GenerateGaussianVariates(gaussianVariates, numberOfSimulations, is_sobol, 2);

	double F_S = S_t*exp((r-repoS)*T);
	double F_L = L_t*exp((r-repoL)*T);

	for (size_t i = 1; i <numberOfSimulations;i++)
	{
		double X_S = gaussianVariates.at(i);
		double X_L = rho*X_S+sqrt(1-rho)*gaussianVariates.at(numberOfSimulations-1+i);
		double S_T = F_S*exp((-0.5*SigmaS*SigmaS)*T+SigmaS*X_S*sqrt(T));
		double L_T = F_L*exp((-0.5*SigmaL*SigmaL)*T+SigmaL*X_L*sqrt(T));
		double pathwisePayoff = ((S_T-L_T)-K>0.0?(S_T-L_T)-K:0.0);
		payoffSum+=pathwisePayoff;
		payoffSumSquared+=pathwisePayoff*pathwisePayoff;
	}

	double optionValue = exp(-r*T)*(payoffSum/(numberOfSimulations-1));

	return optionValue;
}

double SolveForCallValue(bool is_sobol, size_t numberOfSimulations)
{
	double payoffSum = 0.0;
	double payoffSumSquared = 0.0;

	vector<double> gaussianVariates;
	GenerateGaussianVariates(gaussianVariates, numberOfSimulations, is_sobol);

	for (size_t i = 1; i <numberOfSimulations;i++)
	{
		double S_T = S_t*exp((r-0.5*SigmaS*SigmaS)*T+SigmaS*gaussianVariates.at(i)*sqrt(T));
		double pathwisePayoff = (S_T-K>0.0?S_T-K:0.0);
		payoffSum+=pathwisePayoff;
		payoffSumSquared+=pathwisePayoff*pathwisePayoff;
	}

	double optionValue = exp(-r*T)*(payoffSum/(numberOfSimulations-1));

	return optionValue;
}

double SolveForSpreadOptionCallValueWithKirkFormula()
{
	BoostRandom random(666);
	double F_S = S_t*exp((r-repoS)*T);
	double F_L = L_t*exp((r-repoL)*T);

	double sigma = sqrt(SigmaS*SigmaS+SigmaL*SigmaL*(F_L/(F_L+K))*(F_L/(F_L+K)-2*rho*SigmaS*SigmaL*(F_L/(F_L+K))));
	double S_A = (F_S/(F_L+K));
	double d1 = (log(S_A)+(0.5*sigma*sigma*T))/(sigma*sqrt(T));
	double d2 = d1 - sigma*sqrt(T);
	double value = exp(-r*T)*(F_S*random.NormalCdf(d1)-(F_L+K)*random.NormalCdf(d2));
	return value;
}

double SolveForCallValueWithBlackScholes()
{
	double d1 = (log(S_t/K)+(r+0.5*SigmaS*SigmaS)*T)/(SigmaS*sqrt(T));
	double d2 = d1 - SigmaS*sqrt(T);

	BoostRandom random(666);

	return S_t*random.NormalCdf(d1)-random.NormalCdf(d2)*K*exp(-r*T);
}

int main(int argc, char* argv[])
{
	size_t numberOfSimulations = boost::lexical_cast<size_t>(argv[1]);
	bool isSpreadOption = boost::lexical_cast<bool>(argv[2]);
	rho = boost::lexical_cast<double>(argv[3]);

	if (!isSpreadOption)
	{
		cout<<endl<<"Call Value with BS    : "<<SolveForCallValueWithBlackScholes()<<endl<<endl;
		for (size_t i = 10000;i<numberOfSimulations;i+=10000)
		{
			cout<<i<<" paths Randoms: "<<SolveForCallValue(false,i)<<endl;
		}
		cout<<endl;
		for (size_t i = 10000;i<numberOfSimulations;i+=10000)
		{
			cout<<i<<" paths Sobol  : "<<SolveForCallValue(true, i)<<endl;
		}

	}
	else
	{
		cout<<endl<<"Call Value with Kirk Formula    : "<<SolveForSpreadOptionCallValueWithKirkFormula()<<endl<<endl;
		for (size_t i = 10000;i<numberOfSimulations;i+=10000)
		{
			cout<<i<<" paths Randoms: "<<SolveForSpreadOptionCallValue(false,i)<<endl;
		}
		cout<<endl;
		for (size_t i = 10000;i<numberOfSimulations;i+=10000)
		{
			cout<<i<<" paths Sobol  : "<<SolveForSpreadOptionCallValue(true, i)<<endl;
		}
	}
	return 1;
}
