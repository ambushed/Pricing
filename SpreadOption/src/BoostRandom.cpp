#include "BoostRandom.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include <boost/format.hpp>

BoostRandom::BoostRandom(int aSeed)
        :SEED(aSeed)
{
}

double BoostRandom::NextGaussian()
{
        using namespace boost::random;
        static mt19937 rng(SEED);
        boost::normal_distribution<> nd(0.0, 1.0);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
        return var_nor();
}

double BoostRandom::NormalCdf(double aValue)
{
        boost::math::normal norm(0.0,1.0);
        return boost::math::cdf(norm,aValue);
}