#ifndef __RANDOM_H__
#define __RANDOM_H__

class BoostRandom 
{
	private:

		int SEED;

	public:

		BoostRandom(int aSeed);
		~BoostRandom(){}

		virtual double NextGaussian();
		double NormalCdf(double);
                int GetSeed() const { return SEED; }
};

#endif
