#ifndef RESOURCES_H_
#define RESOURCES_H_

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>//<boost/random/random_number_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/cstdint.hpp>

//extern unsigned long long accept; //test
//extern unsigned long long flip; //test
//extern unsigned long long stay; //test

//------------------------------------------------------------------------------
//
//	Template for functor for integer types (range is inclusive)
//
//------------------------------------------------------------------------------
class RandNumberBase
{
protected:
	static unsigned int s_unSeedShift;
};

template <typename IntType>
class BoostRandomInt : public RandNumberBase
{
private:
	boost::variate_generator< boost::mt19937_64&, boost::uniform_int<IntType> > *m_pvgRNG;
	boost::mt19937_64 *m_pBoostGenerator;

public:
	BoostRandomInt(IntType itMin, IntType itMax, unsigned int unSeed);
	~BoostRandomInt();
	IntType operator()() const { return (*m_pvgRNG)(); }
	void Seed() { (*m_pBoostGenerator).seed(static_cast<unsigned int>(std::time(0))); };
	void Seed(unsigned int unSeed) { (*m_pBoostGenerator).seed(unSeed); };
};


template <typename IntType>
BoostRandomInt<IntType>::BoostRandomInt(IntType itMin, IntType itMax, unsigned int unSeed=static_cast<unsigned int>(std::time(0)))
{
	m_pBoostGenerator = new boost::mt19937_64(unSeed + s_unSeedShift);
	boost::uniform_int<IntType> Uniform(itMin, itMax);
	m_pvgRNG = new boost::variate_generator< boost::mt19937_64&, boost::uniform_int<IntType> >(*m_pBoostGenerator, Uniform);

	// Increment seed for each new object created to prevent overlap
	s_unSeedShift++;
}
template <typename IntType>
BoostRandomInt<IntType>::~BoostRandomInt()
{
	delete m_pvgRNG;
	delete m_pBoostGenerator;
}

//------------------------------------------------------------------------------
//
//	Template for functor for real types (range is exclusive of max)
//
//------------------------------------------------------------------------------
template <typename RealType>
class BoostRandomReal : public RandNumberBase
{
private:
	boost::variate_generator< boost::mt19937&, boost::uniform_real<RealType> > *m_pvgRNG;
	boost::mt19937 *m_pBoostGenerator;

public:
	BoostRandomReal(RealType rtMin, RealType rtMax, unsigned int unSeed);
	~BoostRandomReal();
	RealType operator()() const { return (*m_pvgRNG)(); }
	void Seed() { (*m_pBoostGenerator).seed(static_cast<unsigned int>(std::time(0))); };
	void Seed(unsigned int unSeed) { (*m_pBoostGenerator).seed(unSeed); };
};


template <typename RealType>
BoostRandomReal<RealType>::BoostRandomReal(RealType rtMin, RealType rtMax, unsigned int unSeed=static_cast<unsigned int>(std::time(0)))
{
	m_pBoostGenerator = new boost::mt19937(unSeed + s_unSeedShift);
	boost::uniform_real<RealType> RealUniform(rtMin, rtMax);
	m_pvgRNG = new boost::variate_generator< boost::mt19937&, boost::uniform_real<RealType> >(*m_pBoostGenerator , RealUniform);

	// Increment seed for each new object created to prevent overlap
	s_unSeedShift++;
}
template <typename RealType>
BoostRandomReal<RealType>::~BoostRandomReal()
{
	delete m_pvgRNG;
	delete m_pBoostGenerator;
}


//------------------------------------------------------------------------------
//
//	Remove function for vectors
//
//------------------------------------------------------------------------------
//template <typename T>
//void remove(std::vector<T>& vec, size_t pos)
//{
//    std::vector<T>::iterator it = vec.begin();
//    std::advance(it, pos);
//    vec.erase(it);
//}

//------------------------------------------------------------------------------
//
//	Memory clearning function for vectors of pointers
//
//------------------------------------------------------------------------------
// Functor for deleting pointers in vector.
template<class T> class DeleteVector
{
    public:
    // Overloaded () operator.
    // This will be called by for_each() function.
    bool operator()(T x) const
    {
        // Delete pointer.
        delete x;
        return true;
    }
};

//------------------------------------------------------------------------------
//
//	Boltzmann function and simple random functions
//
//------------------------------------------------------------------------------
double Boltzmann(double dCE, double dT);
// long stdRandomInt(long nLowBound, long nUpBound);
// double stdRandomDouble(double dLowBound, double dUpBound);
// bool stdRandomBool();


//------------------------------------------------------------------------------
//
//	String splitting function
//
//------------------------------------------------------------------------------
// Courtesy of Evan Teran
std::vector<std::string> &StringSplit(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> StringSplit(const std::string &s, char delim);

#endif /* RESOURCES_H_ */
