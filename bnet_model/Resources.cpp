#include "Resources.h"

using namespace std;

//unsigned long long accept = 0ULL; //test
//unsigned long long flip = 0ULL; //test
//unsigned long long stay = 0ULL; //test

// set initial value of static method of RandNumberBase
unsigned int RandNumberBase::s_unSeedShift = 0U;

//------------------------------------------------------------------------------
//
//	Splits strings
//
//------------------------------------------------------------------------------
vector<string> &StringSplit(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> StringSplit(const string &s, char delim)
{
    vector<string> elems;
    StringSplit(s, delim, elems);
    return elems;
}

// //------------------------------------------------------------------------------
// //
// //  Returns the Boltzmann probability
// //
// //------------------------------------------------------------------------------
double Boltzmann(double dE, double dT)
{
    return exp( -dE / dT );
}


// //------------------------------------------------------------------------------
// //
// //  Returns a random integer within the range specified.
// //
// //------------------------------------------------------------------------------
// long stdRandomInt(const long nLowBound, const long nUpBound)
// {
//     srand(time(0));
//     return static_cast<long>(rand() / ((RAND_MAX + 1.0) / nUpBound)) + nLowBound;
// }


// //------------------------------------------------------------------------------
// //
// //  Returns a random floating point value within the range specified.
// //
// //------------------------------------------------------------------------------
// double stdRandomDouble(const double dLowBound, const double dUpBound)
// {
//     srand(time(0));
//     double dUniform = static_cast<double>(rand()) / (RAND_MAX + 1.0);
//     return dUniform*(dUpBound - dLowBound) + dLowBound;
// }


// //------------------------------------------------------------------------------
// //
// //  Returns a random boolean value.
// //
// //------------------------------------------------------------------------------
// bool stdRandomBool()
// {
//     if (stdRandomInt(0,1)) { return true; }
//     else { return false; }
// }
