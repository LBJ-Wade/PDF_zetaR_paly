#ifndef INCLUDED_SRK32_hpp_
#define INCLUDED_SRK32_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

class SRKintegrater
{
protected:
  double t,t0,
    A0[3][3],A1[3][3],B0[3][3],B1[3][3],Alpha[3],Beta1[3],Beta2[3],C0[3],C1[3];
  vector<double> dW; // dW[alpha]
  vector< vector<double> > xx,xxi,aa[3],H0[3]; // xx[xp][I]
  vector< vector<double> > aIs,vIs; // aIs[I][alpha] only for Idim == noisedim
  vector< vector< vector<double> > > bb[3],Hk[3]; //bb[3][alpha][xp][I]
  int xpdim,Idim,noisedim;

public:
  SRKintegrater(){}
  SRKintegrater(vector< vector<double> > &XPi, double T0, int NoiseDim);
  void SRK2(double dt);
  void coeff(double dt, int step);
  void init_txp();
  void set_txp(double T, vector< vector<double> > &PSV);
  double e1(vector<double> &X, vector<double> &P);
  double return_H();
  double return_V();
  double return_t();
  double return_xp(int xp, int I);
  double return_e1();
  double vielbein(vector< vector<double> > &XP, int I, int alpha);
  double eIsigma(vector< vector<double> > &XP, int I);
  double eIs(vector< vector<double> > &XP, int I, int alpha);
  double eta_perp(vector< vector<double> > &XP); // only for 2-field model
  double return_etaperp();
  virtual double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L);
  virtual double DI(int xp, int I, vector< vector<double> > &psv);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double gIa(int xp, int I, int alpha, vector< vector<double> > &psv);
};

double Uniform();
double rand_normal(double mu, double sigma);

#endif
