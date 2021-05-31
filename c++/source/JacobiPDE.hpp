#ifndef INCLUDED_JacobiPDE_hpp_
#define INCLUDED_JacobiPDE_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class JacobiPDE
{
protected:
  vector< vector< vector<double> > > *site, *hI;
  // site[xp][I][index] : phase-space value at each site
  // hI[xp][I][index] : step size
  vector< vector<int> > *siteNo; // siteNo[xp][I] : site number in xpI
  vector<double> *ff, *f_next; // function to be solved
  
  // ---------------------------------------------------------
  /*
  num_p[num][xp][I] : site number positively next to num in xpI direction
  num_pp[num][xpI][I][xpJ][J] : positively next to positively next in xpI and xpJ direction
  hp[num][xp][I] : site step positively next to num in xpI direction
  set in BoundaryCondition
  */
  vector< vector< vector<int> > > *num_p, *num_m;
  vector< vector< vector< vector< vector<int> > > > > *num_pp, *num_pm, *num_mm;
  vector< vector< vector<double> > > *hp, *hm;
  // ---------------------------------------------------------
  
  double maxstep,tol,rhoc,xpdim,Idim;
  // maxstep : max step number for Jacobi algorithm
  // tol : tolerance error for Jacobi algorithm
  // rhoc : energy density on the end of inflation hypersurface
  // xpdim : xpdim=1 -> slow-roll field-space, xpdim=2 -> full phase-space
  // Idim : number of inflatons
  int volume, funcNo;
  // volume : total number of sites
  // funcNo : number of functions to be solved
  vector<int> xpvol; // xpvol[0] : volume in field space, xpvol[1] : volume in momentum space
  vector<bool> Omega; // Omega[siteNo] : the site is in inflationary range or not

public:
  JacobiPDE(){}
  // Params[0] = maxstep, Params[1] = tol, Params[2] = funcNo, Params[3] = rhoc
  JacobiPDE(vector< vector< vector<double> > > &Site, vector<double> &Params);
  double PDE_1step(int num, int func); // execute 1 Jacobi step at site[num]
  void PDE_solve(int func); // solve func by Jacobi
  int Ind2No(vector< vector<int> > &index); // index set to site number
  int No2Ind(int num, int xp, int I); // site number to index of xpI
  double No2PSV(int num, int xp, int I); // site number to phase-space value of xpI
  int ceilXP(int xp, int I, vector< vector<double> > &psv); // ceil of phase-space value psv in xpI direction
  double Interpolation_f(vector< vector<double> > &psv, int func); // interpolation of func
  void export_fg(string filename); // export func to data file
  //virtual double H(vector<double> &X, vector<double> &P); // Hubble parameter
  //virtual double V(vector<double> &X); // potential
  //virtual double VI(vector<double> &X, int I); // \partial_I V
  //virtual double metric(vector<double> &X, int I, int J); // field-space metric G_IJ
  //virtual double inversemetric(vector<double> &X, int I, int J); // inverse field-space metric G^IJ
  //virtual double affine(vector<double> &X, int I, int J, int K); // Christoffel symbol Gamma^I_JK
  //virtual double derGamma(vector<double> &X, int I, int J, int K, int L); // Gamma^I_{JK,L}
  virtual double DI(int xp, int I, vector< vector<double> > &psv); 
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double CC(int num, vector< vector<double> > &psv, int func);
  /*
    solve PDE :
    (DI(xp,I) \partial_xpI + 1./2 DIJ(xpI,xpJ) \partial_xpI \partial_xpJ) f = CC
  */
  virtual void BoundaryCondition(); // set boundary condition
  virtual bool EndSurface(vector< vector<double> > &psv); // judge phase-space point psv is in the inflationary range or not
};

#endif
