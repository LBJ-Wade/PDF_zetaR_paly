#include "source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "USRzetaR"

#define HH (1e-8)
#define DOTX0 (1./(5e+8)/sqrt(2))
#define DOTXF (1./(2e+14)/sqrt(10)/M_PI)
#define XF ((DOTX0-DOTXF)/3/HH)
//#define NF 40

#define XMIN -20
#define XMAX 30
#define NMIN 0
#define NMAX 8 //40
#define HX 0.1
#define HN 0.1

#define MAXSTEP (1e+6)
#define TOL (1e-10)

double Bf(double NN);


int main(int argc, char** argv)
{
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch

  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while(sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  h = HN, sitev = NMIN;
  while(sitev <= NMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();

  vector<double> params = {MAXSTEP,TOL,1,0,(double)sitepack[0].size(),0,0,0,0,0};

  vector< vector<double> > xpi = {{0,0}};

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params);
  sdn.PDE_solve(0);
  string model = MODEL;
  string str = "Mn_" + model + ".dat";
  sdn.export_fg(str);
 

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


double StocDeltaN::DI(int xp, int I, vector< vector<double> > &psv) {
  if (I == 0) {
    return 0;
  } else {
    return 1;
  }
}

double StocDeltaN::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv) {
  if (I == 0 && J == 0) {
    return 1;
  } else {
    return 0;
  }
}

bool StocDeltaN::EndSurface(vector< vector<double> > &psv) {
  return psv[0][0] <= Bf(psv[0][1]); //&& psv[0][1] <= NF;
}

double Bf(double NN) {
  return 2*M_PI*XF/HH - 2*M_PI*DOTX0/3/HH/HH*(1-exp(-3*NN));
}
