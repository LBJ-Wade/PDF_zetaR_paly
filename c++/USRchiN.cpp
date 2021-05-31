#include "source/JacobiPDEforchi.hpp"
#include <sys/time.h>

#define MODEL "USRchiN"

#define HH (1e-8)
#define DOTX0 (1./(5e+8)/sqrt(2))
#define DOTXF (1./(2e+14)/sqrt(10)/M_PI)
#define XF ((DOTX0-DOTXF)/3/HH)

#define XMIN -20
#define XMAX 30
#define NMIN 0
#define NMAX 8 //40
#define ITMIN 0.1
#define ITMAX 0.1
#define HIT 1
#define HX 0.1
#define HN 0.1

#define MAXSTEP (1e+6)
#define TOL (1e-10)

double Bf(double NN);


int main(int argc, char** argv)
{
  srand(time(NULL));
  
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6; // start stop watch


  double h = HX, sitev = XMIN;
  vector<double> site, it;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  h = HN, sitev = NMIN;
  while (sitev <= NMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();

  h = HIT;
  double itv = ITMIN;
  while (itv <= ITMAX) {
    it.push_back(itv);
    itv += h;
  }
  double funcNo = it.size();
  
  vector<double> params = {MAXSTEP,TOL,funcNo};

  for (int i=0; i<funcNo; i++) {
    cout << it[i] << ' ';
  }
  cout << endl;

  JacobiPDE pde(sitepack,it,params);
  for (int i=0; i<funcNo; i++) {
    pde.PDE_solve(i);
  }
  string model = MODEL;
  string str = "chi_" + model + ".dat";
  pde.export_fg(str);
    

  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
}


double JacobiPDE::DD(vector< vector<double> > &psv, int func) {
  return (*it)[func];
}

double JacobiPDE::DI(int xp, int I, vector< vector<double> > &psv, int func) {
  if (I == 0) {
    return 0;
  } else {
    return 1;
  }
}

double JacobiPDE::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv, int func) {
  if (I == 0 && J == 0) {
    return 1;
  } else {
    return 0;
  }
}

double JacobiPDE::CC(int num, vector< vector<double> > &psv, int func) {
  return 0;
}

void JacobiPDE::BoundaryCondition() // set boundary condition
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int number=0; number<volume; number++) {
    vector< vector<double> > PSV0(xpdim, vector<double>(Idim,0)); // temporal variable for phase-space value

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	PSV0[xp][I] = No2PSV(number,xp,I); // extract phase-space value of site[number]
      }
    }
    
    if (EndSurface(PSV0)) { // if the site is in inflationary region
      Omega[number] = true; // to be solved
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = (rand()%21)/10. - 1; //rand()%1; // set IC for function f randomly
      }
    } else {
      Omega[number] = false; // no to be solved
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = 1; // Set f to be 0. Particularly f should be 0 on the end of inflation hypersurface
      }
    }
    
    
    vector< vector<int> > index(xpdim, vector<int>(Idim,0));
    vector< vector<int> > ind_p, ind_m, ind_pp, ind_pm, ind_mm;

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	index[xp][I] = No2Ind(number,xp,I);
      }
    }

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	ind_p = index;
	ind_m = index;
	
	if (index[xp][I] == 0) {
	  ind_m[xp][I]++;
	  (*hm)[number][xp][I] = (*hI)[xp][I][index[xp][I]];
	} else {
	  ind_m[xp][I]--;
	  (*hm)[number][xp][I] = (*hI)[xp][I][index[xp][I]-1];
	}

	if (index[xp][I] == (*siteNo)[xp][I]-1) {
	  ind_p[xp][I]--;
	  (*hp)[number][xp][I] = (*hI)[xp][I][index[xp][I]-1];
	} else {
	  ind_p[xp][I]++;
	  (*hp)[number][xp][I] = (*hI)[xp][I][index[xp][I]];
	}

	(*num_m)[number][xp][I] = Ind2No(ind_m);
	(*num_p)[number][xp][I] = Ind2No(ind_p);

	
	for (int xptemp=0; xptemp<xpdim; xptemp++) {
	  for (int J=0; J<Idim; J++) {
	    if (xp!=xptemp || I!=J) {
	      ind_pp = index;
	      ind_pm = index;
	      ind_mm = index;
	      
	      if (index[xp][I] == 0) {
		ind_mm[xp][I]++;
	      } else {
		ind_mm[xp][I]--;
	      }

	      if (index[xp][I] == (*siteNo)[xp][I]-1) {
		ind_pp[xp][I]--;
		ind_pm[xp][I]--;
	      } else {
		ind_pp[xp][I]++;
		ind_pm[xp][I]++;
	      }

	      if (index[xptemp][J] == 0) {
		ind_pm[xptemp][J]++;
		ind_mm[xptemp][J]++;
	      } else {
		ind_pm[xptemp][J]--;
		ind_mm[xptemp][J]--;
	      }

	      if (index[xptemp][J] == (*siteNo)[xptemp][J]-1) {
		ind_pp[xptemp][J]--;
	      } else {
		ind_pp[xptemp][J]++;
	      }

	      (*num_pp)[number][xp][I][xptemp][J] = Ind2No(ind_pp);
	      (*num_pm)[number][xp][I][xptemp][J] = Ind2No(ind_pm);
	      (*num_mm)[number][xp][I][xptemp][J] = Ind2No(ind_mm);
	    }
	  }
	}
      }
    }
  }
}

bool JacobiPDE::EndSurface(vector< vector<double> > &psv)
{
  return psv[0][0] <= Bf(psv[0][1]);
}


double Bf(double NN) {
  return 2*M_PI*XF/HH - 2*M_PI*DOTX0/3/HH/HH*(1-exp(-3*NN));
}


