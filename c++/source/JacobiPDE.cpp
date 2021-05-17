#include "JacobiPDE.hpp"

// ------------------- user decision -----------------------
// ---------------------------------------------------------

double JacobiPDE::H(vector<double> &X, vector<double> &P) // Hubble parameter
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

double JacobiPDE::V(vector<double> &X) // potential
{
  double m1 = 0.01;
  double m2 = 0.1;
  
  return 1./2*m1*m1*X[0]*X[0] + 1./2*m2*m2*X[1]*X[1];
}

double JacobiPDE::VI(vector<double> &X, int I) // \partial_I V
{
  double m1 = 0.01;
  double m2 = 0.1;

  if (I == 0) {
    return m1*m1*X[0];
  } else {
    return m2*m2*X[1];
  }
}

double JacobiPDE::metric(vector<double> &X, int I, int J) // field-space metric G_IJ
{
  double MM = 1e-3;

  if (I == 0 && J == 0) {
    return 1 + 2*X[1]*X[1]/MM/MM;
  } else if (I == 1 && J == 1) {
    return 1;
  } else {
    return 0;
  }
}

double JacobiPDE::inversemetric(vector<double> &X, int I, int J) // inverse field-space metric G^IJ
{
  double MM = 1e-3;

  if (I == 0 && J == 0) {
    return 1./(1+2*X[1]*X[1]/MM/MM);
  } else if (I == 1 && J == 1) {
    return 1;
  } else {
    return 0;
  }
}

double JacobiPDE::affine(vector<double> &X, int I, int J, int K) // Christoffesl symbol Gamma^I_JK
{
  double MM = 1e-3;

  if (I == 0 && ((J == 0 && K == 1) || (J == 1 && K == 0))) {
    return 2*X[1] / (2*X[1]*X[1] + MM*MM);
  } else if (I == 1 && J == 0 && K == 0) {
    return -2*X[1]/MM/MM;
  } else {
    return 0;
  }
}

double JacobiPDE::derGamma(vector<double> &X, int I, int J, int K, int L) // Gamma^I_{JK,L}
{
  double MM = 1e-3;

  if (L == 1) {
    if (I == 0 && ((J == 0 && K == 1) || (J == 1 && K == 0))) {
      return 2*(MM*MM-2*X[1]*X[1]) / (MM*MM+2*X[1]*X[1]) / (MM*MM+2*X[1]*X[1]);
    } else if (I == 1 && J == 0 && K == 0) {
      return -2./MM/MM;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

// ---------------------------------------------------------
/* 
solve (DI(xp,I) \partial_xpI + 1./2 DIJ(xpI,xpJ) \partial_xpI \partial_xpJ) f = CC
func swithes f.
*/
double JacobiPDE::DI(int xp, int I, vector< vector<double> > &psv)
{
  double DI = 0;

  if (xpdim == 1) { // slow-roll field-space
    for (int J=0; J<Idim; J++) {
      DI -= inversemetric(psv[0],I,J)*VI(psv[0],J)/V(psv[0]);

      for (int K=0; K<Idim; K++) {
	DI -= 1./2*affine(psv[0],I,J,K)*DIJ(0,J,0,K,psv);
      }
    }
  } else if (xpdim == 2) { // full phase-space
    if (xp == 0) {
      DI = 0;
      
      for (int J=0; J<Idim; J++) {
	DI += inversemetric(psv[0],I,J)*psv[1][J]/H(psv[0],psv[1]);
	
	for (int K=0; K<Idim; K++) {
	  DI -= 1./2*affine(psv[0],I,J,K)*DIJ(0,J,0,K,psv);
	}
      }
    } else {
      double Hubble = H(psv[0],psv[1]);
      DI = -3*psv[1][I]-VI(psv[0],I)/Hubble;
      
      for (int J=0; J<Idim; J++) {
	for (int K=0; K<Idim; K++) {
	  DI += affine(psv[0],J,I,K)*DIJ(0,K,1,J,psv);
	  
	  for (int S=0; S<Idim; S++) {
	    DI += affine(psv[0],S,I,J)*inversemetric(psv[0],J,K)*psv[1][K]*psv[1][S]/Hubble
	      +1./2*derGamma(psv[0],S,I,J,K)*psv[1][S]*DIJ(0,J,0,K,psv);
	    
	    for (int R=0; R<Idim; R++) {
	      DI -= 1./2*(affine(psv[0],R,J,K)*affine(psv[0],S,I,R)
			  + affine(psv[0],R,I,J)*affine(psv[0],S,K,R))
		*psv[1][S]*DIJ(0,J,0,K,psv);
	    }
	  }
	}
      }
    }
  }
  
  return DI;
}

double JacobiPDE::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv)
{
  double DDIJ;
  
  if (xpdim == 1) { // slow-roll field-space
    DDIJ = V(psv[0])/12./M_PI/M_PI * inversemetric(psv[0],I,J);
  } else if (xpdim == 2) { // full phase-space
    if (xpI == 0 && xpJ == 0) { 
      DDIJ = H(psv[0],psv[1])*H(psv[0],psv[1])/4./M_PI/M_PI * inversemetric(psv[0],I,J);
    } else if (xpI == 1 && xpJ == 1) {
      DDIJ = 0;
      for (int K=0; K<Idim; K++) {
	for (int L=0; L<Idim; L++) {
	  for (int M=0; M<Idim; M++) {
	    for (int N=0; N<Idim; N++) {
	      DDIJ += affine(psv[0],K,I,L)*psv[1][K]*affine(psv[0],M,J,N)*psv[1][M]
		*DIJ(0,L,0,N,psv);
	    }
	  }
	}
      }
    } else if (xpI == 0) {
      DDIJ = 0;
      
      for (int K=0; K<Idim; K++) {
	for (int L=0; L<Idim; L++) {
	  DDIJ += affine(psv[0],K,J,L)*psv[1][K]*DIJ(0,I,0,L,psv);
	}
      }
    } else {
      DDIJ = DIJ(xpJ,J,xpI,I,psv);
    }
  }
  
  return DDIJ;
}

double JacobiPDE::CC(int num, vector< vector<double> > &psv, int func)
{
  double CC = 0;
  
  if (func == 0) { // for <N>
    CC = -1;
  } else if (func == 1) { // for <delta N^2>
    for (int xpI=0; xpI<xpdim; xpI++) {
      for (int I=0; I<Idim; I++) {
	for (int xpJ=0; xpJ<xpdim; xpJ++) {
	  for (int J=0; J<Idim; J++) {
	    CC -= DIJ(xpI,I,xpJ,J,psv)
	      *(ff[0][(*num_p)[num][xpI][I]] - ff[0][(*num_m)[num][xpI][I]])
	      *(ff[0][(*num_p)[num][xpJ][J]] - ff[0][(*num_m)[num][xpJ][J]])
	      /((*hp)[num][xpI][I]+(*hm)[num][xpI][I])/((*hp)[num][xpJ][J]+(*hm)[num][xpJ][J]);
	  }
	}
      }
    }
  }

  return CC;
}
// ---------------------------------------------------------

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
	ff[func][number] = rand()%1; // set IC for function f randomly
      }
    } else {
      Omega[number] = false; // no to be solved
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = 0; // Set f to be 0. Particularly f should be 0 on the end of inflation hypersurface
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
  if (xpdim == 1) {
    return V(psv[0]) >= rhoc;
  } else if (xpdim == 2) {
    return 3*H(psv[0],psv[1])*H(psv[0],psv[1]) >= rhoc;
  } else {
    return 0;
  }
}

// ---------------------------------------------------------
// ---------------------------------------------------------



JacobiPDE::JacobiPDE(vector< vector< vector<double> > > &Site, vector<double> &Params)
{
  srand((unsigned)time(NULL)); // initialize random seed

#ifdef _OPENMP
  cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#endif

  site = new vector< vector< vector<double> > >;
  hI = new vector< vector< vector<double> > >;
  (*site) = Site;
  (*hI) = Site;
  
  maxstep = Params[0];
  tol = Params[1];
  funcNo = Params[2];
  rhoc = Params[3];

  xpdim = (*site).size();
  Idim = (*site)[0].size();

  siteNo = new vector< vector<int> >;
  (*siteNo) = vector< vector<int> >(xpdim, vector<int>(Idim,0));
  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      (*siteNo)[xp][I] = (*site)[xp][I].size();

      for (int index=0; index<(*site)[xp][I].size()-1; index++) {
	(*hI)[xp][I][index] = (*site)[xp][I][index+1] - (*site)[xp][I][index];
      }
    }
  }

  volume = 1;
  xpvol = vector<int>(xpdim,1);
  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      volume *= (*siteNo)[xp][I];
      xpvol[xp] *= (*siteNo)[xp][I];
    }
  }

  cout << "total # of sites : " << volume << endl;

  ff = new vector<double>[funcNo];
  f_next = new vector<double>[funcNo];

  for (int i=0; i<funcNo; i++) {
    ff[i] = vector<double>(volume,0);
    f_next[i] = vector<double>(volume,0);
  }
  
  Omega = vector<bool>(volume,true);


  num_p = new vector< vector< vector<int> > >;
  num_m = new vector< vector< vector<int> > >;
  num_pp = new vector< vector< vector< vector< vector<int> > > > >;
  num_pm = new vector< vector< vector< vector< vector<int> > > > >;
  num_mm = new vector< vector< vector< vector< vector<int> > > > >;
  
  (*num_p) = vector< vector< vector<int> > >(volume, vector< vector<int> >(xpdim,
									vector<int>(Idim,0)));
  (*num_m) = (*num_p);
  (*num_pp) = vector< vector< vector< vector< vector<int> > > > >(volume, vector< vector< vector< vector<int> > > >(xpdim, vector< vector< vector<int> > >(Idim, vector< vector<int> >(xpdim, vector<int>(Idim,0)))));
  (*num_pm) = (*num_pp);
  (*num_mm) = (*num_pp);

  hp = new vector< vector< vector<double> > >;
  hm = new vector< vector< vector<double> > >;
  
  (*hp) = vector< vector< vector<double> > >(volume,
					  vector< vector<double> >(xpdim,
								   vector<double>(Idim,0)));
  (*hm) = (*hp);

  BoundaryCondition();
}

double JacobiPDE::PDE_1step(int num, int func)
{
  vector< vector<double> > PSV0(xpdim, vector<double>(Idim,0));

  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      PSV0[xp][I] = No2PSV(num,xp,I);
    }
  }

  double uu = CC(num,PSV0,func), coeff = 0;

  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      double DItemp = DI(xp,I,PSV0);
      
      if (DItemp < 0) {
	uu += DItemp / (*hm)[num][xp][I] * ff[func][(*num_m)[num][xp][I]];
	coeff += DItemp / (*hm)[num][xp][I];
      } else {
	uu -= DItemp / (*hp)[num][xp][I] * ff[func][(*num_p)[num][xp][I]];
	coeff -= DItemp / (*hp)[num][xp][I];
      }
      
      uu -= DIJ(xp,I,xp,I,PSV0)
	*(ff[func][(*num_p)[num][xp][I]]*(*hm)[num][xp][I]+ff[func][(*num_m)[num][xp][I]]*(*hp)[num][xp][I])
	/(*hp)[num][xp][I]/(*hm)[num][xp][I]/((*hp)[num][xp][I]+(*hm)[num][xp][I]);
      coeff -= DIJ(xp,I,xp,I,PSV0)/(*hp)[num][xp][I]/(*hm)[num][xp][I];

      for (int xptemp=0; xptemp<xpdim; xptemp++) {
	for (int J=0; J<Idim; J++) {
	  if (xp!=xptemp || I!=J) {
	    uu -= 1./2*DIJ(xp,I,xptemp,J,PSV0)
	      *(ff[func][(*num_pp)[num][xp][I][xptemp][J]]-ff[func][(*num_pm)[num][xp][I][xptemp][J]]
		-ff[func][(*num_pm)[num][xptemp][J][xp][I]]+ff[func][(*num_mm)[num][xp][I][xptemp][J]])
	      /((*hp)[num][xp][I]+(*hm)[num][xp][I])/((*hp)[num][xptemp][J]+(*hm)[num][xptemp][J]);
	  }
	}
      }
    }
  }

  uu /= coeff;
  
  return uu;
}

void JacobiPDE::PDE_solve(int func)
{
  double u_norm, err;

  for (int step=0; step<maxstep; step++) {
    u_norm = 0;
    err = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:u_norm, err)
#endif
    {
#ifdef _OPENMP
#pragma omp for
#endif
      for (int num=0; num<volume; num++) {
	if (Omega[num]) { 
	  u_norm += ff[func][num]*ff[func][num];
	  f_next[func][num] = PDE_1step(num,func);
	  err += (ff[func][num]-f_next[func][num])*(ff[func][num]-f_next[func][num]);
	}
      }
    }

    for (int num=0; num<volume; num++) {
      if (Omega[num]) {
	ff[func][num] = f_next[func][num];
      }
    }
    err = sqrt(err)/sqrt(u_norm);
    cout << "\rerr" << func+1 << " : " << setw(11) << left << err << "  step : " << step << flush;

    if (err < tol) {
      break;
    }
  }

  cout << endl;
}

int JacobiPDE::Ind2No(vector< vector<int> > &index)
{
  int number = 0, temp;

  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      temp = index[xp][I];
      for (int J=0; J<I; J++) {
	temp *= (*siteNo)[xp][J];
      }
      for (int xptemp=0; xptemp<xp; xptemp++) {
	temp *= xpvol[xptemp];
      }
      number += temp;
    }
  }

  return number;
}

int JacobiPDE::No2Ind(int num, int xp, int I)
{
  int index = num;

  for (int xptemp=0; xptemp<xp; xptemp++) {
    index /= xpvol[xptemp];
  }

  index %= xpvol[xp];

  for (int J=0; J<I; J++) {
    index /= (*siteNo)[xp][J];
  }

  index %= (*siteNo)[xp][I];

  return index;
}

double JacobiPDE::No2PSV(int num, int xp, int I)
{
  return (*site)[xp][I][No2Ind(num,xp,I)];
}

int JacobiPDE::ceilXP(int xp, int I, vector< vector<double> > &psv)
{
  int index = 0;

  while ((*site)[xp][I][index] <= psv[xp][I]) {
    index++;
  }

  return index;
}

double JacobiPDE::Interpolation_f(vector< vector<double> > &psv, int func)
{
  int pmcheck, xpsite = pow(2,Idim);
  double weight, intf = 0;
  vector< vector<int> > index(xpdim,vector<int>(Idim,0));
  
  for (int num=0; num<pow(2,Idim*xpdim); num++) {
    weight = 1;
    
    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	pmcheck = num;

	for (int xptemp=0; xptemp<xp; xptemp++) {
	  pmcheck /= xpsite;
	}
	pmcheck %= xpsite;

	for (int J=0; J<I; J++) {
	  pmcheck /= 2;
	}
	pmcheck %= 2;

	if (pmcheck == 1) {
	  index[xp][I] = ceilXP(xp,I,psv);
	  weight *= (psv[xp][I]-(*site)[xp][I][ceilXP(xp,I,psv)-1])/(*hI)[xp][I][ceilXP(xp,I,psv)-1];
	} else {
	  index[xp][I] = ceilXP(xp,I,psv)-1;
	  weight *= ((*site)[xp][I][ceilXP(xp,I,psv)]-psv[xp][I])/(*hI)[xp][I][ceilXP(xp,I,psv)-1];
	}
      }
    }

    intf += ff[func][Ind2No(index)]*weight;
  }

  return intf;
}

void JacobiPDE::export_fg(string filename)
{
  ofstream ofs(filename);

  for (int number=0; number<volume; number++) {
    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	ofs << No2PSV(number,xp,I) << ' ';
      }
    }

    for (int func=0; func<funcNo; func++) {
      ofs << ff[func][number] << ' ';
    }
    ofs << endl;
  }
}

