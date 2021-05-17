#include "SRK32.hpp"
#include "MT.h"

// ------------------- user decision -----------------------
// ---------------------------------------------------------

double SRKintegrater::H(vector<double> &X, vector<double> &P)
{
  double rho = V(X);

  for (int I=0; I<X.size(); I++) {
    for (int J=0; J<X.size(); J++) {
      rho += 1./2*inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  return sqrt(rho/3.);
}

double SRKintegrater::V(vector<double> &X) 
{
  double m1 = 0.01;
  double m2 = 0.1;

  return 1./2*m1*m1*X[0]*X[0] + 1./2*m2*m2*X[1]*X[1];
}

double SRKintegrater::VI(vector<double> &X, int I) 
{
  double m1 = 0.01;
  double m2 = 0.1;

  if (I == 0) {
    return m1*m1*X[0];
  } else {
    return m2*m2*X[1];
  }
}

double SRKintegrater::metric(vector<double> &X, int I, int J)
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

double SRKintegrater::inversemetric(vector<double> &X, int I, int J)
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

double SRKintegrater::affine(vector<double> &X, int I, int J, int K)
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

double SRKintegrater::derGamma(vector<double> &X, int I, int J, int K, int L) {
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

double SRKintegrater::DI(int xp, int I, vector< vector<double> > &psv)
{
  double DI = 0;

  if (xpdim == 1) {
    for (int J=0; J<Idim; J++) {
      DI -= inversemetric(psv[0],I,J)*VI(psv[0],J)/V(psv[0]);

      for (int K=0; K<Idim; K++) {
	DI -= 1./2*affine(psv[0],I,J,K)*DIJ(0,J,0,K,psv);
      }
    }
  } else if (xpdim == 2) {
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

double SRKintegrater::DIJ(int xpI, int I, int xpJ, int J,
			  vector< vector<double> > &psv)
{
  double DDIJ;
  
  if (xpdim == 1) {
    DDIJ = V(psv[0])/12./M_PI/M_PI * inversemetric(psv[0],I,J);
  } else if (xpdim == 2) {
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

double SRKintegrater::gIa(int xp, int I, int alpha, vector< vector<double> > &psv)
{
  double ggIa;
  
  if (xpdim == 1) {
    ggIa = sqrt(V(psv[0])/3.)/2./M_PI * vielbein(psv,I,alpha);
  } else if (xpdim == 2) {
    if (xp == 0) {
      ggIa = H(psv[0],psv[1])/2./M_PI * vielbein(psv,I,alpha);
    } else if (xp == 1) {
      ggIa = 0;
      
      for (int J=0; J<Idim; J++) {
	for (int K=0; K<Idim; K++) {
	  ggIa += affine(psv[0],K,I,J)*psv[1][K]*gIa(0,J,alpha,psv);
	}
      }
    }
  }

  return ggIa;
}

// ---------------------------------------------------------
// ---------------------------------------------------------



SRKintegrater::SRKintegrater(vector< vector<double> > &XPi, double T0, int NoiseDim)
{
  init_genrand((unsigned)time(NULL));

  t = T0;
  xx = XPi;
  t0 = T0;
  xxi = XPi;

  xpdim = xx.size();
  Idim = xx[0].size();
  noisedim = NoiseDim;

  dW = vector<double>(noisedim,0);
  aIs = vector< vector<double> >(Idim, vector<double>(Idim,0));
  vIs = aIs;
  for (int I=0; I<Idim; I++) {
    for (int J=0; J<Idim; J++) {
      aIs[I][J] = rand_normal(0,1);
    }
  }

  for (int i=0; i<3; i++) {
    C0[i] = 0;
    C1[i] = 0;
    Alpha[i] = 0;
    Beta1[i] = 0;
    Beta2[i] = 0;
    aa[i] = vector< vector<double> >(xpdim, vector<double>(Idim,0));
    H0[i] = aa[i];
    bb[i] = vector< vector< vector<double> > >(noisedim,vector< vector<double> >(xpdim,vector<double>(Idim,0)));
    Hk[i] = bb[i];

    for (int j=0; j<3; j++) {
      A0[i][j] = 0;
      A1[i][j] = 0;
      B0[i][j] = 0;
      B1[i][j] = 0;
    }
  }

  A0[1][0] = 1;
  A1[1][0] = 1;
  A1[2][0] = 1;
  B1[1][0] = 1;
  B1[2][0] = -1;
  C0[1] = 1;
  C1[1] = 1;
  C1[2] = 1;
  Alpha[0] = 1./2;
  Alpha[1] = 1./2;
  Beta1[0] = 1;
  Beta2[1] = 1./2;
  Beta2[2] = -1./2;
}

void SRKintegrater::SRK2(double dt)
{
  for (int alpha=0; alpha<noisedim; alpha++) {
    dW[alpha] = rand_normal(0.,1.)*sqrt(dt);
  }

  for (int i=0; i<3; i++) {
    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	H0[i][xp][I] = xx[xp][I];

	for (int j=0; j<3; j++) {
	  H0[i][xp][I] += A0[i][j]*aa[j][xp][I]*dt;

	  for (int alpha=0; alpha<noisedim; alpha++) {
	    H0[i][xp][I] += B0[i][j]*bb[j][alpha][xp][I]*dW[alpha];
	  }
	}

	for (int beta=0; beta<noisedim; beta++) {
	  Hk[i][beta][xp][I] = xx[xp][I];

	  for (int j=0; j<3; j++) {
	    Hk[i][beta][xp][I] += A1[i][j]*aa[j][xp][I]*dt;

	    for (int alpha=0; alpha<noisedim; alpha++) {
	      Hk[i][beta][xp][I] += B1[i][j]*bb[j][alpha][xp][I]*dW[alpha]*dW[beta]
		/2./sqrt(dt);

	      if (alpha == beta) {
		Hk[i][beta][xp][I] -= B1[i][j]*bb[j][alpha][xp][I]*sqrt(dt)/2.;
	      }
	    }
	  }
	}
      }
    }

    coeff(dt,i);
  }

  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      for (int i=0; i<3; i++) {
	xx[xp][I] += Alpha[i]*aa[i][xp][I]*dt;

	for (int alpha=0; alpha<noisedim; alpha++) {
	  xx[xp][I] += (Beta1[i]*dW[alpha] + Beta2[i]*sqrt(dt))*bb[i][alpha][xp][I];
	}
      }
    }
  }

  t += dt;
}

void SRKintegrater::coeff(double dt, int step)
{
  double T = t + C0[step]*dt;
  vector< vector<double> > XP = H0[step];
  
  for (int xp=0; xp<xpdim; xp++) {
    for (int I=0; I<Idim; I++) {
      aa[step][xp][I] = DI(xp,I,XP);
    }
  }

  T = t + C1[step]*dt;

  for (int alpha=0; alpha<noisedim; alpha++) {
    XP = Hk[step][alpha];

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	bb[step][alpha][xp][I] = gIa(xp,I,alpha,XP);
      }
    }
  }
}

void SRKintegrater::init_txp()
{
  t = t0;
  xx = xxi;
}

void SRKintegrater::set_txp(double T, vector< vector<double> > &PSV)
{
  t = T;
  xx = PSV;
}

double SRKintegrater::e1(vector<double> &X, vector<double> &P)
{
  double e1 = 0;

  for (int I=0; I<Idim; I++) {
    for (int J=0; J<Idim; J++) {
      e1 += inversemetric(X,I,J)*P[I]*P[J];
    }
  }

  e1 /= 2*H(X,P)*H(X,P);
  
  return e1;
}

double SRKintegrater::return_H()
{
  return H(xx[0],xx[1]);
}

double SRKintegrater::return_V()
{
  return V(xx[0]);
}

double SRKintegrater::return_t()
{
  return t;
}

double SRKintegrater::return_xp(int xp, int I)
{
  return xx[xp][I];
}

double SRKintegrater::return_e1()
{
  return e1(xx[0],xx[1]);
}

double SRKintegrater::vielbein(vector< vector<double> > &XP, int I, int alpha)
{
  if (alpha == 0) {
    return eIsigma(XP,I);
  } else {
    return eIs(XP,I,alpha);
  }
}

double SRKintegrater::eIsigma(vector< vector<double> > &XP, int I)
{
  double NormN = 0, eIsigma = 0;
  
  if (xpdim == 1) {
    for (int J=0; J<Idim; J++) {
      for (int K=0; K<Idim; K++) {
	NormN += inversemetric(XP[0],J,K)*VI(XP[0],J)*VI(XP[0],K);
      }
    }

    for (int J=0; J<Idim; J++) {
      eIsigma += inversemetric(XP[0],I,J)*VI(XP[0],J);
    }
  } else if (xpdim == 2) {
    for (int J=0; J<Idim; J++) {
      for (int K=0; K<Idim; K++) {
	NormN += inversemetric(XP[0],J,K)*XP[1][J]*XP[1][K];
      }
    }

    for (int J=0; J<Idim; J++) {
      eIsigma += inversemetric(XP[0],I,J)*XP[1][J];
    }
  }
  
  return eIsigma/sqrt(NormN);
}

double SRKintegrater::eIs(vector< vector<double> > &XP, int I, int alpha)
{
  vIs = aIs;

  double Norm;
  for (int al=0; al<alpha; al++) {
    Norm = 0;
    for (int J=0; J<Idim; J++) {
      for (int K=0; K<Idim; K++) {
	Norm += metric(XP[0],J,K)*vielbein(XP,J,al)*aIs[K][alpha];
      }
    }

    for (int J=0; J<Idim; J++) {
      vIs[J][alpha] -= Norm*vielbein(XP,J,al);
    }
  }

  Norm = 0;
  for (int J=0; J<Idim; J++) {
    for (int K=0; K<Idim; K++) {
      Norm += metric(XP[0],J,K)*vIs[J][alpha]*vIs[K][alpha];
    }
  }
  
  return vIs[I][alpha]/sqrt(Norm);
}

// only for 2-field model
double SRKintegrater::eta_perp(vector< vector<double> > &XP)
{
  double Vs = 0;

  for (int I=0; I<Idim; I++) {
    Vs += VI(XP[0],I)*eIs(XP,I,1);
  }

  double sigmadot = 0;

  if (xpdim == 1) {
    for (int I=0; I<Idim; I++) {
      for (int J=0; J<Idim; J++) {
	sigmadot += inversemetric(XP[0],I,J)*VI(XP[0],I)*VI(XP[0],J)/3./V(XP[0]);
      }
    }
  } else {
    for (int I=0; I<Idim; I++) {
      for (int J=0; J<Idim; J++) {
	sigmadot += inversemetric(XP[0],I,J)*XP[1][I]*XP[1][J]; ////////////
      }
    }
  }

  sigmadot = sqrt(sigmadot);

  if (xpdim == 1) {
    return Vs/sqrt(V(XP[0])/3.)/sigmadot;
  } else {
    return Vs/H(XP[0],XP[1])/sigmadot;
  }
}

double SRKintegrater::return_etaperp()
{
  return eta_perp(xx);
}


double Uniform()
{
  return genrand_real3();
}

double rand_normal(double mu, double sigma)
{
  double z = sqrt(-2.*log(Uniform())) * sin(2.*M_PI*Uniform());
  return mu + sigma*z;
}
