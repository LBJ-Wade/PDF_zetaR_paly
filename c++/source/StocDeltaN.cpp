#include "StocDeltaN.hpp"
#include "matplotlibcpp.hpp"

StocDeltaN::StocDeltaN(string Model, vector< vector< vector<double> > > &Site,
		       vector< vector<double> > &XPi, double T0, vector<double> &Params):
  JacobiPDE(Site,Params), SRKintegrater(XPi,T0,Params[4])
{
  model = Model;
  timestep = Params[5];
  Nmax = Params[6];
  deltaN = Params[7];
  recursion = Params[8];
  xpdim = JacobiPDE::xpdim;
  Idim = JacobiPDE::Idim;
  //BoundaryCondition(Params[9]);
  BoundaryCondition();

  cout << "model : " << model << endl;
}

/*
void StocDeltaN::BoundaryCondition(double Ncut)
{
  double setnum = 0;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int number=0; number<volume; number++) {
    vector< vector<double> > PSV0(xpdim, vector<double>(Idim,0));

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	PSV0[xp][I] = No2PSV(number,xp,I);
      }
    }
    
    if (EndSurface(PSV0)) {
      Omega[number] = true;
      for (int func=0; func<funcNo; func++) {
	if (func==0) {
	  StocDeltaN ssdn = *this;
	  ssdn.set_txp(0,PSV0);

	  while (EndSurface(ssdn.xx)) {
	    ssdn.SRK2(timestep);

	    if (ssdn.return_t() > Ncut) {
	      Omega[number] = false;
	      break;
	    }
	  }

	  if (Omega[number]) {
	    ff[func][number] = ssdn.return_t();
	  } else {
	    ff[func][number] = 0;
	  }
	} else {
	  ff[func][number] = rand()%1;
	}
      }
    } else {
      Omega[number] = false;
      for (int func=0; func<funcNo; func++) {
	ff[func][number] = 0;
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

#pragma omp critical
    {
      setnum++;
      cout << "\rBox initializing... " << (int)(setnum/volume*100) << "%" << flush;
    }
  }
  cout << endl;
}
*/

void StocDeltaN::solve()
{
  PDE_solve(0);
  PDE_solve(1);
  
  string str = "Mn_" + model + ".dat";

  export_fg(str);

  vector< vector<double> > dN2List[recursion];
  str = "traj_" + model + ".dat";
  ofstream trajfile(str);
  double dt = timestep;
  int recNo = 0;

  if (xpdim == 1) {
    if (Idim == 2) {
      x1traj.clear();
      x2traj.clear();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      x1traj.clear();
      p1traj.clear();
    } else if (Idim == 2) {
      x1traj.clear();
      x2traj.clear();
    }
  }

  init_txp();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<recursion; i++) {
    StocDeltaN ssdn = *this;
    vector< vector<double> > dN2data;
    double dN2, dataNo;

    while (true) {
      ssdn.SRK2(dt);

      if (i == 0) {
	trajfile << ssdn.return_t() << ' ';
	for (int xp=0; xp<xpdim; xp++) {
	  for (int I=0; I<Idim; I++) {
	    trajfile << ssdn.return_xp(xp,I) << ' ';
	  }
	}
	for (int func=0; func<funcNo; func++) {
	  trajfile << ssdn.return_intf(func) << ' ';
	}
	trajfile << endl;

	if (xpdim == 1) {
	  if (Idim == 2) {
	    x1traj.push_back(ssdn.return_xp(0,0));
	    x2traj.push_back(ssdn.return_xp(0,1));
	  }
	} else if (xpdim == 2) {
	  if (Idim == 1) {
	    x1traj.push_back(ssdn.return_xp(0,0));
	    p1traj.push_back(ssdn.return_xp(1,0));
	  } else if (Idim == 2) {
	    x1traj.push_back(ssdn.return_xp(0,0));
	    x2traj.push_back(ssdn.return_xp(0,1));
	  }
	}
      }

      dN2data.push_back({ssdn.return_intf(0),ssdn.return_intf(1)});

      vector< vector<double> > PSV0(xpdim, vector<double>(Idim,0));
      for (int xp=0; xp<xpdim; xp++) {
	for (int I=0; I<Idim; I++) {
	  PSV0[xp][I] = ssdn.return_xp(xp,I);
	}
      }

      if (!EndSurface(PSV0)) {
	break;
      }
    }

    for (double N=0; N<Nmax; N+=deltaN) {
      dN2 = 0;
      dataNo = 0;

      for (int list=0; list<dN2data.size(); list++) {
	if (N <= dN2data[list][0] && dN2data[list][0] < N+deltaN) {
	  dN2 += dN2data[list][1];
	  dataNo++;
	}
      }

      if (dataNo != 0) {
	dN2 /= dataNo;
      }
      dN2List[i].push_back({N,dN2});
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      recNo++;
      cout << "\r" << recNo << "/" << recursion << flush;
    }
  }
  cout << endl;

  double meandN2, predN2 = 0;
  str = "calP_" + model + ".dat";
  ofstream calPfile(str);

  for (int list=0; list<dN2List[0].size(); list++) {
    meandN2 = 0;
    recNo = 0;

    for (int i=0; i<recursion; i++) {
      meandN2 += dN2List[i][list][1];
      if (dN2List[i][list][1] != 0) {
	recNo++;
      }
    }

    meandN2 /= recNo;
    calPfile << dN2List[0][list][0] << ' ' << meandN2 << ' ' << (meandN2-predN2)/deltaN << endl;
    Ndata.push_back(dN2List[0][list][0]);
    calPdata.push_back((meandN2-predN2)/deltaN);
    predN2 = meandN2;
  }
}

void StocDeltaN::sample()
{
  init_txp();

  vector< vector<double> > xmax = xx, xmin = xx;
  double Vi, Hi;
  if (xpdim == 1) {
    Vi = return_V();
  } else if (xpdim == 2) {
    Hi = return_H();
  }
  
  string str = "sample_" + model + ".dat";
  ofstream ofs(str);
  double dt = timestep;

  while (true) {
    SRK2(dt);

    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	if (xx[xp][I] < xmin[xp][I]) {
	  xmin[xp][I] = xx[xp][I];
	}
	if (xx[xp][I] > xmax[xp][I]) {
	  xmax[xp][I] = xx[xp][I];
	}
      }
    }

    ofs << setprecision(6) << return_t() << ' ';
    for (int xp=0; xp<xpdim; xp++) {
      for (int I=0; I<Idim; I++) {
	ofs << return_xp(xp,I) << ' ';
      }
    }
    ofs << setprecision(17);
    if (xpdim == 1) {
      ofs << return_V() << ' ';
    } else if (xpdim == 2) {
      ofs << return_H() << ' ';
    }
    if (Idim == 2) {
      ofs << return_etaperp();
    }
    ofs << endl;

    if (xpdim == 1) {
      if (Idim == 1) {
	Ntraj.push_back(return_t());
	x1traj.push_back(return_xp(0,0));
      } else if (Idim == 2) {
	x1traj.push_back(return_xp(0,0));
	x2traj.push_back(return_xp(0,1));
      }
    } else if (xpdim == 2) {
      if (Idim == 1) {
	x1traj.push_back(return_xp(0,0));
	p1traj.push_back(return_xp(1,0));
      } else if (Idim == 2) {
	x1traj.push_back(return_xp(0,0));
	x2traj.push_back(return_xp(0,1));
      }
    }
    
    if (!EndSurface(xx)) {
      break;
    }
  }

  cout << "[xi, xf, xmin, xmax] : " << endl;
  for (int I=0; I<Idim; I++) {
    cout << "[" << xxi[0][I] << ", " << xx[0][I] << ", " << xmin[0][I] << ", "
	 << xmax[0][I] << "]" << endl;
  }
  cout << endl;

  if (xpdim == 2) {
    cout << "[pi, pf, pmin, pmax] : " << endl;
    for (int I=0; I<Idim; I++) {
      cout << "[" << xxi[1][I] << ", " << xx[1][I] << ", " << xmin[1][I] << ", "
	   << xmax[1][I] << "]" << endl;
    }
    cout << endl;
  }

  cout << "N = " << return_t() << endl;
  cout << setprecision(17);
  if (xpdim == 1) {
    cout << "Vi = " << Vi << ",  Vf = " << return_V();
  } else if (xpdim == 2) {
    cout << "Hi = " << Hi << ",  Hf = " << return_H();
  }
  cout << setprecision(6) << endl;
  if (xpdim == 2) {
    cout << "eH = " << return_e1() << endl;
  }
}

void StocDeltaN::sample_plot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (xpdim == 1) {
    if (Idim == 1) {
      g.xlabel(string("$N$"));
      g.ylabel(string("$\\phi$"));
      g.plot(Ntraj,x1traj,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$\\phi^2$"));
      g.plot(x1traj,x2traj,1,string("b"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      g.xlabel(string("$\\phi$"));
      g.ylabel(string("$\\pi$"));
      g.plot(x1traj,p1traj,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$\\phi^2$"));
      g.plot(x1traj,x2traj,1,string("b"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::sample_logplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (xpdim == 1) {
    if (Idim == 1) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$N$"));
      g.ylabel(string("$|\\phi|$"));
      g.ylog();
      g.plot(Ntraj,x1abs,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      vector<double> x2abs;
      for (auto& x2 : x2traj) {
	x2abs.push_back(fabs(x2));
      }
      
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.ylog();
      g.plot(x1traj,x2abs,1,string("b"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> p1abs;
      for (auto& p1 : p1traj) {
	p1abs.push_back(fabs(p1));
      }
      
      g.xlabel(string("$\\phi$"));
      g.ylabel(string("$|\\pi|$"));
      g.ylog();
      g.plot(x1traj,p1abs,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      vector<double> x2abs;
      for (auto& x2 : x2traj) {
	x2abs.push_back(fabs(x2));
      }
      
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.ylog();
      g.plot(x1traj,x2abs,1,string("b"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::sample_loglinearplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (xpdim == 1) {
    if (Idim == 2) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$\\phi^2$"));
      g.xlog();
      g.plot(x1abs,x2traj,1,string("b"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$|\\phi|$"));
      g.ylabel(string("$\\pi$"));
      g.xlog();
      g.plot(x1abs,p1traj,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$\\phi^2$"));
      g.xlog();
      g.plot(x1abs,x2traj,1,string("b"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::sample_loglogplot()
{
  string filename = "sample_" + model + ".pdf";
  matplotlibcpp g;
  g.open();

  if (xpdim == 1) {
    if (Idim == 2) {
      vector<double> x1abs, x2abs;
      for (int i=0; i<x1traj.size(); i++) {
	x1abs.push_back(fabs(x1traj[i]));
	x2abs.push_back(fabs(x2traj[i]));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.xlog();
      g.ylog();
      g.plot(x1abs,x2abs,1,string("b"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> x1abs, p1abs;
      for (int i=0; i<x1traj.size(); i++) {
	x1abs.push_back(fabs(x1traj[i]));
	p1abs.push_back(fabs(p1traj[i]));
      }
      
      g.xlabel(string("$|\\phi|$"));
      g.ylabel(string("$|\\pi|$"));
      g.xlog();
      g.ylog();
      g.plot(x1abs,p1abs,1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      vector<double> x1abs, x2abs;
      for (int i=0; i<x1traj.size(); i++) {
	x1abs.push_back(fabs(x1traj[i]));
	x2abs.push_back(fabs(x2traj[i]));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.xlog();
      g.ylog();
      g.plot(x1abs,x2abs,1,string("b"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::f_plot(int func)
{
  string filename;
  if (func == 0) {
    filename = "N_" + model + ".pdf";
  } else if (func == 1) {
    filename = "dN2_" + model + ".pdf";
  }
  matplotlibcpp g;
  g.open();
  
  if (xpdim == 1) {
    if (Idim == 1) {
      g.xlabel(string("$\\phi$"));
      if (func == 0) {
	g.ylabel(string("$<N>$"));
      } else if (func == 1) {
	g.ylabel(string("$<\\delta N^2>$"));
	g.ylog();
      }
      g.plot((*site)[0][0],ff[func],1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$\\phi^2$"));
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1traj,x2traj,3,string("r"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      g.xlabel(string("$\\phi$"));
      g.ylabel(string("$\\pi$"));
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[1][0],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[1][0],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1traj,p1traj,3,string("r"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::f_logplot(int func)
{
  string filename;
  if (func == 0) {
    filename = "N_" + model + ".pdf";
  } else if (func == 1) {
    filename = "dN2_" + model + ".pdf";
  }
  matplotlibcpp g;
  g.open();
  
  if (xpdim == 1) {
    if (Idim == 2) {
      vector<double> x2abs;
      for (auto& x2 : x2traj) {
	x2abs.push_back(fabs(x2));
      }
      
      g.xlabel(string("$\\phi^1$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.ylog();
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1traj,x2abs,3,string("r"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> p1abs, p1siteabs;
      for (auto& p1 : p1traj) {
	p1abs.push_back(fabs(p1));
      }
      for (auto& p1site : (*site)[1][0]) {
	p1siteabs.push_back(fabs(p1site));
      }
      
      g.xlabel(string("$\\phi$"));
      g.ylabel(string("$|\\pi|$"));
      g.ylog();
      if (func == 0) {
	g.contourf((*site)[0][0],p1siteabs,ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],p1siteabs,ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1traj,p1abs,3,string("r"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::f_loglinearplot(int func)
{
  string filename;
  if (func == 0) {
    filename = "N_" + model + ".pdf";
  } else if (func == 1) {
    filename = "dN2_" + model + ".pdf";
  }
  matplotlibcpp g;
  g.open();
  
  if (xpdim == 1) {
    if (Idim == 1) {
      g.xlabel(string("$\\phi$"));
      if (func == 0) {
	g.ylabel(string("$<N>$"));
      } else if (func == 1) {
	g.ylabel(string("$<\\delta N^2>$"));
	g.ylog();
      }
      g.xlog();
      g.plot((*site)[0][0],ff[func],1,string("b"));
      g.save(filename);
      g.show();
    } else if (Idim == 2) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$\\phi^2$"));
      g.xlog();
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1abs,x2traj,3,string("r"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> x1abs;
      for (auto& x1 : x1traj) {
	x1abs.push_back(fabs(x1));
      }
      
      g.xlabel(string("$|\\phi|$"));
      g.ylabel(string("$\\pi$"));
      g.xlog();
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[1][0],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[1][0],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1abs,p1traj,3,string("r"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::f_loglogplot(int func)
{
  string filename;
  if (func == 0) {
    filename = "N_" + model + ".pdf";
  } else if (func == 1) {
    filename = "dN2_" + model + ".pdf";
  }
  matplotlibcpp g;
  g.open();
  
  if (xpdim == 1) {
    if (Idim == 2) {
      vector<double> x1abs, x2abs;
      for (int i=0; i<x1traj.size(); i++) {
	x1abs.push_back(fabs(x1traj[i]));
	x2abs.push_back(fabs(x2traj[i]));
      }
      
      g.xlabel(string("$|\\phi^1|$"));
      g.ylabel(string("$|\\phi^2|$"));
      g.xlog();
      g.ylog();
      if (func == 0) {
	g.contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],(*site)[0][1],ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1abs,x2abs,3,string("r"));
      g.save(filename);
      g.show();
    }
  } else if (xpdim == 2) {
    if (Idim == 1) {
      vector<double> x1abs, p1abs, p1siteabs;
      for (int i=0; i<x1traj.size(); i++) {
	x1abs.push_back(fabs(x1traj[i]));
	p1abs.push_back(fabs(p1traj[i]));
      }
      for (auto& p1site : (*site)[1][0]) {
	p1siteabs.push_back(fabs(p1site));
      }
      
      g.xlabel(string("$|\\phi|$"));
      g.ylabel(string("$|\\pi|$"));
      g.xlog();
      g.ylog();
      if (func == 0) {
	g.contourf((*site)[0][0],p1siteabs,ff[func],string("$<N>$"));
      } else if (func == 1) {
	g.log_contourf((*site)[0][0],p1siteabs,ff[func],string("$<\\delta N^2>$"));
      }
      g.plot(x1abs,p1abs,3,string("r"));
      g.save(filename);
      g.show();
    }
  }
  g.close();
}

void StocDeltaN::calP_plot()
{
  string filename = "calP_" + model + ".pdf";
  matplotlibcpp g;
  g.open();
  g.xlabel(string("$<N>$"));
  g.ylabel(string("$\\mathcal{P}_\\zeta$"));
  g.ylog();
  g.plot(Ndata,calPdata,1,string("b"));
  g.save(filename);
  g.show();
  g.close();
}

double StocDeltaN::return_intf(int func)
{
  return Interpolation_f(xx,func);
}

