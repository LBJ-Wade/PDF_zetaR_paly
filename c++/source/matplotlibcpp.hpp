#ifndef INCLUDED_matplotlibcpp_hpp_
#define INCLUDED_matplotlibcpp_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

using namespace std;

class matplotlibcpp {
  FILE *p;

public:
  matplotlibcpp(){}
  void open();
  void xlabel(string label);
  void ylabel(string label);
  void plot(vector<double> X, vector<double> Y, double lw, string color);
  void xlog();
  void ylog();
  void contourf(vector<double> X, vector<double> Y, vector<double> Z, string label);
  void log_contourf(vector<double> X, vector<double> Y, vector<double> Z, string label);
  void show();
  void save(string filename);
  void close();
};

void matplotlibcpp::open()
{
  p = popen("python -c 'import code; import os; import sys; sys.stdout = sys.stderr = open(os.devnull, \"w\"); code.InteractiveConsole().interact()'", "w");
  fprintf(p, "import matplotlib.pyplot as plt\n");
}

void matplotlibcpp::xlabel(string label)
{
  fprintf(p, "plt.xlabel(\"%s\")\n", label.c_str());
}

void matplotlibcpp::ylabel(string label)
{
  fprintf(p, "plt.ylabel(\"%s\")\n", label.c_str());
}

void matplotlibcpp::plot(vector<double> X, vector<double> Y, double lw, string color)
{
  if (X.size() == Y.size()) {
    if (!isnan(X[0]) && !isnan(Y[0])) {
      fprintf(p, "plt.plot([%f", X[0]);
    }
    for (int i=1; i<X.size(); i++) {
      if (!isnan(X[i]) && !isnan(Y[i])) {
	fprintf(p, ",%e", X[i]);
      }
    }
    if (!isnan(X[0]) && !isnan(Y[0])) {
      fprintf(p, "],[%f", Y[0]);
    }
    for (int i=1; i<Y.size(); i++) {
      if (!isnan(X[i]) && !isnan(Y[i])) {
	fprintf(p, ",%e", Y[i]);
      }
    }
    fprintf(p, "], lw=%f, color=\"%s\")\n", lw, color.c_str());
  }
}

void matplotlibcpp::xlog()
{
  fprintf(p, "plt.xscale(\"log\")\n");
}

void matplotlibcpp::ylog()
{
  fprintf(p, "plt.yscale(\"log\")\n");
}

void matplotlibcpp::contourf(vector<double> X, vector<double> Y, vector<double> Z, string label)
{
  int xsize = X.size();
  int ysize = Y.size();
  
  fprintf(p, "import numpy as np\n");
  fprintf(p, "x = np.array([%f", X[0]);
  for (int i=1; i<xsize; i++) {
    fprintf(p, ",%f", X[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "y = np.array([%f", Y[0]);
  for (int i=1; i<ysize; i++) {
    fprintf(p, ",%f", Y[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "X, Y = np.meshgrid(x,y)\n");
  fprintf(p, "Z = np.zeros((%d,%d))\n", ysize, xsize);
  for (int i=0; i<Z.size(); i++) {
    fprintf(p, "Z[%d,%d] = %f\n", i/xsize, i%xsize, Z[i]);
  }
  fprintf(p, "CF = plt.contourf(X,Y,Z)\n");
  fprintf(p, "CB = plt.colorbar(CF)\n");
  fprintf(p, "CB.set_label(\"%s\")\n", label.c_str());
}

void matplotlibcpp::log_contourf(vector<double> X, vector<double> Y, vector<double> Z,
				 string label)
{
  int xsize = X.size();
  int ysize = Y.size();
  
  fprintf(p, "import numpy as np\n");
  fprintf(p, "x = np.array([%f", X[0]);
  for (int i=1; i<xsize; i++) {
    fprintf(p, ",%f", X[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "y = np.array([%f", Y[0]);
  for (int i=1; i<ysize; i++) {
    fprintf(p, ",%f", Y[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "X, Y = np.meshgrid(x,y)\n");
  fprintf(p, "Z = np.zeros((%d,%d))\n", ysize, xsize);
  for (int i=0; i<Z.size(); i++) {
    fprintf(p, "Z[%d,%d] = %f\n", i/xsize, i%xsize, log10(Z[i]));
  }
  fprintf(p, "CF = plt.contourf(X,Y,Z)\n");
  fprintf(p, "CB = plt.colorbar(CF)\n");
  fprintf(p, "CB.set_label(\"%s\")\n", label.c_str());
}

void matplotlibcpp::show()
{
  fprintf(p, "plt.show()\n");
}

void matplotlibcpp::save(string filename)
{
  fprintf(p, "plt.savefig(\"%s\")\n", filename.c_str());
}

void matplotlibcpp::close()
{
  fprintf(p, "plt.close()");
  fprintf(p, "quit()");
}
  
#endif
