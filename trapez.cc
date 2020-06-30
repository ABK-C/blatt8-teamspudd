#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

// Testfunktionen als Funktor
class Pol1 {
public:
  double operator()(double x) { return 3 * x + 2; }
};

class Pol2 {
public:
  double operator()(double x) { return -2*pow(x,2)+3*x+1; }
};

class Gauss {
public:
  double operator()(double x) { return 1 / (sqrt(M_PI * 2)) * exp(-x * x / 2); }
};

// berechnet Werte nach Trapezformel von I_0 bis I_N

template <class Functor> std::vector<double> trapez(Functor f, double a, double b, int N) {
  std::vector<double> I(N + 1); // Feld mit N+1 Eintragen
  const double h = b - a;
  I[0] = h / 2 * (f(a) + f(b));
  for (int k = 1; k <= N; ++k) {
    int n = pow(2, k);
    double Summe = 0;
    for (int i = 1; i <= n-1; i++){
      Summe += f(a+i*h/n);
    }
    I[k] = I[0]/n + h/n*Summe;

  }
  return I;
}

// berechnet die Richardsonextrapolation aus I(k-1)  und I(k)
double richardson(double Iprev, double I) { return (4/3.0)*I - (1/3.0)*Iprev; }

// berechet Naeherungen ueber das Romberg-Verfahren
// I: Ergebnis von trapez()
std::vector<std::vector<double>> romberg(std::vector<double> I) {
  const int N = I.size() - 1;
  std::vector<std::vector<double>> R(N + 1);
  for (int k = 0; k <= N; ++k) {
    R[k].push_back(I[k]);
  }
  return R;
}


void testeAufgabe1() {
  Pol1 f;
  std::vector<double> I_f = trapez(f, 0, 3, 3);
  for (double tf : I_f) {
    std::cout << "A1: f:" << tf << " == " << 19.5 << ":" << (tf == 19.5 ? "ja" : "nein") << std::endl;
  }
  Pol2 g;
  std::cout << "A1: g(1) = 2 ?" << (g(1) == 2 ? "ja" : "nein") << std::endl;
  std::vector<double> I_g = trapez(g, 0, 3, 3);
  std::cout << "A1: g0:" << I_g[0] << " == " << -10.5 << ":" << (I_g[0] == -10.5 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g1:" << I_g[1] << " == " << -3.75 << ":" << (I_g[1] == -3.75 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g2:" << I_g[2] << " == " << -2.0625 << ":" << (I_g[2] == -2.0625 ? "ja" : "nein") << std::endl;
  double rich = richardson(I_g[0], I_g[1]);
  std::cout << "A1: Richardson : " << rich << " : " << (rich == -1.5 ? "ja " : "nein") << std::endl;
}

void testeAufgabe2() {
  Pol1 f;
  std::vector<std::vector<double>> Rf = romberg(trapez(f, 0, 3, 3));
  bool alle_richtig = true;
  int entries = 0;
  for (auto row : Rf) {
    for (double val : row) {
      alle_richtig &= val == 19.5;
      ++entries;
    }
  }
  std::cout << "A2: alle Eintraege für f sind 1.5:" << (alle_richtig ? "ja" : "nein") << std::endl;
  std::cout << "A2: korrekte Zahl an Einträgen:" << (entries == 10 ? "ja" : "nein") << std::endl;
  Pol2 g;
  std::vector<std::vector<double>> Rg = romberg(trapez(g, 0, 3, 3));
  std::cout << "A2: R[1][1] und R[2][1] für g gleich -1.5: " << ((Rg[1][1] == -1.5) && (Rg[2][1] == -1.5) ? " ja " : " nein") << std::endl;
}


int main() {
  // Testfunktion:
  Pol1 f;
  Pol2 g;
  Gauss e;
  std::cout << "f(0) = " << f(0) << '\n';
  // berechne Trapezformel fuer f
  std::vector<double> tf = trapez(f, 0., 3., 3);
  std::vector<double> tg = trapez(g, 0., 3., 3);
  std::vector<double> te = trapez(e, 0., 3., 3);
  std::cout
      << "#############################################################\n";
  // Ausgabe:
  std::cout << "Trapez:\n";
  for (unsigned int i = 0; i < tf.size(); ++i) { // Schleife ueber Werte im Feld
    std::cout << "I_" << i << " = " << tf[i] << std::endl;
  }
  std::cout << "Romberg:\n";
  std::vector<std::vector<double> > R = romberg(tf);
  for(int k = 0, l = R.size() ;  k < l ; ++k) {
    for(int n = 0, m = R[k].size() ;  n < m ; ++n) {
      std::cout << R[k][n] << " ";
    }
    std::cout << std::endl;
  }
  using namespace std;
  ofstream fout("Ergebnisse.txt");
  fout << "f(x)           " << "I_0(f)        " << "I_1(f)        " << "I_2(f)        " << "I_3(f)        " << "Integral von 0 bis 3" << std::endl;
  fout << "3x + 2         " << tf[0] << "          " << tf[1] << "           " << tf[2] << "         " << tf[3] << "           " << "19,5" <<std::endl; 
  fout << "-2x² + 3x + 1  " << tg[0] << "         " << tg[1] << "         " << tg[2] << "      " << tg[3] << "        " << "-1,5" <<std::endl; 
  fout << "Gauss          " << te[0] << "      " << te[1] << "       " << te[2] << "     " << te[3] << "       " << "0,49865" <<std::endl; 
  
  fout << "\n"  << "f(x)           " << "4/3I_1 -1/3I_0   " << "4/3I_2 -1/3I_1   " << "4/3I_3 -1/3I_2   " << "Integral von 0 bis 3" << std::endl;  
  fout << "3x + 2         " << richardson(tf[0], tf[1]) << "             " << richardson(tf[1], tf[2])  << "             " << richardson(tf[2], tf[3])  << "             " << "19,5" <<std::endl;
   fout << "-2x² + 3x + 1  " << richardson(tg[0], tg[1]) << "             " << richardson(tg[1], tg[2]) << "             " << richardson(tg[2], tg[3]) << "             " << "-1,5" <<std::endl; 
  fout << "Gauss          " << richardson(te[0], te[1]) << "         " << richardson(te[1], te[2]) << "         " << richardson(te[2], te[3]) << "         " << "0,49865" <<std::endl;

  testeAufgabe1();
  testeAufgabe2();
  
}
