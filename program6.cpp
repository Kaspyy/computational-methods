#include <iostream>
#include <iomanip>
using namespace std;

void wypisz_wektor(double *wektor, int n)
{
    for (int i = 0; i < n; i++)
        cout << "|" << setw(10) << wektor[i] << "|" << endl;
}
void AlgorytmThomasa_macierz(double *L, double *D, double *U, int n)
{
    for (int i = 0; i < n; i++)
    {
        D[i + 1] = D[i + 1] - L[i] * (U[i] / D[i]);
    }
}
void AlgorytmThomasa_wektor(double *b, int n, double *L, double *D)
{
    for (int i = 1; i < n + 1; i++)
        b[i] = b[i] - (L[i - 1] / D[i - 1]) * b[i - 1];
}
void rozwiazanie(double *L, double *D, double *U, double *b, double *x, int n)
{
    x[n] = b[n] / D[n];
    for (int i = n - 1; i >= 0; i--)
        x[i] = (b[i] - U[i] * x[i + 1]) / D[i];
}
int main()
{
    int n = 5;
    double L[n] = {1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0, 1.0 / 11.0};
    double D[n + 1] = {10.0, 20.0, 30.0, 30.0, 20.0, 10.0};
    double U[n] = {1.0 / 2.0, 1.0 / 4.0, 1.0 / 6.0, 1.0 / 8.0, 1.0 / 10.0};
    double b[n + 1] = {31.0, 165.0 / 4.0, 917.0 / 30.0, 851.0 / 28.0, 3637.0 / 90.0, 332.0 / 11.0};
    double x[n + 1];

    cout << "Wektory: l, d, u, uzyskane z trojdiagonalnej macierzy A:" << endl
         << endl;
    cout << "Wektor dolnej przekatnej l: " << endl;
    wypisz_wektor(L, n);
    cout << endl
         << "Wektor srodkowej przekatnej d: " << endl;
    wypisz_wektor(D, n + 1);
    cout << endl
         << "Wektor gornej przekatnej u: " << endl;
    wypisz_wektor(U, n);
    cout << endl
         << "Wektor b: " << endl;
    wypisz_wektor(b, n + 1);
    cout << endl
         << endl;
    cout << "------------------Po wykonaniu algorytmu Thomasa--------------------------\n" << endl;
    AlgorytmThomasa_macierz(L, D, U, n);
    AlgorytmThomasa_wektor(b, n, L, D);
    cout << "Wektor L: " << endl;
    wypisz_wektor(L, n);
    cout << endl
         << "Wektor D: " << endl;
    wypisz_wektor(D, n + 1);
    cout << endl
         << "Wektor U: " << endl;
    wypisz_wektor(U, n);
    cout << endl
         << "Wektor B: " << endl;
    wypisz_wektor(b, n);

    rozwiazanie(L, D, U, b, x, n);
    cout << "\n\n--------------------Rozwiazanie------------------------\n" << endl;
    cout << endl
         << "Wektor x:" << endl;
    wypisz_wektor(x, n + 1);
    return 0;
}
