#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

#define LICZBA_ITERACJI 100
int iteracja = 1;
double tol_x = 0.000000001; //zadana tolerancja błędu
double tol_f = 0.000000001; //zadana tolerancja residuum

double Jacobi[3][3] = {};

double f1(double x, double y, double z) { return x * x + y * y + z * z - 2; }
double f2(double x, double y, double z) { return x * x + y * y - 1; }
double f3(double x, double y, double z) { return x * x - y; }
double f1x(double x, double y, double z) { return 0; }
double f1y(double x, double y, double z) { return 1 / (2 * x * (1 + 2 * y)); }
double f1z(double x, double y, double z) { return y / (x * (1 + 2 * y)); }
double f2x(double x, double y, double z) { return 0; }
double f2y(double x, double y, double z) { return 1 / (1 + 2 * y); }
double f2z(double x, double y, double z) { return -1 / (1 + 2 * y); }
double f3x(double x, double y, double z) { return 1 / (2 * z); }
double f3y(double x, double y, double z) { return 1 / (-2 * z); }
double f3z(double x, double y, double z) { return 0; }

void oblicz_macierz_Jacobiego(double x, double y, double z)
{
    Jacobi[0][0] = f1x(x, y, z);
    Jacobi[0][1] = f1y(x, y, z);
    Jacobi[0][2] = f1z(x, y, z);
    Jacobi[1][0] = f2x(x, y, z);
    Jacobi[1][1] = f2y(x, y, z);
    Jacobi[1][2] = f2z(x, y, z);
    Jacobi[2][0] = f3x(x, y, z);
    Jacobi[2][1] = f3y(x, y, z);
    Jacobi[2][2] = f3z(x, y, z);
}

bool comp(double a, double b)
{
    return (a < b);
}

int main()
{
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;
    double dx, dy, dz; //poprawka
    double x1, y1, z1;
    double fx1, fx2, fx3;
    double ex, ey, ez;
    double r1, r2, r3;
    double estymator_bledu;
    double residuum;

    cout << "iteracja \tx \t  y \t    z \t estymator bledu \t residuum" << endl;
    for (int i = 0; i < LICZBA_ITERACJI; i++)
    {
        oblicz_macierz_Jacobiego(x0, y0, z0);
        fx1 = f1(x0, y0, z0);
        fx2 = f2(x0, y0, z0);
        fx3 = f3(x0, y0, z0);

        //obliczanie poprawki
        dx = Jacobi[0][0] * fx1 + Jacobi[0][1] * fx2 + Jacobi[0][2] * fx3;
        dy = Jacobi[1][0] * fx1 + Jacobi[1][1] * fx2 + Jacobi[1][2] * fx3;
        dz = Jacobi[2][0] * fx1 + Jacobi[2][1] * fx2 + Jacobi[2][2] * fx3;

        //odejmowanie poprawki
        x1 = x0 - dx;
        y1 = y0 - dy;
        z1 = z0 - dz;

        //obliczanie estymatora bledu
        ex = fabs(x1 - x0) / 2;
        ey = fabs(y1 - y0) / 2;
        ez = fabs(z1 - z0) / 2;
        estymator_bledu = max({ex, ey, ez}, comp);

        //obliczanie residuum
        r1 = fabs(f1(x1, y1, z1));
        r2 = fabs(f2(x1, y1, z1));
        r3 = fabs(f3(x1, y1, z1));
        residuum = max({r1, r2, r3}, comp);

        cout << iteracja << "\t\t" << x1 << "  " << y1 << "  " << z1 <<"\t " << estymator_bledu << "\t\t " << residuum << endl;
        if (estymator_bledu <= tol_x)
        {
            cout << "\nSTOP: Estymator bledu mniejszy od zadanej tolerancji\n";
            break;
        }
        if (residuum <= tol_f)
        {
            cout << "\nSTOP: Residuum mniejsze od zadanej tolerancji\n";
            break;
        }
        x0 = x1;
        y0 = y1;
        z0 = z1;
        iteracja++;

    }
}