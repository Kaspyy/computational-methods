#include <iostream>
#include <iomanip>
#include <cmath>
#define liczba_iteracji 30
using namespace std;

double tol_x = 0.000001; //zadana tolerancja błędu
double estymator_bledu;
double tol_f = 0.000001; //zadana tolerancja residuum
double residuum;

double funkcja_1(double x) { return sin(x / 4) * sin(x / 4) - x; }
double funkcja_2(double x) { return tan(2 * x) - x - 1; }
double funkcja_1_pochodna(double x) { return 0.25 * sin(x / 2) - 1; }
double funkcja_2_pochodna(double x) { return 2 / (cos(2 * x) * cos(2 * x)) - 1; }
double phi_picard_1(double x) { return sin(x / 4) * sin(x / 4); }
double phi_picard_2(double x) { return tan(2 * x) - 1; }

void metoda_picarda(double (*phi)(double), double (*f)(double))
{
    double x0, x1;
    x0 = 1;
    if (fabs(phi(x0)) >= 1)
    {
        cout << "STOP: Rozbieznosc funkcji, Phi(x) >= 1" << endl;
        return;
    }
    else
    {
        cout << "\t x:\t\t  Estymator bledu:\t      residuum:\n";
        for (int i = 0; i < liczba_iteracji; i++)
        {
            x1 = phi(x0);
            estymator_bledu = fabs(x1 - x0);
            residuum = f(x1);
            cout << x1 << "\t" << estymator_bledu << "\t" << residuum << endl;
            if (fabs(estymator_bledu) <= tol_x)
            {
                cout << "\nSTOP: Estymator bledu mniejszy od zadanej tolerancji\n";
                break;
            }
            if (fabs(residuum) <= tol_f)
            {
                cout << "\nSTOP: residuum mniejsze od zadanej tolerancji\n";
                break;
            }
            if (i == liczba_iteracji)
            {
                cout << "\nSTOP: przekroczono liczbe iteracji\n";
                break;
            }
            x0 = x1;
        }
    }
}
void metoda_bisekcji(double (*f)(double), double a, double b)
{
    double xn;
    cout << "\t x:\t\t  Estymator bledu:\t      residuum:\n";
    for (int i = 0; i < liczba_iteracji; i++)
    {
        xn = (a + b) / 2;
        estymator_bledu = fabs(a - b) / 2.0;
        residuum = f(xn);
        cout << xn << "\t" << estymator_bledu << "\t" << residuum << endl;
        if (fabs(estymator_bledu) <= tol_x)
        {
            cout << "\nSTOP: Estymator bledu mniejszy od zadanej tolerancji\n";
            break;
        }
        if (fabs(residuum) <= tol_f)
        {
            cout << "\nSTOP: residuum mniejsze od zadanej tolerancji\n";
            break;
        }
        if (i == liczba_iteracji)
        {
            cout << "\nSTOP: przekroczono liczbe iteracji\n";
            break;
        }
        if (f(xn) > 0)
        {
            if (f(a) > 0)
                a = xn;
            else
                b = xn;
        }
        else
        {
            if (f(a) < 0)
                a = xn;
            else
                b = xn;
        }
    }
}
void metoda_newtona(double (*f)(double), double (*fp)(double), double start)
{
    double x0, x1;
    x0 = start;
    cout << "\t x:\t\t  Estymator bledu:\t      residuum:\n";
    for (int i = 0; i < liczba_iteracji; i++)
    {
        x1 = x0 - f(x0) / fp(x0);
        estymator_bledu = fabs(x0 - x1);
        residuum = f(x1);
        cout << x1 << "\t" << estymator_bledu << "\t" << residuum << endl;
        if (fabs(estymator_bledu) <= tol_x)
        {
            cout << "\nSTOP: Estymator bledu mniejszy od zadanej tolerancji\n";
            break;
        }
        if (fabs(residuum) <= tol_f)
        {
            cout << "\nSTOP: residuum mniejsze od zadanej tolerancji\n";
            break;
        }
        if (i == liczba_iteracji)
        {
            cout << "\nSTOP: przekroczono liczbe iteracji\n";
            break;
        }
        x0 = x1;
    }
}
void metoda_siecznych(double (*f)(double), double start, double start2)
{
    double x0, x1, x2;
    x0 = start;
    x1 = start2;
    cout << "\t x:\t\t  Estymator bledu:\t      residuum:\n";
    for (int i = 0; i < liczba_iteracji; i++)
    {
        x2 = x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0));
        estymator_bledu = fabs(x1 - x2);
        residuum = f(x2);
        cout << x2 << "\t" << estymator_bledu << "\t" << residuum << endl;
        if (fabs(estymator_bledu) <= tol_x)
        {
            cout << "\nSTOP: Estymator bledu mniejszy od zadanej tolerancji\n";
            break;
        }
        if (fabs(residuum) <= tol_f)
        {
            cout << "\nSTOP: residuum mniejsze od zadanej tolerancji\n";
            break;
        }
        if (i == liczba_iteracji)
        {
            cout << "\nSTOP: przekroczono liczbe iteracji\n";
            break;
        }
        x0 = x1;
        x1 = x2;
    }
}

int main()
{
    cout.precision(16);
    cout.setf(ios::fixed | ios::showpos);
    cout << "\nMetoda Picarda,funkcja1" << endl;
    metoda_picarda(phi_picard_1, funkcja_1);
    cout << "\nMetoda Picarda,funkcja2" << endl;
    metoda_picarda(phi_picard_2, funkcja_2);
    cout << "\nMetoda Bisekcji,funkcja1" << endl;
    metoda_bisekcji(funkcja_1, -2, 1);
    cout << "\nMetoda Bisekcji,funkcja2" << endl;
    metoda_bisekcji(funkcja_2, 0, 0.7);
    cout << "\nMetoda Newtona,funkcja1" << endl;
    metoda_newtona(funkcja_1, funkcja_1_pochodna, 1);
    cout << "\nMetoda Newtona,funkcja2" << endl;
    metoda_newtona(funkcja_2, funkcja_2_pochodna, 0.7);
    cout << "\nMetoda Siecznych,funkcja1" << endl;
    metoda_siecznych(funkcja_1, 1, 2);
    cout << "\nMetoda Siecznych,funkcja2" << endl;
    metoda_siecznych(funkcja_2, -0.5, 0);
    return 0;
}