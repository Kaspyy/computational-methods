#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

double tol_x = 0.000000001; //zadana tolerancja błędu
double tol_f = 0.000000001; //zadana tolerancja residuum
int iter = 90;
const int n = 4;

double licz_residuum(double A[n][n], double x[n], double b[n])
{
    double residuum = 0.0;
    double max_residuum = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            residuum += A[i][j] * x[j];
        residuum = fabs(b[i] - residuum);
        if (residuum > max_residuum)
            max_residuum = residuum;
        residuum = 0.0;
    }
    return max_residuum;
}
double licz_estymator(double x0[n], double x1[n])
{
    double estymator = 0.0;
    double max_estymator = 0.0;
    for (int i = 0; i < n; i++)
    {
        estymator = fabs(x1[i] - x0[i]);
        if (estymator > max_estymator)
            max_estymator = estymator;
        estymator = 0.0;
    }
    return max_estymator;
}

void metoda_Jacobiego(double A[n][n], double b[n], double x0[n])
{

    cout << "-----------------------------------------------Metoda Jacobiego-------------------------------------------\n"
         << endl;
    int i, j, k;

    double x_next[n];
    double LiUsuma;
    for (k = 0; k < iter; k++)
    {
        for (i = 0; i < n; i++)
        {
            LiUsuma = 0.0;
            for (j = 0; j < n; j++)
                if (j != i)
                    LiUsuma += A[i][j] * x0[j];
            x_next[i] = 1 / A[i][i] * (b[i] - LiUsuma);
        }
        double estymator_bledu = licz_estymator(x0, x_next);
        double residuum = licz_residuum(A, x_next, b);

        for (i = 0; i < n; i++)
            x0[i] = x_next[i];
        cout << setw(10) << "iteracja" << setw(16) << "x0" << setw(24) << "x1" << setw(24) << "x2" << setw(24) << "x3" << setw(31) << "estymator bledu" << setw(17) << "residuum" << endl;
        cout << setw(10) << k + 1 << setw(24) << x_next[0] << setw(24) << x_next[1] << setw(24) << x_next[2] << setw(24) << x_next[3] << setw(24) << estymator_bledu << setw(24) << residuum << endl;

        if (estymator_bledu <= tol_x && residuum <= tol_f)
        {
            cout << "\nSTOP: Estymator bledu i residuum mniejsze od zadanej tolerancji\n";
            break;
        }
    }
}

void metoda_Gaussa_Seidla(double** A, double* b, double* x0)
{

    cout << "-----------------------------------------------Metoda Gaussa-Seidla-------------------------------------------\n"
         << endl;

    double x_prev[n];
    double Usuma;
    double Lsuma;

    for (int k = 0; k < iter; k++)
    {
        for (int i = 0; i < n; i++)
        {
            Usuma = 0.0;
            for (int j = i + 1; j < n; j++)
                Usuma += A[i][j] * x0[j]; //Ux
            Lsuma = 0.0;
            for (int j = 0; j < i; j++)
                Lsuma += A[i][j] * x0[j]; //Lx
            x_prev[i] = x0[i];
            x0[i] = 1 / A[i][i] * (b[i] - Usuma - Lsuma);
        }
        double estymator_bledu = licz_estymator(x0, x_prev);
        double residuum = licz_residuum(A, x0, b);
        cout << setw(10) << "iteracja" << setw(16) << "x0" << setw(24) << "x1" << setw(24) << "x2" << setw(24) << "x3" << setw(31) << "estymator bledu" << setw(17) << "residuum" << endl;
        cout << setw(10) << k + 1 << setw(24) << x0[0] << setw(24) << x0[1] << setw(24) << x0[2] << setw(24) << x0[3] << setw(24) << estymator_bledu << setw(24) << residuum << endl;

        if (estymator_bledu <= tol_x && residuum <= tol_f)
        {
            cout << "\nSTOP: Estymator bledu i residuum mniejsze od zadanej tolerancji\n";
            break;
        }
    }
}

void metoda_SOR(double A[n][n], double b[n], double x0[n])
{
    cout << "-----------------------------------------------Metoda SOR-------------------------------------------\n"
         << endl;

    double x_next[n];
    double x_prev[n];
    double Usuma;
    double Lsuma;
    double omega = 0.5;

    for (int k = 0; k < iter; k++)
    {
        for (int i = 0; i < n; i++)
        {
            Usuma = 0.0;
            for (int j = i + 1; j < n; j++)
                Usuma += A[i][j] * x0[j]; //Ux
            Lsuma = 0.0;
            for (int j = 0; j < i; j++)
                Lsuma += A[i][j] * x0[j]; //Lx
            x_prev[i] = x0[i];
            x_next[i] = (1 - omega) * x0[i] + omega / A[i][i] * (b[i] - Usuma - Lsuma);
            x0[i] = x_next[i];
        }
        double estymator_bledu = licz_estymator(x0, x_prev);
        double residuum = licz_residuum(A, x0, b);
        cout << setw(10) << "iteracja" << setw(16) << "x0" << setw(24) << "x1" << setw(24) << "x2" << setw(24) << "x3" << setw(31) << "estymator bledu" << setw(17) << "residuum" << endl;
        cout << setw(10) << k + 1 << setw(24) << x0[0] << setw(24) << x0[1] << setw(24) << x0[2] << setw(24) << x0[3] << setw(24) << estymator_bledu << setw(24) << residuum << endl;

        if (estymator_bledu <= tol_x && residuum <= tol_f)
        {
            cout << "\nSTOP: Estymator bledu i residuum mniejsze od zadanej tolerancji\n";
            break;
        }
    }
}
int main()
{
    cout.precision(14);
    cout.setf(ios::fixed);
    double A[n][n] =
        {
            100., -1., 2., -3.,
            1., 200., -4., 5.,
            -2., 4., 300., -6.,
            3., -5., 6., 400.};
    double b[n] = {
        116.,
        -226.,
        912.,
        -1174.};
    double x0[n] = {
        2.,
        2.,
        2.,
        2.};

    //metoda_Jacobiego(A, b, x0);
    //metoda_Gaussa_Seidla(A, b, x0);
    metoda_SOR(A, b, x0);
    return 0;
}