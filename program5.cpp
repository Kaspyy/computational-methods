#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace std;

void wypisz_LU(double A[4][4], int n, int indeks[])
{
    int i, j;
    cout << "macierz U:" << endl;
    for (i = 0; i < n; i++)
    {
        cout << setw(2) << "|";
        for (j = 0; j < n; j++)
        {
            if (j < i)
                cout << setw(10) << 0;
            else
                cout << setw(10) << A[indeks[i]][j];
        }
        cout << setw(1) << "|";
        cout << endl;
    }
    cout << endl
         << endl;
    cout << "macierz L:" << endl;
    for (i = 0; i < n; i++)
    {
        cout << setw(2) << "|";
        for (j = 0; j < n; j++)
        {
            if (j > i)
                cout << setw(10) << 0;
            else
            {
                if (j == i)
                    cout << setw(10) << 1;
                else
                    cout << setw(10) << A[indeks[i]][j];
            }
        }
        cout << setw(1) << "|";
        cout << endl;
    }
}
void dek_LU(double macierz[4][4], int n, int indeks[])
{
    int l, indeks_max, i, j, k, temp;
    double max = 0;
    for (k = 0; k < n; k++)
    { // sterowanie wielkością macierzy do eliminacji Gaussa
        for (i = k + 1; i < n; i++)
        {
            for (j = n - 1; j >= k; j--)
            {
                double e = macierz[indeks[i]][k] / macierz[indeks[k]][k]; //wspolczynnik
                if (macierz[indeks[k]][k] != 0)
                { //wypelnianie macierzy U kolejnymi wyrazami
                    macierz[indeks[i]][j] = macierz[indeks[i]][j] - macierz[indeks[k]][j] * (macierz[indeks[i]][k] / macierz[indeks[k]][k]);
                    if (j == k)
                        macierz[indeks[i]][j] = e; //wypelnienie macierzy wspolczynnikiem
                }
                else
                {
                    max = macierz[indeks[k]][k];
                    l = k + 1;
                    while (l < n)
                    {
                        if (fabs(macierz[indeks[l]][k]) > max)
                        { //czesciowy wyboru elementu podstawowego
                            max = macierz[indeks[l]][k];
                            indeks_max = l;
                        }
                        l++;
                    }
                    if (indeks_max == n)
                    {
                        cout << "Brak rozwiazan, lub nieskonczenie wiele rozwiazan.\n";
                        system("pause");
                    }
                    temp = indeks[k];
                    indeks[k] = indeks[indeks_max];
                    indeks[indeks_max] = temp;
                }
            }
        }
    }
    wypisz_LU(macierz, n, indeks);
}

void rozwiazanie(double A[4][4], int n, double b[4], double x[4], double y[4], int indeks[])
{
    int i, j;
    y[0] = b[indeks[0]] / A[indeks[0]][0];
    double suma = 0;
    for (i = 1; i < n; i++)
    {
        for (j = 0; j < i; j++)
            suma = suma + A[indeks[i]][j] * y[j];
        y[i] = (b[indeks[i]] - suma);
        suma = 0;
    }
    x[n - 1] = y[n - 1] / A[indeks[n - 1]][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        for (int j = i + 1; j < n; j++)
            suma = suma + A[indeks[i]][j] * x[j];
        x[i] = (y[i] - suma) / A[indeks[i]][i];
        suma = 0;
    }
}

void wypisz_macierz(double macierz[4][4], int n, int indeks[])
{
    int i, j;
    for (i = 0; i < n; i++) //usuwanie wierszy w macierzy
    {
        cout << setw(2) << "|";
        for (j = 0; j < n; j++) //usuwanie kolumn w macierzy
        {
            cout << setw(10) << macierz[indeks[i]][j]; //setw()- ustawia odstepy
        }
        cout << "|" << endl;
    }
}
void wypisz_wektor(double wektor[4], int n, int indeks[4])
{
    int i;
    for (i = 0; i < n; i++)
        cout << setw(20) << "|" << setw(9) << wektor[indeks[i]] << "|" << endl;
}

int main()
{
    double A[4][4] = {
        {1.0, -20.0, 30.0, -4.0},
        {2.0, -40.0, -6.0, 50.0},
        {9.0, -180.0, 11.0, -12.0},
        {-16.0, 15.0, -140.0, 13.0}};
    double b[4] =
        {
            35.0,
            104.0,
            -366.0,
            -354.0};

    double y[4], x[4];
    int indeks[4];
    int n = 4, i, j;
    for (i = 0; i < n; i++) //wektor indeksów
        indeks[i] = i;
    cout << "macierz A:" << endl;
    wypisz_macierz(A, n, indeks);
    cout << endl;
    cout << "wektor B:" << endl;
    wypisz_wektor(b, n, indeks);
    cout << endl;
    dek_LU(A, n, indeks);
    rozwiazanie(A, n, b, x, y, indeks);

    cout << "\nwektor X:" << endl;
    wypisz_wektor(x, n, indeks);

    cout << "\nwektor Y:" << endl;
    wypisz_wektor(y, n, indeks);

    return 0;
}