#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;
const int iter = 100;

template <typename T>
T f(T x)
{
    return T(sin(x));
}

template <typename T>
T df(T x)
{
    return T(cos(x));
}
template <typename T> void roznicaBlad()
{
    T tab[iter][5];
    T blad_bezwgledny[iter][6];
    T h = 0.1;
    T poczatek_przedzialu = 0.0;
    T srodek_przedzialu = M_PI / 4.0;
    T koniec_przedzialu = M_PI / 2.0;
    ofstream output;
    output.open("bledy.txt");
    T err;
    for (int i = 0; i < iter; i++)
    {
        //roznica wsteczna dwupunktowa
        tab[i][0] = (f(poczatek_przedzialu + h) - f(poczatek_przedzialu)) / h;
        blad_bezwgledny[i][0] = fabs(tab[i][0] - df(poczatek_przedzialu));
        //roznica centralna dwupunktowa
        tab[i][1] = (f(srodek_przedzialu + h) - f(srodek_przedzialu - h)) / (2.0 * h);
        blad_bezwgledny[i][1] = fabs(tab[i][1] - df(srodek_przedzialu));
        //roznica progresywna dwupunktowa
        tab[i][2] = (f(koniec_przedzialu) - f(koniec_przedzialu - h)) / h;
        blad_bezwgledny[i][2] = fabs(tab[i][2] - df(koniec_przedzialu));
        //roznica progresywna trzypunktowa
        tab[i][3] = ((-3.0 / 2.0 * f(poczatek_przedzialu)) + (2.0 * f(poczatek_przedzialu + h)) - (1.0 / 2.0 * f(poczatek_przedzialu + 2 * h))) / h;
        blad_bezwgledny[i][3] = fabs(tab[i][3] - df(poczatek_przedzialu));
        //roznica wsteczna trzypunktowa
        tab[i][4] = ((3.0 / 2.0 * f(koniec_przedzialu)) - (2.0 * f(koniec_przedzialu - h)) + (1.0 / 2.0 * f(koniec_przedzialu - 2 * h))) / h;
        blad_bezwgledny[i][4] = fabs(tab[i][4] - df(koniec_przedzialu));
        blad_bezwgledny[i][5] = h;
        h *= 0.6;
        output << fixed << T(log10(blad_bezwgledny[i][5])) << " ";
        for (int j = 0; j < 5; j++)
        {
            output << fixed << T(log10(blad_bezwgledny[i][j])) << " ";
        }
        output << endl;
    }
    output.close();
    for (int j = 0; j < 5; j++)
    {
        err = (log10(blad_bezwgledny[1][j]) - log10(blad_bezwgledny[0][j])) / (log10(blad_bezwgledny[1][5]) - log10(blad_bezwgledny[0][5]));
        cout << "typ " << typeid(T).name() << " " << fixed << err << endl;
    }
}
int main()
{
    roznicaBlad<double>();
    //roznicaBlad<float>();
    return 0;
}