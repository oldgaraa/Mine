#include "stdafx.h"
#include <iostream>
#define USE_MATH_DEFINES
#include <cmath>
#include "math.h"
#include <iomanip>
#include <fstream>
using namespace std;


#define Rz 6356767                   // Условный радиус Земли
#define pc 101325                    // Нормальное давление
#define Tc 288.15                    // Нормальная температура
#define gc 9.80665                   // Нормальное ускорение свободного падения
#define k 1.4                        // Показатель адиабаты 
#define R 287.05287                  // Газовая постоянная
#define S 110.4                      // Эмперические коэффициенты
#define bettaS 1.458E-6

// Массив высот согласно варианту
const int Natm = 9;

// Массивы параметров для границ слоев атмосферы
double Hp[Natm] = { -2E3,         0E3,       11E3,          20E3,          32E3,         47E3,         51E3,         71E3,         85E3 };
double hp[Natm] = { -1.99937E3,   0.0E3,     11.01907E3,    20.06312E3,    32.16190E3,   47.35009E3,   51.41248E3,   71.80197E3,   86.15199E3 };
double Tp[Natm] = { 301.15,       288.15,    216.65,        216.65,        228.65,       270.65,       270.65,       214.65,       186.65 };
double betta_GOST[Natm] = { -6.5E-3,      -6.5E-3,   0.0E-3,        1.0E-3,        2.8E-3,       0.0E-3,       -2.8E-3,      -2.0E-3,      0.0E-3 };
double pp[Natm] = { 0 };



double fH(double h)
{
	return(Rz*h / (Rz + h));
}
double fh(double H)
{
	return (Rz*H / (Rz - H));
}
int findi(double h)
{
	double H = fH(h);
	for (int i = 1; i <= 9; i++)
	{
		if ((H > Hp[i - 1]) && (H <= Hp[i]))
		{
			return (i - 1);
			break;
		}
	}
}
double fT(double h)
{
	int i = findi(h);
	double H = fH(h);
	return (Tp[i] + betta_GOST[i]*(H - Hp[i]));
}
double fp(double h)
{
	int i = findi(h);
	double H = fH(h);
	if (betta_GOST[i] == 0) return (pp[i]*exp(-gc*(H - Hp[i]) / (R*fT(h))));
	else return(pp[i]*pow(1 + betta_GOST[i]*(H - Hp[i]) / Tp[i], -gc / (R*betta_GOST[i])));
}
double fro(double h)
{
	return (fp(h) / (R*fT(h)));
}
double fa(double h)
{
	return sqrt(k*R*fT(h));
}
double fmu(double h)
{
	return (bettaS*pow(fT(h), 3.0 / 2.0) / (fT(h) + S));
}
double ftetta(double h)
{
	return(fmu(h) / fro(h));
}
// Функцию ANM_GOST() необходимо запустить перед работой с остальными функциями
void Atm_GOST()
{
	pp[0] = pc / (pow(1 + betta_GOST[0] * (0 - Hp[0]) / Tp[0], -gc / (R*betta_GOST[0])));
	for (int i = 1; i <= Natm - 1; i++)
	{
		if (betta_GOST[i-1] == 0) pp[i]= (pp[i-1] * exp(-gc*(Hp[i] - Hp[i-1]) / (R*Tp[i])));
		else pp[i]=(pp[i-1] * pow(1 + betta_GOST[i-1] * (Hp[i] - Hp[i-1]) / Tp[i-1], -gc / (R*betta_GOST[i-1])));
	}

}

