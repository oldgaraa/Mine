// Интегрирование движение мины, последующая ее импульсная коррекция для попадания в заданную цель через расчет угла Пеленга
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include "math.h"
#include "GOST4401-81.h"
using namespace std;
ofstream txtfile2("2.txt");
ofstream txtfile("1.txt");

#define dtdt  K[0]
#define dxgdt  K[1]
#define dygdt  K[2]
#define dzgdt  K[3]
#define dVxgdt K[4]
#define dVygdt K[5]
#define dVzgdt K[6]
#define domegaxdt K[7]
#define domegaydt K[8]
#define domegazdt K[9]
#define droRGdt K[10]
#define dlamdaRGdt K[11]
#define dmuRGdt K[12]
#define dnuRGdt K[13]

#define t  par[0]
#define xg  par[1]
#define yg  par[2]
#define zg  par[3]
#define Vxg par[4]
#define Vyg par[5]
#define Vzg par[6]
#define omegax par[7]
#define omegay par[8]
#define omegaz par[9]
#define roRG par[10]
#define lamdaRG par[11]
#define muRG par[12]
#define nuRG par[13]

#define t_RK  parbuf[0]
#define xg_RK  parbuf[1]
#define yg_RK  parbuf[2]
#define zg_RK  parbuf[3]
#define Vxg_RK parbuf[4]
#define Vyg_RK parbuf[5]
#define Vzg_RK parbuf[6]
#define omegax_RK parbuf[7]
#define omegay_RK parbuf[8]
#define omegaz_RK parbuf[9]
#define roRG_RK parbuf[10]
#define lamdaRG_RK parbuf[11]
#define muRG_RK parbuf[12]
#define nuRG_RK parbuf[13]

// Объявление функций
double linInterp(double x1, double x2, double y1, double y2, double x);
double Cx(double M, double alpha);
double Cy_alpha(double M, double alpha);
double Cz_betta(double M, double betta);
double mx_omegax(double M, double alpha);
double mz_omegaz(double M, double alpha);
double mz_alpha(double M, double alpha);
double my_omegay(double M, double betta);
double my_betta(double M, double betta);
double mx(double M, double alpha);
double Cx_update(double M, double alpha);
double Cy_alpha_update(double M, double alpha);
double Cz_betta_update(double M, double betta);
double mx_omegax_update(double M, double alpha);
double mz_omegaz_update(double M, double alpha);
double mz_alpha_update(double M, double alpha);
double my_omegay_update(double M, double betta);
double my_betta_update(double M, double betta);
double mx_update(double M, double alpha);
void MATR(double ro, double lamda, double mu, double nu, double(*A)[3][3], double(*AT)[3][3]);
double fromPItoPI(double angle);

const double g = 9.80665;
double K[99], par[99], parbuf[99], par_drob[99],
A[3][3] = { 0 }, AT[3][3] = { 0 }, B[3][3] = { 0 };
double tang, risk, kren, m, V0, q, M, r, V, X, Y, Z, Mx, My, Mz, Mstab, alpha, betta, alpha_full, vx, vy, vz, vr, vfi, vhi,
delta_tang, delta_risk, delta_kren, x_t, y_t, z_t, Vx_t = 0, Vy_t = 0, Vz_t = 0, Ksi, hi, x, y, z, rg, ky = 0, kz = 0, t_on = 1E10,
Ix = 0.041,
Iy = 0.753,
Iz = 0.753,
dm = 0.003 * 11 + 0.11,
Sm = M_PI*dm*dm / 4.0,
L = 0.90,
Kvr = 0.1,
dhi = 0.03*M_PI / 180,
t_nominal = 1e100,
t_GO = t_nominal - 3.5;
bool ENGINE_WORKS = false;
int i_engine;

int main()
{
	Atm_GOST();

			//Начальные условия
		double N = 11,
			tang = 47 * M_PI / 180.0;
			risk = 0;
			kren = 0;
			m = 19;
			V0 = 440;
			int colOfEngine = 7;
			double deltahi = 2 * M_PI / colOfEngine;
			double hi_engine[100]; bool engine_done[100] = { false };
			for (int i = 0; i <= colOfEngine - 1; i++) hi_engine[i] = i*deltahi;
			double P_ud = 5000;
			double t_engine = 0.04;
			double Ksi_min = 5*M_PI / 180;
			double Ksi_max = 13 * M_PI / 180;
			double dt = 0.001;
			// номинал
			// x_t = 4851  z = -2
			x_t = 4873;
			y_t = 0;
			z_t = -7;
			t = 0;
			xg = 0;
			yg = 1e-3;
			zg = 0;
			Vxg = V0*cos(tang);
			Vyg = V0*sin(tang);
			Vzg = 0;
			omegax = 0;
			omegay = 0;
			omegaz = 0;
			roRG = cos(risk*0.5)*cos(tang*0.5)*cos(kren*0.5) - sin(risk*0.5)*sin(kren*0.5)*sin(tang*0.5);
			lamdaRG = sin(risk*0.5)*sin(tang*0.5)*cos(kren*0.5) + cos(risk*0.5)*cos(tang*0.5)*sin(kren*0.5);
			muRG = sin(risk*0.5)*cos(tang*0.5)*cos(kren*0.5) + cos(risk*0.5)*sin(tang*0.5)*sin(kren*0.5);
			nuRG = cos(risk*0.5)*sin(tang*0.5)*cos(kren*0.5) - sin(risk*0.5)*cos(tang*0.5)*sin(kren*0.5);
			
			txtfile << "t,s" << ";" << "x,m" << ";" << "y,m" << ";" << "z,m" << ";" << "Vx,m/s" << ";" << "Vy,m/s" << ";" << "Vz,m/s" << ";" << "V,m/s" << ";"
				<< "tang,grad" << ";" << "risk,grad" << ";" << "kren,grad" << ";" << "omega_x,grad/s" << ";" << "omega_y,grad/s" << ";" << "omega_z,grad/s" << ";"
				<< "alpha,grad" << ";" << "betta,grad" << ";" << "Ksi, grad" << ";" << "hi, grad" << ";" << "rg, m" << ";" << "ENGINE_WORKS" << ";" << "engine_done" << endl;
			// Интегрирование траектории до сброса ГО
			int step = 0;
			for (;;) {

				// Выход из цикла
				if (abs(t - t_GO) < 1e-5) break;
				if (abs(yg) < 1e-6) break;
				// Дробление шага
				if (t - t_GO > 0) {
				if(yg<0){
					dt *= 0.5;
					for (int i = 0; i <= 13; i++) par[i] = par_drob[i];
				}
				else {
					for (int i = 0; i <= 13; i++) par_drob[i] = par[i];

					// Расчет промежуточных параметров на шаге
					MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
					tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
					risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
					kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
					r = pow((x_t - xg)*(x_t - xg) + (y_t - yg)*(y_t - yg) + (z_t - zg)*(z_t - zg), 0.5),
						V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
						q = fro(yg)*V*V*0.5,
						vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
						vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
						vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
						alpha = -atan2(vy, vx);
					betta = asin(vz / V);
					alpha_full = sqrt(alpha*alpha + betta*betta),
					M = V / fa(yg);
					X = -Cx(M, alpha_full)*q*Sm;
					Y = (Cy_alpha(M, alpha)*alpha)*q*Sm;
					Z = (Cz_betta(M, betta)*betta)*q*Sm;
					Mx = (mx_omegax(M, alpha_full)*omegax*dm / V + Kvr*mx(M, alpha_full))*q*Sm*dm;
					My = (my_omegay(M, betta)*omegay*dm / V + my_betta(M, betta)*betta)*q*Sm*L;
					Mz = (mz_omegaz(M, alpha)*omegaz*dm / V + mz_alpha(M, alpha)*alpha)*q*Sm*L;
					cout << yg << endl;

					//вывод в файл
					//int st = 5 / dt;
					int st = 0.1 / dt;
					if (step % st == 0)
						txtfile
						<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
						<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
						<< omegax*180/M_PI  << ";" << omegay * 180 / M_PI << ";" << omegaz * 180 / M_PI << ";"
						<< alpha *180/M_PI << ";" << betta * 180 / M_PI << ";"
						<< endl;
				}
				step++;

					//Производные
					dtdt = 1 * dt,
						dxgdt = Vxg*dt,
						dygdt = Vyg*dt,
						dzgdt = Vzg*dt,
						dVxgdt = 1 / m*(X*A[0][0] + Y*A[0][1] + Z*A[0][2])*dt,
						dVygdt = 1 / m*(X*A[1][0] + Y*A[1][1] + Z*A[1][2] - m*g)*dt,
						dVzgdt = 1 / m*(X*A[2][0] + Y*A[2][1] + Z*A[2][2])*dt,
						domegaxdt = (Mx / Ix - omegay*omegaz*(Iz - Iy) / Ix)*dt,
						domegaydt = (My / Iy - omegax*omegaz*(Ix - Iz) / Iy)*dt,
						domegazdt = (Mz / Iz - omegay*omegax*(Iy - Ix) / Iz)*dt,
						droRGdt = -(omegax*lamdaRG + omegay*muRG + omegaz*nuRG)*0.5*dt,
						dlamdaRGdt = (omegax*roRG - omegay*nuRG + omegaz*muRG)*0.5*dt,
						dmuRGdt = (omegax*nuRG + omegay*roRG - omegaz*lamdaRG)*0.5*dt,
						dnuRGdt = (-omegax*muRG + omegay*lamdaRG + omegaz*roRG)*0.5*dt;

				// Приращения
				for (int i = 0; i <= 13; i++) par[i] += K[i];

				// Нормировка параметров РГ
				double norm_RG = pow(roRG*roRG + lamdaRG*lamdaRG + nuRG*nuRG + muRG*muRG, 0.5);
				roRG /= norm_RG;
				lamdaRG /= norm_RG;
				nuRG /= norm_RG;
				muRG /= norm_RG;
			}  
			// Конец интегрирования траектории до сброса ГО

			m = 19;
			L = 0.78;
			dt = 0.0001;
			// ПОСЛЕ СБРОСА ГО
			for (;;) {
				// Выход из цикла
				if (yg < 1E-5) break;
				// Дробление шага
				if (yg < 0) {
					dt *= 0.5;
					for (int i = 0; i <= 13; i++) par[i] = par_drob[i];
				}
				else {
					for (int i = 0; i <= 13; i++) par_drob[i] = par[i];

					//Выключение двигателя
					if ((ENGINE_WORKS == true) && (t >= t_on + t_engine)) {
						ENGINE_WORKS = false;
						ky = 0;
						kz = 0;
						t_on = 1E10;
					}
					//Пересчет k1 k2
					if (ENGINE_WORKS == true) {
						ky = cos(hi_engine[i_engine] + kren);
						kz = sin(hi_engine[i_engine] + kren);
					}

					// Расчет промежуточных параметров на шаге
					MATR(roRG, lamdaRG, muRG, nuRG, &A, &AT);
					tang = asin(2 * (roRG*nuRG + lamdaRG*muRG));
					risk = atan2(2 * (roRG*muRG - lamdaRG*nuRG), roRG*roRG + lamdaRG*lamdaRG - nuRG*nuRG - muRG*muRG);
					kren = atan2(2 * (roRG*lamdaRG - muRG*nuRG), roRG*roRG - lamdaRG*lamdaRG - nuRG*nuRG + muRG*muRG);
					V = pow(Vxg*Vxg + Vyg*Vyg + Vzg*Vzg, 0.5),
						q = fro(yg)*V*V*0.5,
						vx = Vxg*AT[0][0] + Vyg*AT[0][1] + Vzg*AT[0][2],
						vy = Vxg*AT[1][0] + Vyg*AT[1][1] + Vzg*AT[1][2],
						vz = Vxg*AT[2][0] + Vyg*AT[2][1] + Vzg*AT[2][2],
						alpha = -atan2(vy, vx);
					betta = asin(vz / V);
					alpha_full = sqrt(alpha*alpha + betta*betta);
					M = V / fa(yg);
					X = -Cx_update(M, alpha_full)*q*Sm;
					Y = (Cy_alpha_update(M, alpha)*alpha)*q*Sm + ky*P_ud*g;
					Z = (Cz_betta_update(M, betta)*betta)*q*Sm + kz*P_ud*g;
					Mx = (mx_omegax_update(M, alpha_full)*omegax*dm / V + Kvr*mx_update(M, alpha_full))*q*Sm*dm;
					My = (my_omegay_update(M, betta)*omegay*dm / V + my_betta_update(M, betta)*betta)*q*Sm*L;
					Mz = (mz_omegaz_update(M, alpha)*omegaz*dm / V + mz_alpha_update(M, alpha)*alpha)*q*Sm*L;
					cout << yg << endl;

					// КОРРЕКТИРОВКА***************************************
					rg = pow((x_t - xg)*(x_t - xg) + (y_t - yg)*(y_t - yg) + (z_t - zg)*(z_t - zg), 0.5);
					x = (x_t - xg)*AT[0][0] + (y_t - yg)*AT[0][1] + (z_t - zg)*AT[0][2];
					y = (x_t - xg)*AT[1][0] + (y_t - yg)*AT[1][1] + (z_t - zg)*AT[1][2];
					z = (x_t - xg)*AT[2][0] + (y_t - yg)*AT[2][1] + (z_t - zg)*AT[2][2];
					r = sqrt(y*y + z*z);
					Ksi = atan2(r, x);
					Ksi = fromPItoPI(Ksi);
					hi = atan2(z, y);
					hi = fromPItoPI(hi);
					if (rg > 50)
						if ((ENGINE_WORKS == false) && (t >= t_nominal - 3))
							if ((Ksi < Ksi_max) && (Ksi > Ksi_min))
								for (int i = 0; i <= colOfEngine-1; i++) 
									if ((abs(fromPItoPI(hi - omegax*t_engine*0.5) - fromPItoPI(hi_engine[i] + kren)) <= dhi) && (engine_done[i] == false)) {
										ENGINE_WORKS = true;
										engine_done[i] = true;
										i_engine = i;
										ky = cos(hi_engine[i_engine] + kren);
										kz = sin(hi_engine[i_engine] + kren);
										t_on = t;
										break;
									}

					//вывод в файл
					//int st = 0.2 / dt;
					int st = 0.01 / dt;
					if (step % st == 0) {
						txtfile
							<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
							<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
							<< omegax * 180 / M_PI << ";" << omegay * 180 / M_PI << ";" << omegaz * 180 / M_PI << ";"
							<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
							<< Ksi * 180 / M_PI << ";" << hi * 180 / M_PI << ";" << rg << ";" << ENGINE_WORKS << ";";
						for (int i = 0; i < 8; i++) txtfile << engine_done[i] << ";";
						txtfile << endl;
					}
				}
				step++;

				//Производные
				dtdt = 1 * dt,
					dxgdt = Vxg*dt,
					dygdt = Vyg*dt,
					dzgdt = Vzg*dt,
					dVxgdt = 1 / m*(X*A[0][0] + Y*A[0][1] + Z*A[0][2])*dt,
					dVygdt = 1 / m*(X*A[1][0] + Y*A[1][1] + Z*A[1][2] - m*g)*dt,
					dVzgdt = 1 / m*(X*A[2][0] + Y*A[2][1] + Z*A[2][2])*dt,
					domegaxdt = (Mx / Ix - omegay*omegaz*(Iz - Iy) / Ix)*dt,
					domegaydt = (My / Iy - omegax*omegaz*(Ix - Iz) / Iy)*dt,
					domegazdt = (Mz / Iz - omegay*omegax*(Iy - Ix) / Iz)*dt,
					droRGdt = -(omegax*lamdaRG + omegay*muRG + omegaz*nuRG)*0.5*dt,
					dlamdaRGdt = (omegax*roRG - omegay*nuRG + omegaz*muRG)*0.5*dt,
					dmuRGdt = (omegax*nuRG + omegay*roRG - omegaz*lamdaRG)*0.5*dt,
					dnuRGdt = (-omegax*muRG + omegay*lamdaRG + omegaz*roRG)*0.5*dt;

				// Приращения
				for (int i = 0; i <= 13; i++) par[i] += K[i];

				// Нормировка параметров РГ
				double norm_RG = pow(roRG*roRG + lamdaRG*lamdaRG + nuRG*nuRG + muRG*muRG, 0.5);
				roRG /= norm_RG;
				lamdaRG /= norm_RG;
				nuRG /= norm_RG;
				muRG /= norm_RG;
			}
			// КОНЕЦ ИНТЕГРИРОВАНИЯ

			   // Вывод последнего шага
			txtfile
				<< t << ";" << xg << ";" << yg << ";" << zg << ";" << Vxg << ";" << Vyg << ";" << Vzg << ";" << V << ";"
				<< tang * 180 / M_PI << ";" << risk * 180 / M_PI << ";" << kren * 180 / M_PI << ";"
				<< omegax * 180 / M_PI << ";" << omegay * 180 / M_PI << ";" << omegaz * 180 / M_PI << ";"
				<< alpha * 180 / M_PI << ";" << betta * 180 / M_PI << ";"
				<< Ksi * 180 / M_PI << ";" << hi * 180 / M_PI << ";" << rg << ";" << ENGINE_WORKS << ";";
			for (int i = 0; i < 8; i++) txtfile << engine_done[i] << ";";
			txtfile << endl;

			cout << rg << endl;
	cout << "END!" << endl;
	system("pause");
	return 0;
}

//Матрицы
void MATR(double ro, double lamda, double mu, double nu, double(*A)[3][3], double(*AT)[3][3]) {
	double a[3][3] = { 0 };
	// Матрица из нормальной земной в связную
	a[0][0] = ro*ro + lamda*lamda - mu*mu - nu*nu;
	a[0][1] = 2 * (-ro*nu + lamda*mu);
	a[0][2] = 2 * (ro*mu + lamda*nu);

	a[1][0] = 2 * (ro*nu + lamda*mu);
	a[1][1] = ro*ro - lamda*lamda + mu*mu - nu*nu;
	a[1][2] = 2 * (-ro*lamda + nu*mu);

	a[2][0] = 2 * (-ro*mu + lamda*nu);
	a[2][1] = 2 * (ro*lamda + nu*mu);
	a[2][2] = ro*ro - lamda*lamda - mu*mu + nu*nu;

	for (int i = 0; i <= 2; i++) for (int j = 0; j <= 2; j++)
		(*A)[i][j] = a[i][j];
	//Из связной в нормальную
	for (int i = 0; i <= 2; i++) for (int j = 0; j <= 2; j++)
		(*AT)[i][j] = a[j][i];
}

double fromPItoPI(double angle) {
	if (angle < -M_PI) angle += 2 * M_PI;
	if (angle > M_PI) angle -= 2 * M_PI;
	return angle;
}

double linInterp(double x1, double x2, double y1, double y2, double x) {
double y = (y2*(x - x1) + y1 * (x2 - x)) / (x2 - x1);
return(y);
}

double Cx(double mach, double alpha) {
	double Cx[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.40 },
		{ 1.000, 0.325, 0.325, 0.331, 0.338, 0.368, 0.442, 0.619, 0.646, 0.701, 0.650, 0.628 },
		{ 2.000, 0.326, 0.326, 0.334, 0.341, 0.378, 0.458, 0.623, 0.652, 0.710, 0.662, 0.632 },
		{ 3.000, 0.334, 0.334, 0.338, 0.345, 0.382, 0.462, 0.632, 0.663, 0.722, 0.673, 0.641 },
		{ 4.000, 0.349, 0.349, 0.355, 0.362, 0.405, 0.478, 0.661, 0.675, 0.730, 0.682, 0.650 },
		{ 5.000, 0.396, 0.396, 0.402, 0.410, 0.441, 0.493, 0.682, 0.702, 0.754, 0.695, 0.664 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > Cx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > Cx[i][0]) && (i<5); i++);


	double cx1 = linInterp(Cx[i - 1][0], Cx[i][0], Cx[i - 1][j - 1], Cx[i][j - 1], alphaToDeg);
	double cx2 = linInterp(Cx[i - 1][0], Cx[i][0], Cx[i - 1][j], Cx[i][j], alphaToDeg);

	double cx = linInterp(Cx[0][j - 1], Cx[0][j], cx1, cx2, mach);
	//std::cout « i « '\t' « j « '\t' « cx « '\n';
	if (cx < 0.2) return 0.2;
	if (cx > 0.8) return 0.8;
	return cx;
}

double Cy_alpha(double mach, double alpha) {
	double Cy[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, 0.156, 0.156, 0.181, 0.188, 0.192, 0.214, 0.226, 0.218, 0.206, 0.198, 0.187 },
		{ 2.000, 0.155, 0.155, 0.180, 0.187, 0.191, 0.213, 0.225, 0.217, 0.205, 0.197, 0.186 },
		{ 3.000, 0.155, 0.155, 0.180, 0.187, 0.191, 0.213, 0.225, 0.216, 0.205, 0.197, 0.186 },
		{ 4.000, 0.154, 0.154, 0.179, 0.186, 0.190, 0.212, 0.224, 0.216, 0.204, 0.197, 0.186 },
		{ 5.000, 0.154, 0.154, 0.179, 0.186, 0.190, 0.212, 0.224, 0.215, 0.204, 0.196, 0.185 }

	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > Cy[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > Cy[i][0]) && (i < 5); i++);
	double cy1 = linInterp(Cy[i - 1][0], Cy[i][0], Cy[i - 1][j - 1], Cy[i][j - 1], alphaToDeg);
	double cy2 = linInterp(Cy[i - 1][0], Cy[i][0], Cy[i - 1][j], Cy[i][j], alphaToDeg);
	double cy = linInterp(Cy[0][j - 1], Cy[0][j], cy1, cy2, mach);
	// std::cout « i « '\t' « j « '\t' « cy « '\n';
	return cy;
}

double mz_alpha(double mach, double alpha) {
	double mz[6][12] = {
		{ 0.000, 0.00000, 0.50000, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -0.0218, -0.0218, -0.0286, -0.0307, -0.0319, -0.0326, -0.0377, -0.0368, -0.0331, -0.0289, -0.0260 },
		{ 2.000, -0.0217, -0.0217, -0.0285, -0.0306, -0.0318, -0.0325, -0.0376, -0.0367, -0.0330, -0.0288, -0.0259 },
		{ 3.000, -0.0217, -0.0217, -0.0285, -0.0306, -0.0318, -0.0325, -0.0376, -0.0366, -0.0330, -0.0288, -0.0259 },
		{ 4.000, -0.0216, -0.0216, -0.0284, -0.0305, -0.0317, -0.0324, -0.0375, -0.0366, -0.0329, -0.0288, -0.0259 },
		{ 5.000, -0.0216, -0.0216, -0.0284, -0.0305, -0.0317, -0.0324, -0.0375, -0.0365, -0.0329, -0.0287, -0.0258 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mz[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mz[i][0]) && (i < 5); i++);
	double mz1 = linInterp(mz[i - 1][0], mz[i][0], mz[i - 1][j - 1], mz[i][j - 1], alphaToDeg);
	double mz2 = linInterp(mz[i - 1][0], mz[i][0], mz[i - 1][j], mz[i][j], alphaToDeg);
	double mz_res = linInterp(mz[0][j - 1], mz[0][j], mz1, mz2, mach);
	// std::cout « i « '\t' « j « '\t' « mz_res « '\n';
	return mz_res;
}

double mz_omegaz(double mach, double alpha) {
	double mzWz[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -6.317, -6.317, -7.176, -7.434, -7.520, -8.337, -9.067, -9.368, -9.076, -8.938, -8.509 },
		{ 2.000, -6.309, -6.309, -7.167, -7.424, -7.511, -8.326, -9.055, -9.356, -9.065, -8.927, -8.498 },
		{ 3.000, -6.290, -6.290, -7.156, -7.413, -7.498, -8.314, -9.042, -9.342, -9.051, -8.913, -8.485 },
		{ 4.000, -6.250, -6.250, -7.141, -7.397, -7.483, -8.296, -9.022, -9.322, -9.031, -8.849, -8.467 },
		{ 5.000, -6.241, -6.241,

		-7.104, -7.359, -7.444, -8.253, -8.976, -9.274, -8.969, -8.848, -8.423 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mzWz[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mzWz[i][0]) && (i < 5); i++);
	double mz1 = linInterp(mzWz[i - 1][0], mzWz[i][0], mzWz[i - 1][j - 1], mzWz[i][j - 1], alphaToDeg);
	double mz2 = linInterp(mzWz[i - 1][0], mzWz[i][0], mzWz[i - 1][j], mzWz[i][j], alphaToDeg);
	double mz = linInterp(mzWz[0][j - 1], mzWz[0][j], mz1, mz2, mach);
	//std::cout « i « '\t' « j « '\t' « "mzWz" « '\n';
	return mz;
}

double mx(double mach, double alpha) {
	double mx[6][12] = {
		{ 0.000, 0.0000, 0.5000, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -0.203, -0.203, -0.244, -0.256, -0.261, -0.294, -0.310, -0.289, -0.261, -0.244, -0.223 },
		{ 2.000, -0.203, -0.203, -0.243, -0.255, -0.261, -0.293, -0.309, -0.289, -0.260, -0.243, -0.223 },
		{ 3.000, -0.202, -0.202, -0.243, -0.255, -0.260, -0.293, -0.309, -0.288, -0.260, -0.243, -0.222 },
		{ 4.000, -0.202, -0.202, -0.242, -0.254, -0.260, -0.292, -0.308, -0.288, -0.260, -0.242, -0.222 },
		{ 5.000, -0.202, -0.202, -0.242, -0.254, -0.259, -0.292, -0.308, -0.287, -0.259, -0.242, -0.221 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mx[i][0]) && (i < 5); i++);
	double mx1 = linInterp(mx[i - 1][0], mx[i][0], mx[i - 1][j - 1], mx[i][j - 1], alphaToDeg);
	double mx2 = linInterp(mx[i - 1][0], mx[i][0], mx[i - 1][j], mx[i][j], alphaToDeg);
	double mx_res = linInterp(mx[0][j - 1], mx[0][j], mx1, mx2, mach);
	// std::cout « i « '\t' « j « '\t' « alphaToDeg « '\t' « mx_res « '\n';
	return mx_res;
}

double mx_omegax(double mach, double alpha) {
	double mxWx[6][12] = {
		{ 0.000, 0.00000, 0.50000, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -11.288, -11.288, -13.596, -14.282, -14.514, -16.358, -16.276, -16.127, -14.514, -13.596, -12.438 },
		{ 2.000, -11.274, -11.274, -13.579, -14.261, -14.496, -16.338, -16.256, -16.107, -14.496, -13.579, -12.422 },
		{ 3.000, -11.257, -11.257, -13.557, -14.242, -14.474, -16.313, -16.231, -16.082, -14.474, -13.558, -12.403 },
		{ 4.000, -11.233, -11.233, -13.528, -14.212, -14.443, -16.278, -16.196, -16.048, -14.443, -13.529, -12.377 },
		{ 5.000, -11.175, -11.175, -13.460, -14.139, -14.368, -16.194, -16.113, -16.065, -14.368, -13.460, -12.313 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mxWx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mxWx[i][0]) && (i < 5); i++);
	double mx1 = linInterp(mxWx[i - 1][0], mxWx[i][0], mxWx[i - 1][j - 1], mxWx[i][j - 1], alphaToDeg);
	double mx2 = linInterp(mxWx[i - 1][0], mxWx[i][0], mxWx[i - 1][j], mxWx[i][j], alphaToDeg);
	double mx = linInterp(mxWx[0][j - 1], mxWx[0][j], mx1, mx2, mach);
	//std::cout « i « '\t' « j « '\t' « "mxWx" « '\n';
	return mx;
}

double Cz_betta(double M, double betta) {
	return -Cy_alpha(M, betta);
}

double my_omegay(double M, double betta) {
	return mz_omegaz(M, betta);
}

double my_betta(double M, double betta) {
	return mz_alpha(M, betta);
}

//////////////////////////////////////ПОСЛЕ СБРОСА ГО/////////////////////////////////////////////////////////////

double Cx_update(double mach, double alpha) {
	double Cx[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, 0.334, 0.334, 0.336, 0.341, 0.418, 0.520, 0.787, 0.818, 0.901, 0.893, 0.884 },
		{ 2.000, 0.335, 0.335, 0.337, 0.342, 0.420, 0.522, 0.790, 0.821, 0.905, 0.896, 0.887 },
		{ 3.000, 0.337, 0.337, 0.339, 0.344, 0.421, 0.524, 0.793, 0.824, 0.908, 0.900, 0.891 },
		{ 4.000, 0.338, 0.338, 0.340, 0.345, 0.423, 0.526, 0.796, 0.828, 0.911, 0.903, 0.894 },
		{ 5.000, 0.339, 0.339, 0.341, 0.346, 0.425, 0.528, 0.799, 0.831, 0.915, 0.907, 0.898 }

	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > Cx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > Cx[i][0]) && (i < 5); i++);
	double cx1 = linInterp(Cx[i - 1][0], Cx[i][0], Cx[i - 1][j - 1], Cx[i][j - 1], alphaToDeg);
	double cx2 = linInterp(Cx[i - 1][0], Cx[i][0], Cx[i - 1][j], Cx[i][j], alphaToDeg);
	double cx = linInterp(Cx[0][j - 1], Cx[0][j], cx1, cx2, mach);
	if (cx < 0.2) return 0.2;
	if (cx > 0.8) return 0.8;
	return cx;
}

double Cy_alpha_update(double mach, double alpha) {
	double Cy[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, 0.142, 0.142, 0.168, 0.174, 0.177, 0.198, 0.209, 0.200, 0.185, 0.177, 0.165 },
		{ 2.000, 0.142, 0.142, 0.167, 0.174, 0.176, 0.197, 0.208, 0.199, 0.184, 0.176, 0.165 },
		{ 3.000, 0.142, 0.142, 0.167, 0.173, 0.176, 0.197, 0.208, 0.199, 0.184, 0.176, 0.164 },
		{ 4.000, 0.141, 0.141, 0.167, 0.173, 0.175, 0.197, 0.207, 0.198, 0.183, 0.175, 0.164 },
		{ 5.000, 0.141, 0.141, 0.166, 0.172, 0.175, 0.196, 0.206, 0.198, 0.183, 0.175, 0.163 }

	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > Cy[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > Cy[i][0]) && (i < 5); i++);
	double cy1 = linInterp(Cy[i - 1][0], Cy[i][0], Cy[i - 1][j - 1], Cy[i][j - 1], alphaToDeg);
	double cy2 = linInterp(Cy[i - 1][0], Cy[i][0], Cy[i - 1][j], Cy[i][j], alphaToDeg);
	double cy = linInterp(Cy[0][j - 1], Cy[0][j], cy1, cy2, mach);
	return cy;
}

double mz_alpha_update(double mach, double alpha) {
	double mz[6][12] = {
		{ 0.000, 0.00000, 0.50000, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -0.0291, -0.0291, -0.0371, -0.0390, -0.0394, -0.0450, -0.0489, -0.0466, -0.0405, -0.0367, -0.0340 },
		{ 2.000, -0.0291, -0.0291, -0.0370, -0.0389, -0.0394, -0.0450, -0.0488, -0.0465, -0.0404, -0.0366, -0.0340 },
		{ 3.000, -0.0290, -0.0290, -0.0370, -0.0389, -0.0394, -0.0449, -0.0487, -0.0464, -0.0404, -0.0365, -0.0339 },
		{ 4.000, -0.0290, -0.0290, -0.0370, -0.0388, -0.0393, -0.0448, -0.0487, -0.0464, -0.0403, -0.0365, -0.0339 },
		{ 5.000, -0.0290, -0.0290, -0.0369, -0.0388, -0.0392, -0.0448, -0.0486, -0.0463, -0.0403, -0.0364, -0.0338 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mz[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mz[i][0]) && (i < 5); i++);
	double mz1 = linInterp(mz[i - 1][0], mz[i][0], mz[i - 1][j - 1], mz[i][j - 1], alphaToDeg);
	double mz2 = linInterp(mz[i - 1][0], mz[i][0], mz[i - 1][j], mz[i][j], alphaToDeg);
	double mz_res = linInterp(mz[0][j - 1], mz[0][j], mz1, mz2, mach);
	return mz_res;
}

double mz_omegaz_update(double mach, double alpha) {
	double mzWz[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -5.475, -5.475, -6.368, -6.592, -6.629, -7.486, -8.119, -8.231, -7.632, -7.225, -6.704 },
		{ 2.000, -5.392, -5.392, -6.245, -6.467, -6.499, -7.340, -7.960, -8.070, -7.483, -7.084, -6.573 },
		{ 3.000, -5.338, -5.338, -6.224, -6.427, -6.463, -7.298, -7.916, -8.025, -7.441, -7.044, -6.536 },
		{ 4.000, -5.283, -5.283, -6.145, -6.352, -6.397, -7.224, -7.835, -7.943, -7.365, -6.972, -6.469 },
		{ 5.000, -5.201, -5.201, -6.065, -6.264, -6.297, -7.117, -7.713, -7.819, -7.250, -6.864, -6.389 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mzWz[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mzWz[i][0]) && (i < 5); i++);
	double mz1 = linInterp(mzWz[i - 1][0], mzWz[i][0], mzWz[i - 1][j - 1], mzWz[i][j - 1], alphaToDeg);
	double mz2 = linInterp(mzWz[i - 1][0], mzWz[i][0], mzWz[i - 1][j], mzWz[i][j], alphaToDeg);
	double mz = linInterp(mzWz[0][j - 1], mzWz[0][j], mz1, mz2, mach);
	return mz;
}

double mx_update(double mach, double alpha) {
	double mx[6][12] = {
		{ 0.000, 0.000, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -0.203, -0.203, -0.243, -0.256, -0.261, -0.294, -0.310, -0.289, -0.261, -0.244, -0.223 },
		{ 2.000, -0.203, -0.203, -0.243, -0.255, -0.261, -0.293, -0.309, -0.289, -0.260, -0.243, -0.223 },
		{ 3.000, -0.202, -0.202, -0.243, -0.255, -0.260, -0.293, -0.309, -0.288, -0.260, -0.243, -0.222 },
		{ 4.000, -0.202, -0.202, -0.242, -0.254, -0.260, -0.292, -0.308, -0.288, -0.260, -0.242, -0.222 },
		{ 5.000, -0.202, -0.202, -0.242, -0.254, -0.259, -0.292, -0.308, -0.287, -0.259, -0.242, -0.221 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mx[i][0]) && (i < 5); i++);
	double mx1 =

		linInterp(mx[i - 1][0], mx[i][0], mx[i - 1][j - 1], mx[i][j - 1], alphaToDeg);
	double mx2 = linInterp(mx[i - 1][0], mx[i][0], mx[i - 1][j], mx[i][j], alphaToDeg);
	double mx_res = linInterp(mx[0][j - 1], mx[0][j], mx1, mx2, mach);
	return mx_res;
}

double mx_omegax_update(double mach, double alpha) {
	double mxWx[6][12] = {
		{ 0.000, 0.00000, 0.50000, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400 },
		{ 1.000, -12.352, -12.352, -14.776, -15.688, -15.829, -17.638, -17.481, -17.346, -15.552, -14.517, -13.304 },
		{ 2.000, -12.111, -12.111, -14.487, -15.382, -15.520, -17.294, -17.140, -17.007, -15.248, -14.234, -13.044 },
		{ 3.000, -12.043, -12.043, -14.406, -15.296, -15.433, -17.197, -17.043, -16.912, -15.163, -14.154, -12.971 },
		{ 4.000, -11.919, -11.919, -14.259, -15.138, -15.275, -17.021, -16.859, -16.739, -15.007, -14.009, -12.838 },
		{ 5.000, -11.734, -11.734, -14.037, -15.003, -15.037, -16.754, -16.607, -16.479, -14.774, -13.791, -12.639 }
	};
	double alphaToDeg = abs(alpha) * 180/M_PI;
	int j;
	for (j = 2; mach > mxWx[0][j]; j++);
	int i;
	for (i = 2; (alphaToDeg > mxWx[i][0]) && (i < 5); i++);
	double mx1 = linInterp(mxWx[i - 1][0], mxWx[i][0], mxWx[i - 1][j - 1], mxWx[i][j - 1], alphaToDeg);
	double mx2 = linInterp(mxWx[i - 1][0], mxWx[i][0], mxWx[i - 1][j], mxWx[i][j], alphaToDeg);
	double mx = linInterp(mxWx[0][j - 1], mxWx[0][j], mx1, mx2, mach);
	return mx;
}

double Cz_betta_update(double M, double betta) {
	return Cy_alpha(M, betta);
}

double my_omegay_update(double M, double betta) {
	return mz_omegaz(M, betta);
}

double my_betta_update(double M, double betta) {
	return mz_alpha(M, betta);
}
