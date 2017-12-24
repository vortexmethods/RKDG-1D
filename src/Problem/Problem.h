#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "defs.h"

class Problem
{
protected:
	//Указатель на объект с параметрами задачи
	const BaseParams* ptrprm;

public:
	//Размерность задачи (количество консервативных переменных)
	int dim;

	//Количество базисных функций
	int nshape;

	function<double(double)> shapefunc[3];
	double shapefuncnorm2[3];

	//Указатели на функции, содержащие начальные условия (U) либо (\rho, Ux, Uy, Uz, p)
	vector<std::function<double(const double)>> init; 

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Предельное значение (sd = left/right) q-й компоненты решения,
	//рассчитываемое по решению (U) и наклону (V)
	virtual inline double side_val(const vector<vector<double>>& sol, var q, side sd) const = 0;
	virtual inline vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const = 0;


	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Вычисление аналитического потока Flux по заданному решению U
	virtual void getFlux(const vector<double>& U, vector<double>& Flux) const = 0;
	virtual vector<double> getFlux(const vector<double>& U) const = 0;

	//Аналитические потоки, вычисленные по решению:
	//  - в центре ячейки (FLUX), 
	//  - на левом конце (LFLUX), 
	//  - на правом конце (RFLUX)
	// В позициях n и (n+1) - фиктивные ячейки
	vector<vector<double>> FLUX, LFLUX, RFLUX;

	//Заполнение вектора потоков для всех ячеек:
	//  - по решению в центах ячеек: conv_flux(UU)
	//  - по решению в центах ячеек, на левом и правом краях: conv_flux(UU, VV) с учетом фиктивных ячеек
//	void convFlux(const vector<vector<double>>& UU);
	void convFlux(const vector<vector<vector<double>>>& SOL);
		

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Возвращает вектор конс. переменных в точке x, вычисляемый по заданным в объекте Param начальным условиям
	virtual vector<double> initial_var(const double x) const = 0;


	virtual void EigenMatricies(const vector<vector<double>>& sol, \
		vector<vector<double>>& LL, \
		vector<vector<double>>& RR) const = 0;

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ С ПУСТОЙ РЕАЛИЗАЦИЕЙ ПО УМОЛЧАНИЮ
	//Печать в телефайл
	virtual void more_to_file(ostream& str, \
		const vector<vector<vector<double>>>& SOL, int cell) const {	};


	Problem(const BaseParams& prm, int dimension, int nshapefunctions, vector<std::function<double(const double)>> initv);
	~Problem();
};

#endif