#ifndef FLUXGODUNOVTYPE_H_
#define FLUXGODUNOVTYPE_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "defs.h"
#include "Flux.h"
#include "Problem.h"
#include "ProblemTransfer1D.h"


using namespace std;

#define GodunovType FluxGodunovType


//Нелинейные функции-ограничители при расчете численного потока:
//  - поток С.К. Годунова (против потока, монотонна, 1-й порядок)
//  - поток, рассчитанный по полусумме решений в левой и правой ячейках (немонотонна, 2-й порядок)
//  - поток ван-Лира (монотонна, 2-й порядок)

#define Godunov phiGodunov
double phiGodunov(const double r);

#define CentralDiff phiCentralDiff
double phiCentralDiff(const double r);

#define vanLeer phivanLeer
double phivanLeer(const double r);


class FluxGodunovType :
	public Flux
{
private:

	//Малое число - отсечка для наклонов:
	double thres;

	//Ссылка на функцию-ограничитель (phiGodunov, phivanLeer...)
	double(*const phi)(const double);

	//Функция численного потока, 
	//вычисляемого по решению слева (UL) и справа (UR) от разрыва,
	//значению r для переменной val
	double FluxNum(const vector<double>& fluxL, const vector<double>& fluxR, const double r, const var val);

	//Указатели на векторы решения в ячейке через одну слева
	const vector<vector<double>> *ptr_leftleftsol;

protected:
	const vector<vector<double>>& leftleftsol()  { return *ptr_leftleftsol; };

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Установка приватных указателей (чтобы работали вышеперечисленные ссылки): 
	// - на конв. потоки в своей ячейке
	// - на предельные (слева и справа) потоки в фиктивных ячейках
	// - на решение в своей ячейке и в соседних ячейках
	virtual void setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell);


public:
	//Конструктор (prm - объект, содержащий параметры задачи, 
	//            phi - функция-ограничитель)
	FluxGodunovType(const BaseParams& prm, Problem& prb,\
		double(*const phi_)(const double), double thres);

	//Деструктор
	~FluxGodunovType();

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, шаг расчета
	//заполнение приращений решений (DU) и наклонов (DV) на всех ячейках сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	//найденная скорость изменения решения и наклонов умножается на cft
	//(для реализации явного метода Эйлера следует положить cft = ptrprm->tau)
	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);
};



#endif
