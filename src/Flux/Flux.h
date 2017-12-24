#ifndef FLUX_H_
#define FLUX_H_

#include <vector>

#include "defs.h"
#include "Problem.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"

using namespace std;

class Flux
{
private:
	//Указатели на векторы конв. потоков на трех ячейках
	vector<double> *ptr_leftfluxR, *ptr_myflux, *ptr_myfluxL, *ptr_myfluxR, *ptr_rightfluxL;

	//Указатели на векторы решения в трех ячейках
	const vector<vector<double>> *ptr_mysol, *ptr_leftsol, *ptr_rightsol;


protected:
	//Указатель на объект с общими параметрами
	const BaseParams* ptrprm;
			
	//Указатель на объект задача
	Problem* ptrprb;
	ProblemGas1D* ptrprb_toGas;
	ProblemMHD1D* ptrprb_toMHD;

	//Обертка для вызова одноименной функции из класса Problem
    double side_val(const vector<vector<double>>& sol, var q, side sd) const
	{
		return ptrprb->side_val(sol, q, sd);
	}

    vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const
	{
		return ptrprb->gauss_val(sol, Lcoord);
	}

	//Ссылки на векторы конв. потоков на трех ячейках
	vector<double>& myfluxL() { return *ptr_myfluxL; };
	vector<double>& myfluxR() { return *ptr_myfluxR; };
	vector<double>& myflux()  { return *ptr_myflux; };
	vector<double>& leftfluxR()  { return *ptr_leftfluxR; };
	vector<double>& rightfluxL() { return *ptr_rightfluxL; };

	//Ссылки на векторы решения на трех ячейках
	const vector<vector<double>>& mysol() { return *ptr_mysol; };
	const vector<vector<double>>& leftsol()  { return *ptr_leftsol; };
	const vector<vector<double>>& rightsol() { return *ptr_rightsol; };
	

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ С РЕАЛИЗАЦИЕЙ
	//Установка приватных указателей (чтобы работали вышеперечисленные ссылки): 
	// - на конв. потоки в своей ячейке
	// - на предельные (слева и справа) потоки в фиктивных ячейках
	// - на решение в своей ячейке и в соседних ячейках
	virtual void setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell);




public: 

	//Конструктор (prm - объект, содержащий параметры задачи)
	Flux(const BaseParams& prm, Problem& prb);
	
	//Деструктор
	~Flux();

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ 
	//Заполнение приращений решений (DU) и наклонов (DV) на всех ячейках сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	//найденная скорость изменения решения и наклонов умножается на cft
	//(для реализации явного метода Эйлера следует положить cft = ptrprm->tau)
	virtual void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft) = 0;
};

#endif

