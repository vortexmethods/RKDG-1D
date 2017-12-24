#ifndef TIMESTEP_H_
#define TIMESTEP_H_

#include <vector>

#include "Flux.h"
#include "Limiter.h"
#include "Boundary.h"

class Timestep
{
protected:
	//Указатель на объект с параметрами задачи
	const BaseParams* ptrprm;
	
	//Количество консервативных переменных
	int dim;

	//Указатель на ГУ
	const Boundary* ptrbnd;

	//Векторы приращений средних значенй и наклонов:
	vector<vector<vector<double>>> DSOL;
	
	//Снос решения U, V на фиктивную ячейку	
	void ApplyBoundary(const vector<vector<vector<double>>>& SOL);

public:
	//Конструктор
	Timestep(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	
	//Деструктор
	virtual ~Timestep();

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ 
	//Собственно, шаг расчета соответствующим численным методом
	//Вычисление нового решения (Unew) и наклонов (Vnew) на всех ячейках сетки
	//рассчитывается по решению (U) и наклонам (V) на всех ячейках сетки.
	//Шаг по времени tau
	//Для расчета потоков используется численный поток method, лимитер lim.
	//(в реализации для расчета приращений следует использовать method.step(...)
	virtual void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim) = 0;	
};

#endif