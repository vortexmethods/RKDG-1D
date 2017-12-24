#ifndef FLUXLAXFRIEDRICHS_H_
#define FLUXLAXFRIEDRICHS_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"

#define LaxFriedrichs FluxLaxFriedrichs

using namespace std;

class FluxLaxFriedrichs :
    public Flux
{
private:

    //Собственные числа на всей сетке (между i-й и (i-1)-й ячейками) 
    vector<vector<double>> L;

    //Потоки Лакса-Фридрихса на левом и правом конце ячейки
    vector<double> lfL, lfR;

    //Временные переменные: векторы - скачок решения на границах i-й ячейки
    vector<double> LdU, RdU;

	//Способ вычисления скорости на границе ячеек
	SoundVelType SoundVel;

public:

    //Конструктор  (prm - объект, содержащий параметры задачи)
	FluxLaxFriedrichs(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //Деструктор
    ~FluxLaxFriedrichs();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);
};

#endif