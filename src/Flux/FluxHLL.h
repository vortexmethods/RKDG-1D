#ifndef FLUXHLL_H_
#define FLUXHLL_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

using namespace std;

#define HLL FluxHLL

class FluxHLL :
    public Flux
{
private:

    //Собственные числа на всей сетке (между i-й и (i-1)-й ячейками) 
    vector<vector<double>> L;

    //Потоки HLL на левом и правом конце ячейки
    vector<double> hllL, hllR;

    //Временные переменные: векторы - скачок решения на границах i-й ячейки
    vector<double> LdU, RdU;

	//Способ вычисления скорости на границе ячеек
	SoundVelType SoundVel;

public:

    //Конструктор (prm - объект, содержащий параметры задачи)
	FluxHLL(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //Деструктор
    ~FluxHLL();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);

};


#endif