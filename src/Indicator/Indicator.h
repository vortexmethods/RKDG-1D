#ifndef INDICATOR_H_
#define INDICATOR_H_

#include <vector>

#include "defs.h"
#include "Problem.h"

using namespace std;

class Indicator
{
private:

protected:
	//Указатель на объект с параметрами задачи
    const BaseParams* ptrprm;
	const Problem* ptrprb;
public:
	//Конструктор (prm - объект, содержащий параметры задачи, prb - ссылка на решаемую задачу)
	Indicator(const BaseParams& prm, const Problem& prb);
    
	//Деструктор
	~Indicator();
    	
	//ВИРТУАЛЬНАЯ ФУНКЦИЯ 
	//Собственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех ячейках сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
    virtual void calc_indicator(const vector<vector<vector<double>>>& SOL, \
        vector<double>& Ind) const = 0;
};

#endif
