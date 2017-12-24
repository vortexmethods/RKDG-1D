#ifndef INDIKATORKRIVODONOVA_H_
#define INDICATORKRIVODONOVA_H_

#include "Indicator.h"
#include "defs.h"


#define Krivodonova IndicatorKrivodonova

class IndicatorKrivodonova :
    public Indicator
{
private:
	//"Чувствительная" переменная, по которой рассчитывается значение индикаторной функции
	var sens;

public:
	//Конструктор (prm - объект, содержащий параметры задачи, val - "чувствительная" переменная
	//             prb - решаемая задача (т.к. важно направление переноса!)
	IndicatorKrivodonova(const BaseParams& prm, const Problem& prb, const var val);

	//Деструктор
    ~IndicatorKrivodonova();
    
	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех ячейках сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	void calc_indicator(const vector<vector<vector<double>>>& SOL, \
		vector<double>& Ind) const;
};

#endif