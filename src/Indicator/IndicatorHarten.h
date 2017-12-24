#ifndef INDICATORHARTEN_H_
#define INDICATORHARTEN_H_

#include "Indicator.h"
#include "defs.h"

#define Harten IndicatorHarten

class IndicatorHarten :
	public Indicator
{

private:
	double mult; //множитель для сравнения наклонов в соседних ячейках

	//"Чувствительная" переменная, по которой рассчитывается значение индикаторной функции
	var sens;

public:
	//Конструктор 
	//  prm - объект, содержащий параметры задачи, 
	//  val - "чувствительная" переменная
	//  alpha - коэффициент в индикаторе Хартена
	IndicatorHarten(const BaseParams& prm, const Problem& prb, const var val, const double alpha);

	//Деструктор
	~IndicatorHarten();

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех ячейках сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
//	void calc_indicator(const vector<vector<double>>& UU, \
//		const vector<vector<double>>& VV, \
//		vector<double>& Ind) const;

    void calc_indicator(const vector<vector<vector<double>>>& SOL, \
        vector<double>& Ind) const;
};

#endif

