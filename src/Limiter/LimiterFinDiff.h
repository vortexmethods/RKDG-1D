#ifndef LIMITERFINDIFF_H_
#define LIMITERFINDIFF_H_


#include "Limiter.h"

#define FinDiff LimiterFinDiff

class LimiterFinDiff :
    public Limiter
{
private:
	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, расчет значений новых наклонов всех компонент решения
	//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) ячейке сетки
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	//ВЕЗДЕ УСТАНАВЛИВАЮТСЯ НУЛЕВЫЕ НАКЛОНЫ
	void CalculateBound(const vector<vector<vector<double>>>& SOL, \
		const int cell);
public:
	//Конструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию)
	//Должен вызывать индикатор IndEverywhere, срабатывающий всюду!
    LimiterFinDiff(const BaseParams& prm, const Problem& prb, const Indicator& ind);

	//Деструктор
    ~LimiterFinDiff();
};

#endif
