#ifndef LIMITER_H_
#define LIMITER_H_

#include <vector>

#include "Indicator.h"
#include "Problem.h"

using namespace std;

class Limiter
{   
    private:
    protected:

		//Вектор значений индикаторной функции во всех ячейках сетки
        vector<double> Ind;

		//Скорректированные наклоны всех компонент решения во всех ячейках
		vector<vector<vector<double>>> SOLcorr;

		//ВИРТУАЛЬНАЯ ФУНКЦИЯ 
		//Собственно, расчет значений новых наклонов всех компонент решения
		//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) ячейке сетки
		//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
        virtual void CalculateBound(const vector<vector<vector<double>>>& SOL, \
            const int cell) = 0;
    
		//Указатель на объект, содержащий параметры задачи
		const BaseParams* ptrprm;

		//Указатель на объект, содержащий параметры задачи
		const Problem* ptrprb;

		//Указатель на индикатор (вызывается для расчета значений индикаторной функции)
		const Indicator* ptrInd;

    public:

		//Конструктор (prm - объект, содержащий параметры задачи, 
		//dimension - число консервативных переменных,
		//ind - объект, задающий индикаторную функцию)
        Limiter(const BaseParams& prm, const Problem& prb, const Indicator& ind);
        
		//Деструктор
		~Limiter();

		//Функция, вызывающая (при необходимости) расчет лимитированных значений
		//(заполнение массива VV)
		//рассчитываются по решению (UU) и наклонам (VV) на всей сетке
		void Bound(vector<vector<vector<double>>>& SOL);
};

#endif
