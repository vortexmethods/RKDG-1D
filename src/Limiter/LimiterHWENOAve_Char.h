#ifndef LIMITERHWENOAve_Char_H_
#define LIMITERHWENOAve_Char_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define HWENOAve_Char LimiterHWENOAve_Char

using namespace std;

class LimiterHWENOAve_Char :
	public Limiter
{
private:
	double wg; //Степень, в которую возводятся слагаемые при расчете лимитированных моментов

	const double weps = 1e-7; //Малое положительное число

							  //РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
							  //Собственно, расчет значений новых моментов всех компонент решения
	void CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell);

public:
	//Конструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию,
	//degree - параметр лимитера, присваемый переменной wg)
	LimiterHWENOAve_Char(const BaseParams& prm, const Problem& prb, Indicator& ind, double degree);

	//Деструктор
	~LimiterHWENOAve_Char();
};

#endif
