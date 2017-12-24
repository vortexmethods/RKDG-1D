#ifndef LIMITERHWENO_H_
#define LIMITERHWENO_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define HWENO LimiterHWENO

using namespace std;

class LimiterHWENO :
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
	LimiterHWENO(const BaseParams& prm, const Problem& prb, Indicator& ind, double degree);

	//Деструктор
	~LimiterHWENO();

	double centrdiff(double yL, double yR, double h)
	{
		return 0.5 * (yR - yL) / h;
	}

	double diff2(double yL, double y, double yR, double h)
	{
		return (yR - 2.0*y + yL) / (h*h);
	}

	double backbackdiff(double yLL, double yL, double y, double h)
	{
		return 0.5 * (3.0*y - 4.0*yL + yLL) / h;
	}

	double forwforwdiff(double y, double yR, double yRR, double h)
	{
		return 0.5 * (-3.0*y + 4.0*yR - yRR) / h;
	}
};

#endif