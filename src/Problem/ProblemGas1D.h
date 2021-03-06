#ifndef PROBLEMGAS1D_H_
#define PROBLEMGAS1D_H_

#define Gas1D ProblemGas1D

#include "Problem.h"
class ProblemGas1D :
	public Problem
{
private:
	double gamma;

	//скорость звука на границе (left)-й и (right)-й ячеек
	double c_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	double c_semisum(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//скорость звука в i-й ячейке
	double c(const vector<vector<double>>& sol) const;

	//Энтальпия звука на границе (left)-й и (right)-й ячеек
	double h_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//Энтальпия звука на i-й ячейке
	double h(const vector<vector<double>>& sol) const;

	//скорость на границе ячеек
	vector<double> v_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//скорость на i-й ячейке
	vector<double> v(const vector<vector<double>>& sol) const;

	//Квадрат скорости на границе (left)-й и (right)-й ячеек
	double v2_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//Квадрат скорости на i-й ячейке
	double v2(const vector<vector<double>>& sol) const;

	//Уточненная скорость звука на разрыве
	double d_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
   
protected:

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Вычисление конвективного потока по заданному вектору решения
	virtual void getFlux(const vector<double>& U, vector<double>& Flux) const;
	virtual vector<double> getFlux(const vector<double>& U) const;

public:
	ProblemGas1D(const BaseParams& prm, int order, double gam, vector<std::function<double(const double)>> initv);
	~ProblemGas1D();

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Предельное значение (sd = left/right) q-й компоненты решения,
	//рассчитываемое по решению (U) и наклону (V)
	virtual inline double side_val(const vector<vector<double>>& sol, var q, side sd) const;
	virtual inline double val(const vector<vector<double>>& sol, var q) const;
	virtual inline vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const;
	
	//Значения левых (LW) и правых (RW) собственных векторов на всей сетке 
	//(на левых границах ячеек - между (i-1)-й и i-й ячейками
	//имеющие номера, указанные в списке инициализации (по возрастанию)
	//рассчитываются по решениям (UU) и наклонам (VV) на всей сетке
	//(для метода КИР)
	void omega(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& LW, \
		vector<vector<vector<double>>>& RW, \
		const initializer_list<int>& list) const;
	// Матрицы собственных векторов на ячейке
	void EigenMatricies(const vector<vector<double>>& sol, \
		vector<vector<double>>& LL, \
		vector<vector<double>>& RR) const;
	
	//Значения собственных чисел (LL) на всей сетке  
	//имеющие номера, указанные в списке инициализации (по возрастанию)
	//(на левых границах ячеек - между (i-1)-й и i-й ячейками
	//рассчитываются по решениям (UU) и наклонам (VV) на всей сетке	
	void lambda(const vector<vector<vector<double>>>& SOL, \
		const SoundVelType soundveltype, \
		vector<vector<double>>& LL, \
		const initializer_list<int>& list) const;

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Возвращает вектор конс. переменных в точке x, вычисляемый по заданным в объекте Param начальным условиям
	virtual vector<double> initial_var(const double x) const;

	//ВИРТУАЛЬНАЯ ФУНКЦИЯ
	//Печать в телефайл
	virtual void more_to_file(ostream& str, \
		const vector<vector<vector<double>>>& SOL, int cell) const;
};

#endif
