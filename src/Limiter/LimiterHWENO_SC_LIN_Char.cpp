#include "LimiterHWENO_SC_LIN_Char.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера HWENO_SC_LIN_Char
LimiterHWENO_SC_LIN_Char::LimiterHWENO_SC_LIN_Char(const BaseParams& prm, const Problem& prb, \
	const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENO_SC_LIN_Char::~LimiterHWENO_SC_LIN_Char()
{}


void LimiterHWENO_SC_LIN_Char::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
{
	// Некоторые базовые параметры.
	int nx = ptrprm->nx;
	double h = ptrprm->h;
	int dim = ptrprb->dim;
	int nshape = ptrprb->nshape;

	//Определяем соседей проблемной ячейки. При необходимости, посредством фиктивных ячеек.
	const vector<vector<double>>& leftSol = (cell > 0) ? SOL[cell - 1] : SOL[nx];
	const vector<vector<double>>& troublSol = SOL[cell];
	const vector<vector<double>>& rightSol = (cell < nx - 1) ? SOL[cell + 1] : SOL[nx + 1];

	// Линейные веса.
	double ak = 0.01;
	double gamma0 = ak;
	double gamma1 = 1.0 - 2.0*ak;
	double gamma2 = ak;

	//Матрицы для собственных векторов.
	vector<vector<double>> LL;
	vector<vector<double>> RR;
	LL.resize(dim);
	RR.resize(dim);
	for (int j = 0; j < dim; ++j)
	{
		LL[j].resize(dim);
		RR[j].resize(dim);
	}// for j

	 // Заполняем матрицы собственных векторов.
	ptrprb->EigenMatricies(troublSol, LL, RR);

	// Массив не- и лимитированных полиномов в характеристическом пространстве.
	vector<vector<double>> auxSOL_L, auxSOL_T, auxSOL_R;
	vector<vector<double>> auxSOLcorr;
	auxSOLcorr.resize(nshape);
	auxSOL_L.resize(nshape); auxSOL_T.resize(nshape); auxSOL_R.resize(nshape);
	for (int j = 0; j < nshape; ++j)
	{
		auxSOLcorr[j].resize(dim);
		auxSOL_L[j].resize(dim);
		auxSOL_T[j].resize(dim);
		auxSOL_R[j].resize(dim);
	}// for j

	 // Подготовка к расчёту по методу HWENO_SC_LIN_Char --- проекция решения на базис собственных векторов.
	for (int j = 0; j < nshape; ++j)
	{
		prodMatrVec(LL, leftSol[j],   auxSOL_L[j]);
		prodMatrVec(LL, troublSol[j], auxSOL_T[j]);
		prodMatrVec(LL, rightSol[j],  auxSOL_R[j]);
	}


	// Цикл по консервативным переменным 
	if (nshape > 2)
	{
		cout << "Приносим свои извинения, данный метод не обобщён на квадратичный случай." << endl;
	}
	else
	{
		for (int val = 0; val < dim; ++val)
		{
			// Переобозначение моментов решения с ячеек шаблона.
			const double& ul = auxSOL_L[0][val];
			const double& u  = auxSOL_T[0][val];
			const double& ur = auxSOL_R[0][val];

			const double& vl = auxSOL_L[1][val];
			const double& v  = auxSOL_T[1][val];
			const double& vr = auxSOL_R[1][val];


			// Модифицированные полиномы соседних ячеек (по МНК)
			double pLv = (vl + 6.0 * (u - ul)) / 13.0;
			double p = v;
			double pRv = (vr + 6.0 * (ur - u)) / 13.0;

			// Вычисляем индикаторы гладкости полиномов решения в проблемной ячейке.
			double beta0 = 4.0 * pLv*pLv;
			double beta1 = p*p * 4.0;
			double beta2 = 4.0 * pRv*pRv;

			// Вычисляем нелинейные веса. Нормировка будет осуществлена ниже.
			double w0 = gamma0 / pow(weps + beta0, wg);
			double w1 = gamma1 / pow(weps + beta1, wg);
			double w2 = gamma2 / pow(weps + beta2, wg);
			// Нормировка весов.
			double WWW = (w0 + w1 + w2);
			w0 = w0 / WWW; w1 = w1 / WWW; w2 = w2 / WWW;

			//// Честно сделанный пересчёт моментов с учётом привязки базисных функций к ячейкам.
			double t1, t2, t3, t4;

			t1 = w0*pLv;
			t2 = w1*p;
			t3 = w2*pRv;
			t4 = t1 + t2 + t3;
			auxSOLcorr[1][val] = t4;
		}// for val

		prodMatrVec(RR, auxSOLcorr[1], SOLcorr[cell][1]);
	}
}
