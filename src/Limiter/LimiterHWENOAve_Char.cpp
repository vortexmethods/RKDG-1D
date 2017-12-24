#include "LimiterHWENOAve_Char.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера HWENOAve_Char
LimiterHWENOAve_Char::LimiterHWENOAve_Char(const BaseParams& prm, const Problem& prb, \
	Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENOAve_Char::~LimiterHWENOAve_Char()
{}


void LimiterHWENOAve_Char::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
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

	 // Подготовка к расчёту по методу HWENO_SC --- проекция решения на базис собственных векторов.
	for (int j = 0; j < nshape; ++j)
	{
		prodMatrVec(LL, leftSol[j],   auxSOL_L[j]);
		prodMatrVec(LL, troublSol[j], auxSOL_T[j]);
		prodMatrVec(LL, rightSol[j],  auxSOL_R[j]);
	}

	// Линейные веса. //позднее:для правой точки
	double gamma0 = 9.0 / 80.0;
	double gamma1 = 29.0 / 80.0;
	double gamma2 = 42.0 / 80.0;

	if (nshape == 3)
		// Цикл по консервативным переменным 
		for (int val = 0; val < dim; ++val)
		{
			cout << "Приносим свои извинения, данный метод не обобщён на квадратичный случай." << endl;
		}
	else
	{
		for (int val = 0; val < dim; ++val)
		{
			// Denoting slopes and averages of the solution at the 3-point stencil.

			const double& vl = auxSOL_L[1][val] * 2.0 / h;
			const double& v = auxSOL_T[1][val] * 2.0 / h;
			const double& vr = auxSOL_R[1][val] * 2.0 / h;

			const double& yl = auxSOL_L[0][val];
			const double& y = auxSOL_T[0][val];
			const double& yr = auxSOL_R[0][val];

			// Evaluating hermite polynomials of the degree 2.
			// For the right.
			double p0p = (13.0*y - 7.0*yl - 4.0*h*vl) / 6.0;
			double p1p = (5.0*y + 2.0*yr - yl) / 6.0;
			double p2p = (y + 5.0*yr - 2.0*h*vr) / 6.0;
			// And for the left.
			double p0m = (2.0*h*vl + y + 5.0*yl) / 6.0;
			double p1m = (5.0*y + 2.0*yl - yr) / 6.0;
			double p2m = (4.0*h*vr + 13.0*y - 7.0*yr) / 6.0;

			// Calculating the smothness indicators for modified polynomials at the troubled cell. 
			double beta0 = (16.0*vl*vl*h*h + 25.0*(y - yl)*(y - yl) + 38.0*vl*h*(yl - y)) / 3.0;
			double beta1 = (4.0*yl*yl - 13.0*yl*y + 13.0*y*y + 5.0*yl*yr - 13.0*y*yr + 4.0*yr*yr) / 3.0;
			double beta2 = (16.0*vr*vr*h*h + 25.0*(y - yr)*(y - yr) + 38.0*vr*h*(y - yr)) / 3.0;

			// Calculating the nonlinear weights. Normalization of weights will be carried out (applied?) during the final weightened interpolation.
			//For the left.
			double w0m = gamma2 / pow(weps + beta0, wg);
			double w1m = gamma1 / pow(weps + beta1, wg);
			double w2m = gamma0 / pow(weps + beta2, wg);
			// And for the right.
			double w0p = gamma0 / pow(weps + beta0, wg);
			double w1p = gamma1 / pow(weps + beta1, wg);
			double w2p = gamma2 / pow(weps + beta2, wg);

			// And the final weightened interpolation of the troubled cell boundary vaules. 
			double Uplus = (p0p*w0p + p1p*w1p + p2p*w2p) / (w0p + w1p + w2p); // For the right (i+1/2).
			double Uminus = (p0m*w0m + p1m*w1m + p2m*w2m) / (w0m + w1m + w2m); // And for the left (i-1/2).

			// The resulting slope of solution is evaluated as the linear combination of the existed slopes.
			auxSOLcorr[1][val] = 0.5 * (Uplus - Uminus);
		}//for val

		prodMatrVec(RR, auxSOLcorr[1], SOLcorr[cell][1]);
	}
}
