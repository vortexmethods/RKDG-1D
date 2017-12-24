#include "LimiterHWENO.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера HWENO
LimiterHWENO::LimiterHWENO(const BaseParams& prm, const Problem& prb, \
	Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENO::~LimiterHWENO()
{}


void LimiterHWENO::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
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
	double gamma0m = 9.0 / 16.0;//0.09;//0.3;//
	double gamma1m = 6.0 / 16.0;//0.9; //0.6;//
	double gamma2m = 1.0 / 16.0;//0.01;//0.1;//

	double gamma0p = 1.0 / 16.0;//0.01;//0.1;//
	double gamma1p = 6.0 / 16.0;//0.9; //0.6;//
	double gamma2p = 9.0 / 16.0;//0.09;//0.3;//

	if (nshape == 3)
		// Цикл по консервативным переменным 
		for (int val = 0; val < dim; ++val)
		{
			cout << "Приносим свои извинения, данный метод не обобщён на квадратичный случай." << endl;
		}
	else
		for (int val = 0; val < dim; ++val)
		{
			// Denoting slopes and averages of the solution at the 3-point stencil.

			const double& dyl =  2.0*leftSol[1][val] / h;
			const double& dy = 2.0*troublSol[1][val] / h;
			const double& dyr = 2.0*rightSol[1][val] / h;
			
			const double& yl = leftSol[0][val];
			const double& y = troublSol[0][val];
			const double& yr = rightSol[0][val];

			// Evaluating hermite polynomials of the degree 2.
			// For the right (p) and left (m)
			double p0m = (h*dyl + 3.0*yl + y) / 4.0;
			double p0p = (-3.0*dyl*h - 5.0*yl + 9.0*y) / 4.0;

			double p1m = 0.125*(3.0*yl + 6.0*y - yr);
			double p1p = 0.125*(-yl + 6.0*y + 3.0*yr);

			double p2m = (3.0*dyr*h - 5.0*yr + 9.0*y) / 4.0;
			double p2p = (-dyr*h + 3.0*yr + y) / 4.0;
			
			// Calculating the smothness indicators for modified polynomials at the troubled cell. 
			double beta0 = (16.0*dyl*dyl*h*h + 25.0*(y - yl)*(y - yl) + 38.0*dyl*h*(yl - y)) / 3.0;
			double beta1 = (4.0*yl*yl - 13.0*yl*y + 13.0*y*y + 5.0*yl*yr - 13.0*y*yr + 4.0*yr*yr) / 3.0;
			double beta2 = (16.0*dyr*dyr*h*h + 25.0*(y - yr)*(y - yr) + 38.0*dyr*h*(y - yr)) / 3.0;

			// Calculating the nonlinear weights. Normalization of weights will be carried out (applied?) during the final weightened interpolation.
			//For the left.
			double w0m = gamma0m / pow(weps + beta0, wg);
			double w1m = gamma1m / pow(weps + beta1, wg);
			double w2m = gamma2m / pow(weps + beta2, wg);
			//And for the right.
			double w0p = gamma0p / pow(weps + beta0, wg);
			double w1p = gamma1p / pow(weps + beta1, wg);
			double w2p = gamma2p / pow(weps + beta2, wg);

			// And the final weightened interpolation of the troubled cell boundary vaules. 
			double Uplus = (p0p*w0p + p1p*w1p + p2p*w2p) / (w0p + w1p + w2p); // For the right (i+1/2).
			double Uminus = (p0m*w0m + p1m*w1m + p2m*w2m) / (w0m + w1m + w2m); // And for the left (i-1/2).

			// The resulting slope of solution is evaluated as the linear combination of the existed slopes.
			SOLcorr[cell][1][val] = 0.5 * (Uplus - Uminus);
		}
}
