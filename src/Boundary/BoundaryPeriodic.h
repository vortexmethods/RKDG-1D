#ifndef BOUNDARYPERIODIC_H_
#define BOUNDARYPERIODIC_H_

#include <vector>

#include "Boundary.h"

#define Periodic BoundaryPeriodic 

class BoundaryPeriodic :
	public Boundary
{
protected:

public:
	BoundaryPeriodic(const BaseParams& prm, const Problem& prb);
	~BoundaryPeriodic();

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, происходит формирование ГУ путем создания решения в фиктивной ячейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных ячейках
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif