#ifndef BOUNDARYWALL_H_
#define BOUNDARYWALL_H_

#include <vector>

#include "Boundary.h"

#define Wall BoundaryWall 

class BoundaryWall :
	public Boundary
{
protected:

public:
	BoundaryWall(const BaseParams& prm, const Problem& prb);
	~BoundaryWall();

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, происходит формирование ГУ путем создания решения в фиктивной ячейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных ячейках
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif