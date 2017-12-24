#ifndef BOUNDARYSOFT_H_
#define BOUNDARYSOFT_H_

#include <vector>

#include "Boundary.h"

#define Soft BoundarySoft 

class BoundarySoft :
	public Boundary
{
protected:

public:
	BoundarySoft(const BaseParams& prm, const Problem& prb);
	~BoundarySoft();

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, происходит формирование ГУ путем создания решения в фиктивной ячейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных ячейках
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif