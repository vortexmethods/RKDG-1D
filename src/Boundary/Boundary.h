#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <vector>

#include "defs.h"
#include "Problem.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"
#include "ProblemTransfer1D.h"

using namespace std;

#define SET_NO_CONST_REF(u, constu) vector<vector<double>>& u = const_cast<vector<vector<double>>&>(constu);

class Boundary
{
protected:
	//Указатель на объект, содержащий параметры задачи
	const BaseParams* ptrprm;

	//Указатель на объект, содержащий решаемую задачу
	const Problem* ptrprb;

	//Копирование решения (U, V, W) -> (copyU, copyV, copyW)
	void CopySol(const vector<vector<double>>& SOL, vector<vector<double>>& copySOL) const;

public:
	//ВИРТУАЛЬНАЯ ФУНКЦИЯ 
	//Собственно, происходит формирование ГУ путем создания решения в фиктивной ячейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных ячейках
	//рассчитывается по решению (UU) и наклонам (VV) на всех ячейках сетки
	virtual void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const = 0;

	//Конструктор (prm - объект, содержащий параметры задачи
	//prb - решаемая задача
	Boundary(const BaseParams& prm, const Problem& prb);

	//Деструктор
	~Boundary();
};

#endif
