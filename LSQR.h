#pragma once
#include "stdafx.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>


class LSQR
{
public:
	Eigen::VectorXd x;
	LSQR(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, const Eigen::VectorXd &x_ini, const double &gamma);
	Eigen::VectorXd SolutionX();
};

