#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

class CMatrix {
private:
	double** matrix;
	short row;
	short col;

	double* determinant;
	double** same;
	double** unitMatrix;

public:
	CMatrix(short i, short j);
	CMatrix(const CMatrix& copy);

	CMatrix& operator = (const CMatrix& what);
	CMatrix operator + (const CMatrix& what) const;
	CMatrix operator - (const CMatrix& what) const;
	CMatrix operator * (double k) const;
	friend CMatrix operator * (double k, const CMatrix& what);
	CMatrix operator * (const CMatrix& what) const;

	CMatrix trn() const;
	double det();
	CMatrix inverse();

	friend std::istream& operator >> (std::istream& is, const CMatrix& what);
	friend std::ostream& operator << (std::ostream& os, const CMatrix& what);

	~CMatrix();

};

#endif // MATRIX_H
