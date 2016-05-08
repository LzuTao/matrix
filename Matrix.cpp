#include <iomanip>
#include <cstring>
#include "Matrix.h"
#include "Vector.h"

const short WIDTH = 10;
const short PRECISION = 3;

CMatrix operator * (const CMatrix& m_Mat, const CVector& m_Vec);

CMatrix::CMatrix(): matrix(NULL), row(0), col(0),
	determinant(NULL), same(NULL), unitMatrix(NULL){
}

CMatrix::CMatrix(short i = 0, short j = 0): matrix(NULL), row(i), col(j),
	determinant(NULL), same(NULL), unitMatrix(NULL) {

	if ((row <= 0) || (col <= 0)) {
		row = col = 0;
	} else {
		matrix = new double*[row];
		short i;
		for (i = 0; i < row; ++i) {
			matrix[i] = new double[col];

			memset(matrix[i], 0, col * sizeof(double));
		}
	}
}

CMatrix::CMatrix(const CMatrix& copy): matrix(NULL), row(copy.row),
	col(copy.col), determinant(NULL), same(NULL), unitMatrix(NULL) {

	if((row == col) && (copy.determinant != NULL)) {
		determinant = new double;
		*determinant = *(copy.determinant);

		if (copy.same != NULL && copy.unitMatrix != NULL) {
			same = new double*[col];
			unitMatrix = new double*[col];

			for (short i = 0; i < col; ++i) {
				same[i] = new double[col];
				unitMatrix[i] = new double[col];

				memcpy(same[i], copy.same[i], col * sizeof(double));
				memcpy(unitMatrix[i], copy.unitMatrix[i], col * sizeof(double));
			}
		}
	}

	if (copy.matrix != NULL) {
		matrix = new double*[row];
		short i;

		for (i = 0; i < row; ++i) {
			matrix[i] = new double[col];

			memcpy(matrix[i], copy.matrix[i], col * sizeof(double));

		}
	}
}

CMatrix& CMatrix::operator = (const CMatrix& m_Mat) {
	if (this != &m_Mat) {
		// Delete and create a new Matrix
		short i, j;
		if (!((row == m_Mat.row) && (col == m_Mat.col))) {
			for (i = 0; i < row; ++i) {
				delete[] matrix[i];
			}

			delete[] matrix;

			row = m_Mat.row;
			col = m_Mat.col;

			matrix = new double*[row];
			for (i = 0; i < row; ++i) {
				matrix[i] = new double[col];
			}

		}

		for (i = 0; i < row; ++i) {
			for (j = 0; j < col; ++j) {
				matrix[i][j] = m_Mat.matrix[i][j];
			}
		}

	}

	return *this;
}

CMatrix CMatrix::operator + (const CMatrix& what) const {
	CMatrix m_Mat(row, col);
	if ((row == what.row) && (col == what.col)) {
		short i, j;
		for (i = 0; i < row; ++i) {
			for (j = 0; j < col; ++j) {
				m_Mat.matrix[i][j] = matrix[i][j] + what.matrix[i][j];
			}

		}
	}

	return m_Mat;
}

CMatrix CMatrix::operator - (const CMatrix& what) const {
	CMatrix m_Mat(row, col);
	if ((row == what.row) && (col == what.col)) {
		short i, j;
		for (i = 0; i < row; ++i) {
			for (j = 0; j < col; ++j) {
				m_Mat.matrix[i][j] = matrix[i][j] - what.matrix[i][j];
			}

		}
	}

	return m_Mat;
}

CMatrix CMatrix::operator * (double k) const {
	short i, j;
	CMatrix m_Mat(row, col);
	for (i = 0; i < row; ++i) {
		for (j = 0; j < col; ++j) {
			m_Mat.matrix[i][j] = matrix[i][j] * k;
		}

	}

	return m_Mat;
}

CMatrix operator * (double k, const CMatrix& what) {
	return what * k;
}

CMatrix CMatrix::operator * (const CMatrix& what) const {
	CMatrix m_Mat(row, what.col);
	if (col == what.row) {
		short i, j, k;
		for (i = 0; i < row; ++i) {
			for (k = 0; k < what.col; ++k) {
				for (j = 0; j < col; ++j) {
					m_Mat.matrix[i][k] += matrix[i][j] * what.matrix[j][k];
				}
			}
		}
	}

	return m_Mat;
}

CMatrix CMatrix::trn() const {
	CMatrix m_Mat(col, row);
	short i, j;
	for (i = 0; i < col; ++i) {
		for (j = 0; j < row; ++j) {
			m_Mat.matrix[i][j] = matrix[j][i];
		}
	}

	return m_Mat;
}

double CMatrix::det() {
	if (determinant != NULL) {
		return *determinant;
	}

	if ((matrix == NULL) || (row != col)) {
		std::cout << "Error" << std::flush;
		return 0;
	}

	if (col == 2) {
		determinant = new double;
		*determinant = matrix[0][0] * matrix[1][1]
					 - matrix[1][0] * matrix[0][1];

		return *determinant;
	}

	if (col == 3) {
		determinant = new double;
		*determinant =
		 matrix[0][0] * matrix[1][1] * matrix[2][2] +
		 matrix[0][1] * matrix[1][2] * matrix[2][0] +
		 matrix[0][2] * matrix[1][0] * matrix[2][1] -
		 matrix[2][0] * matrix[1][1] * matrix[0][2] -
		 matrix[2][1] * matrix[1][2] * matrix[0][0] -
		 matrix[2][2] * matrix[1][0] * matrix[0][1];

		return *determinant;
	}

	short i, j, k;


	// Create the same matrix of class
	if (same == NULL) {
		same = new double*[col];
		for (i = 0; i < col; ++i) {
			same[i] = new double[col];
			memcpy(same[i], matrix[i], col * sizeof(double));
		}
	}
	// Create unitMatrix with the same level
	unitMatrix = new double*[col];
	for (i = 0; i < col; ++i) {
		unitMatrix[i] = new double[col];
		memset(unitMatrix[i], 0, col * sizeof(double));
		unitMatrix[i][i] = 1;
	}

	// Create triangular matrix
	double m_determinant = 1, temp;
	for (k = 1; k < col; ++k) {
		if (!same[k - 1][k - 1]) {
			for (i = k; i < col; ++i) {
				if (same[i][k]) {
					for (j = 0; j < col; ++j) {
						same[k - 1][j]  += same[col - 1][j];

						unitMatrix[k - 1][j]  += unitMatrix[col - 1][j];
					}

					break;
				} else {
					m_determinant = 0;
					break;
				}
			}

			if (m_determinant == 0) {
				break;
			}

			m_determinant = - m_determinant;
		} else {
			for(i = k; i < col; ++i) {
				temp = -same[i][k-1] / same[k-1][k-1];
				for(j = k - 1; j < col; ++j) {
					same[i][j] += temp * same[k-1][j];
				}

				for (j = 0; j < col; ++j) {
					unitMatrix[i][j] += temp * unitMatrix[k - 1][j];
				}
			}
		}
	}


	for (i = 0; i < col; ++i) {
		m_determinant *= same[i][i];
	}

	determinant = new double(m_determinant);

	if (m_determinant == 0) {
		// If determinant == 0, it will not have inverse matrix
		for (i = 0; i < row; ++i) {
			delete[] same[i];
			delete[] unitMatrix[i];
		}

		delete[] same;
		delete[] unitMatrix;

		return 0;
	}

	return m_determinant;
}

CMatrix CMatrix::inverse() {
	this->det();

	if (determinant != NULL) {
		CMatrix m_Mat(row, row);
		if (*determinant == 0) {
			return m_Mat;
		}


		if (col == 2) {
			m_Mat.matrix[0][0] =  matrix[1][1] / *determinant;
			m_Mat.matrix[0][1] = -matrix[0][1] / *determinant;

			m_Mat.matrix[1][0] = -matrix[1][0] / *determinant;
			m_Mat.matrix[1][1] =  matrix[0][0] / *determinant;

			return m_Mat;
		}

		if (col == 3) {
			m_Mat.matrix[0][0] =  (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]) / *determinant;
			m_Mat.matrix[0][1] =  (matrix[0][2] * matrix[2][1] - matrix[2][2] * matrix[0][1]) / *determinant;
			m_Mat.matrix[0][2] =  (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) / *determinant;

			m_Mat.matrix[1][0] =  (matrix[1][2] * matrix[2][0] - matrix[2][2] * matrix[1][0]) / *determinant;
			m_Mat.matrix[1][1] =  (matrix[0][0] * matrix[2][2] - matrix[2][0] * matrix[0][2]) / *determinant;
			m_Mat.matrix[1][2] =  (matrix[0][2] * matrix[1][0] - matrix[1][2] * matrix[0][0]) / *determinant;

			m_Mat.matrix[2][0] =  (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]) / *determinant;
			m_Mat.matrix[2][1] =  (matrix[0][1] * matrix[2][0] - matrix[2][1] * matrix[0][0]) / *determinant;
			m_Mat.matrix[2][2] =  (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) / *determinant;

			return m_Mat;
		}

		short i, j;
		double temp;

		// Creating No. 1 on the main diagonal of the matrix
		if (!isUnitMatrix()) {
			for (i = 0; i < col; ++i) {
				temp = same[i][i];
				if (temp != 1) {
					for (j = col - 1; j >= i; --j) {
						same[i][j] /= temp;
					}

					for (j = 0; j < col; ++j) {
						unitMatrix[i][j] /= temp;
					}
				}
			}

			// Starting Jordan's process
			short k;
			for (k = col - 1; k > 0; --k) {
				for (i = k - 1; i >= 0; --i) {
					temp = - same[i][k] / same[k][k];
					same[i][k] = 0;

					for (j = col - 1; j >= 0; --j) {
						unitMatrix[i][j] += temp * unitMatrix[k][j];
					}
				}
			}
		}

		// memcpy(m_Mat.matrix, unitMatrix, col * col * sizeof(double));
		for (i = 0; i < col; ++i) {
			memcpy(m_Mat.matrix[i], unitMatrix[i], col * sizeof(double));
		}

		return m_Mat;
	}

	return CMatrix(row, row);
}

bool CMatrix::isUnitMatrix() {
	if (row != col && row <= 0) {
		return false;
	}

	// Create the same matrix
	short i;
	if (same == NULL) {
		same = new double*[col];
		for (i = 0; i < col; ++i) {
			same[i] = new double[col];
			memcpy(same[i], matrix[i], col * sizeof(double));
		}
	}

	short j;
	for (i = 0; i < col; ++i) {
		for (j = 0; j < col; ++j) {
			if (i != j && same[i][j] != 0) {
				return false;
			}
		}

		if (same[i][i] != 1) {
			return false;
		}
	}

	return true;
}

std::istream& operator >> (std::istream& is, const CMatrix& what) {
	short i, j;
	for (i = 0; i < what.row; ++i) {
		for (j = 0; j < what.col; ++j) {
			is >> what.matrix[i][j];
		}

	}

	return is;
}

std::ostream& operator << (std::ostream& os, const CMatrix& what) {
	if (what.matrix != NULL) {
		short val1 = what.col * WIDTH;
		short val2 = what.col - 1;

		os << "--" << std::setfill(' ') << std::setw(val1) << "--" << std::endl;
		os << std::fixed << std::setprecision(PRECISION);
		short i, j;
		for (i = 0; i < what.row; ++i) {
			os << "|";
			for (j = 0; j < val2; ++j) {
				os << std::setw(WIDTH) << what.matrix[i][j];
			}

			os << std::setw(WIDTH) << what.matrix[i][val2] << '|' << std::endl;
		}

		os << "--" << std::setfill(' ') << std::setw(val1) << "--" << std::endl;
	}

	return os;
}

CMatrix::~CMatrix() {
	short i;
	if (matrix != NULL) {
		for (i = 0; i < row; ++i) {
			delete[] matrix[i];
		}

		delete[] matrix;
	}

	if (same != NULL) {
		for (i = 0; i < row; ++i) {
			delete[] same[i];
		}

		delete[] same;
	}

	if (unitMatrix != NULL) {
		for (i = 0; i < row; ++i) {
			delete[] unitMatrix[i];
		}

		delete[] unitMatrix;
	}
	delete determinant;
}

double CMatrix::getMat(short i, short j) const {
	if (matrix != NULL && i >= 0 && j >= 0 && i < row && j < col) {
		return matrix[i][j];
	}

	std::cout << "Error" << std::flush;
	return 0;
}

short CMatrix::getRow() const {
	return row;
}

short CMatrix::getCol() const {
	return col;
}

void CMatrix::setMat(short i, short j, double k) {
	if (matrix != NULL && i >= 0 && j >= 0 && i < row && j < col) {
		matrix[i][j] = k;
	}
}

CMatrix operator * (const CMatrix& m_Mat, const CVector& m_Vec) {
	CMatrix newMat(m_Mat.getRow(), m_Vec.getDim());
	if (m_Mat.getCol() == 1) {
		short i, k;
		for (i = 0; i < m_Mat.getRow(); ++i) {
			for (k = 0; k < m_Vec.getDim(); ++k) {
				newMat.setMat(i, k, newMat.getMat(i, k) + m_Mat.getMat(i, 0) * m_Vec.getVec(k));
			}
		}
	}
	return newMat;
}
