#include <iostream>
#include <stdexcept>

#include "SparseMatrix.h"

template <typename VALUE_TYPE>
static VALUE_TYPE **new_matrix(size_t n_rows, size_t n_columns)
{
	VALUE_TYPE **ppMatrix = new VALUE_TYPE *[n_rows];
	ppMatrix[0] = new VALUE_TYPE[n_rows * n_columns];

	for (int i = 1; i < n_rows; i++)
	{
		ppMatrix[i] = ppMatrix[i - 1] + n_columns;
	}

	return ppMatrix;
}

/**
 * @memo    free a matrix
 */
template <typename VALUE_TYPE>
static void delete_matrix(VALUE_TYPE **ppMatrix)
{
	delete[] ppMatrix[0];
	delete[] ppMatrix;
}

void sparse_matrix::resize(size_t nRows, size_t nColumns, size_t nStoredColumns)
{
	this->delete_sm();

	this->_m_values = new_matrix<double>(nRows, nStoredColumns);
	this->_m_indices = new_matrix<size_t>(nRows, nStoredColumns);

	this->_n_rows = nRows;
	this->_n_columns = nColumns;
	this->_n_stored_columns = nStoredColumns;

	for (size_t nI = 0; nI < nRows; nI++)
	{
		for (size_t nK = 0; nK < nStoredColumns; nK++)
		{
			this->_m_values[nI][nK] = 0.0;
			this->_m_indices[nI][nK] = NOT_INDEX;
		}
	}
}

void sparse_matrix::delete_sm()
{
	if (this->_m_values != nullptr)
	{
		delete_matrix(this->_m_values);
	}
	if (this->_m_indices != nullptr)
	{
		delete_matrix(this->_m_indices);
	}
}

double sparse_matrix::operator()(int rowInd, int columnInd) const
{
	if (rowInd > this->_n_rows)
		throw std::out_of_range("Row index is out of range!");

	for (size_t nK = 0; nK < this->_n_stored_columns; nK++)
	{
		if (this->_m_indices[rowInd][nK] == columnInd)
		{
			return this->_m_values[rowInd][nK];
		}
	}

	return 0;
}

double &sparse_matrix::at(size_t nI, size_t nJ)
{
	if (nI > this->n_rows())
		throw std::runtime_error("Index out of range: " + std::to_string(nI));

	for (int nK = 0; nK < this->_n_stored_columns; nK++)
	{
		size_t nSecondIdx = this->_m_indices[nI][nK];
		if (nSecondIdx == NOT_INDEX || nSecondIdx == nJ)
		{
			this->_m_indices[nI][nK] = nJ;
			return this->_m_values[nI][nK];
		}
	}

	throw std::runtime_error("incorrect n_stored_columns of the matrix");
}

vector sparse_matrix::vector_multiply(const vector &vIn) const
{
	vector vOut(vIn.size(), 0);

	for (size_t nI = 0; nI < this->_n_rows; nI++)
	{
		for (size_t nJ = 0; nJ < this->_n_stored_columns; nJ++)
		{
			size_t nIdx = this->_m_indices[nI][nJ];
			if (NOT_INDEX == nIdx)
				break;

			vOut[nI] += this->_m_values[nI][nJ] * vIn[nIdx];
		}
	}

	return vOut;
}

double vector::operator*(const vector &vVec) const
{
	if (vVec.size() != this->size())
		throw std::runtime_error("Wrong vector size: " + std::to_string(vVec.size()));

	double rDotProduct = 0.0;
	for (size_t nI = 0; nI < size(); ++nI)
	{
		rDotProduct += (*this)[nI] * vVec[nI];
	}

	return rDotProduct;
} // end of operator*

vector &vector::operator*=(double rScale)
{
	for (auto &row : *this)
	{
		row *= rScale;
	}

	return *this;
} // end of operator*=

void vector::add_scale_vector(const vector &vVec, double rScale)
{
	for (size_t nI = 0; nI < this->size(); ++nI)
	{
		(*this)[nI] += vVec[nI] * rScale;
	}
}

vector vector::operator-(const vector &vVec) const
{
	vector r(*this);
	for (size_t nI = 0; nI < vVec.size(); ++nI)
	{
		r[nI] -= vVec[nI];
	}
	return r;
}
// end of add_scale_vector






