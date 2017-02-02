#pragma once

#include <valarray>

template<typename VALUE_TYPE>
VALUE_TYPE **new_matrix(size_t n_rows, size_t n_columns)
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
template<typename VALUE_TYPE>
void delete_matrix(VALUE_TYPE **ppMatrix)
{
	delete[] ppMatrix[0];
	delete[] ppMatrix;
}


class vector;

class sparse_matrix
{
private:
	double **_m_values;
	size_t **_m_indices;

	size_t _n_rows;
	size_t _n_columns;
	size_t _n_stored_columns;

private:
	void null_pointers()
	{
		this->_m_values = nullptr;
		this->_m_indices = nullptr;
	}

	void delete_sm();

public:
	enum
	{
		NOT_INDEX = -1
	};

	sparse_matrix()
	{ null_pointers(); }

	sparse_matrix(size_t nRows, size_t nStoredColumns)
	{
		null_pointers();
		this->resize(nRows, nRows, nStoredColumns);
	}

	sparse_matrix(size_t nRows, size_t nColumns, size_t nStoredColumns)
	{
		null_pointers();
		this->resize(nRows, nColumns, nStoredColumns);
	}

	~sparse_matrix()
	{
		delete_sm();
	}

	void resize(size_t nRows, size_t nColumns, size_t nStoredColumns);

	size_t n_rows() const
	{ return _n_rows; }

	size_t n_columns() const
	{ return _n_columns; }

	size_t n_stored_columns() const
	{ return _n_stored_columns; }

	double operator()(int rowInd, int columnInd) const;

	double &at(size_t nI, size_t nJ);

	void del_row(int nRow);

	vector vector_multiply(const vector &vIn) const;

	class iterator
	{
		sparse_matrix &m_refMatrix;
		int m_nRow;
		int m_nNonzero;


	public:
		iterator(sparse_matrix &refMatrix)
				: m_refMatrix(refMatrix), m_nRow(), m_nNonzero()
		{}

		// access interface
		double operator*() const
		{
			return m_refMatrix._m_values[m_nRow][m_nNonzero];
		}

		double &operator*()
		{
			return m_refMatrix._m_values[m_nRow][m_nNonzero];
		}

		int first_index() const
		{
			return m_nRow;
		}

		int second_index() const
		{
			return m_refMatrix._m_indices[m_nRow][m_nNonzero];
		}

		// iterator interface
		void First1()
		{
			m_nRow = 0;
			return;
		}

		bool IsDone1() const
		{
			return m_nRow >= m_refMatrix.n_rows();
		}

		void Next1()
		{
			++m_nRow;
		}

		void First2()
		{
			m_nNonzero = 0;
		}

		void First2(int nRow)
		{
			m_nRow = nRow;
			m_nNonzero = 0;
		}

		bool IsDone2() const
		{
			return m_nNonzero >= m_refMatrix.n_stored_columns() ||
				   sparse_matrix::NOT_INDEX == m_refMatrix._m_indices[m_nRow][m_nNonzero];
		}

		void Next2()
		{
			++m_nNonzero;
		}
	};
}; // end of class sparse_matrix declaration

class vector : public std::valarray<double>
{
	typedef std::valarray<double> base_type;

public:
	explicit vector(size_t nSize = 0)
			: base_type(nSize)
	{}

	vector(size_t nSize, double rDefault)
			: base_type(rDefault, nSize)
	{}

	double operator*(const vector &vVec) const;

	vector &operator*=(double rScale);

	vector &operator=(double rRight)
	{
		base_type::operator=(rRight);
		return *this;
	}

	double norm_inf()
	{
		return sqrt((*this) * (*this));
	}

	void add_scale_vector(const vector &vVec, double rScale);
}; // end of data_vector