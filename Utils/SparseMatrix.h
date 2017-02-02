#pragma once

#include <valarray>

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
	{ return this->_n_rows; }

	size_t n_columns() const
	{ return this->_n_columns; }

	size_t n_stored_columns() const
	{ return this->_n_stored_columns; }

	double operator()(int rowInd, int columnInd) const;

	double &at(size_t nI, size_t nJ);

	vector vector_multiply(const vector &vIn) const;
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
		return this->apply(std::fabs).max();
	}

	void add_scale_vector(const vector &vVec, double rScale);
}; // end of data_vector