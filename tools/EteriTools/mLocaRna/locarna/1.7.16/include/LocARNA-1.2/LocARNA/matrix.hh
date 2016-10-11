#ifndef LOCARNA_MATRIX_HH
#define LOCARNA_MATRIX_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* @file Define simple, generic matrix class (with templated element
   type)
 */

#include <iostream>
#include <vector>
#include <assert.h>

#include <algorithm>

namespace LocARNA {

    /*
      Define classes for the dynamic programming
      matrices.

      The structures support offsets and moving
      of the offset, while maintaining the entries
      in the overlapping sub-matrix
    */

    //! simple 2D matrix class, provides access via operator (int,int) 
    template <class T>
    class Matrix {
    public:
	typedef T elem_t; //!< type of elements
	typedef typename std::vector<elem_t>::size_type size_type; //!< size type (from underlying vector)
	
	typedef std::pair<size_type,size_type> size_pair_type; //!< type for pair of sizes
    
    protected:
	std::vector<elem_t> mat_; //!< vector storing the matrix entries
	size_type xdim_; //!< first dimension
	size_type ydim_; //!< second dimension
    
	/** 
	 * Computes address/index in 1D vector from 2D matrix indices
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return index in vector
	 * @note this method is used for all internal access to the vector mat_
	 */
	size_type addr(size_type i, size_type j) const {
	    assert(0<=i && i<this->xdim_);
	    assert(0<=j && j<this->ydim_);
	    return i*ydim_+j;
	}

    public:
	/** 
	 * Empty constructor
	 * 
	 */
	Matrix() 
	    : mat_(),xdim_(0),ydim_(0) {
	}
    
	/** 
	 * Construct with dimensions, optionally initialize from array  
	 * 
	 * @param xdim first dimension of matrix
	 * @param ydim second dimension of matrix
	 * @param from pointer to array of elements
	 *
	 * @note if from given and !=0 initialize from array from
	 *
	 */
	Matrix(size_type xdim, size_type ydim, const elem_t *from=0L)
	    : mat_(xdim*ydim),xdim_(xdim),ydim_(ydim) {
	    if (from!=0L) {
		for (size_type i=0; i<xdim_; i++) {
		    for (size_type j=0; j<ydim_; j++) {
			(*this)(i,j)=from[i*ydim+j];
		    }
		}
	    }
	}
    
	/** 
	 * Access size
	 * 
	 * @return size of matrix as pair of dimensions
	 */
	size_pair_type sizes() const {
	    return size_pair_type(xdim_,ydim_);
	}

	/** 
	 * Resize both dimensions
	 *
	 * @param xdim first dimension
	 * @param ydim second dimension
	 */
	void
	resize(size_type xdim, size_type ydim) {
	    xdim_=xdim;
	    ydim_=ydim;
	
	    mat_.resize(xdim_*ydim_);
	}
    
	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 */
	const elem_t & operator() (size_type i,size_type j) const {
	    return mat_[addr(i,j)];
	}
    
	/** 
	 * Read/write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return reference to entry (i,j)
	 */
	elem_t & operator() (size_type i,size_type j) {
	    return mat_[addr(i,j)];
	}

	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 */
	const elem_t get(size_type i,size_type j) const {
	    return mat_[addr(i,j)];
	}

	/** 
	 * Write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * @param x element value
	 * 
	 */
	void
	set(size_type i,size_type j, const elem_t &x) {
	    mat_[addr(i,j)]=x;
	}

	/** 
	 * \brief Fill the whole matrix with the given value 
	 * 
	 * @param val value assigned to each entry 
	 * @post all matrix entries are set to val
	 */
	void 
	fill(const elem_t &val) {
	    for (size_type i=0; i<xdim_*ydim_; ++i)
		mat_[i]=val;
	}

	/** 
	 * Clear the matrix
	 * @post the matrix is resized to dimensions (0,0) 
	 * @note behaves like std::vector::clear()
	 */
	void
	clear() {
	    resize(0,0);
	    mat_.clear();
	}
    
	/** 
	 * Transform matrix in place due to applying a given function to each element
	 * 
	 * @param f function object
	 * 
	 * @post All matrix entries are changed from x to f(x)
	 * @note applies f via in place std::transform to all matrix entries
	 */
	template<class UnaryOperator>
	void transform(UnaryOperator f) {
	    std::transform(mat_.begin(),mat_.end(),mat_.begin(),f);
	}
    };

    /** 
     * Output operator for writing (templated) matrix to output stream
     * 
     * @param out the output stream
     * @param mat the matrix to be written
     * 
     * @return output stream after writing matrix mat
     */
    template <class T>
    std::ostream & operator << (std::ostream &out, Matrix<T> mat) {
	typename Matrix<T>::size_pair_type sizes = mat.sizes();
    
	for (typename Matrix<T>::size_type i=0; i<sizes.first; i++) {
	    for (typename Matrix<T>::size_type j=0; j<sizes.second; j++) {
		out << mat(i,j) << " ";
	    }
	    out << std::endl;
	}
	return out;
    }

    /** 
     * Input operator for reading (templated) matrix from input stream
     * 
     * @param in the input stream
     * @param[out] mat the matrix to be read
     * 
     * @return input stream after reading matrix mat
     */
    template <class T>
    std::istream & operator >> (std::istream &in, Matrix<T> &mat) {
	typename Matrix<T>::size_pair_type sizes = mat.sizes();
	for (typename Matrix<T>::size_type i=0; i<=mat.sizes().first; i++) {
	    for (typename Matrix<T>::size_type j=0; j<=mat.sizes().second; j++) {
		in >> mat(i,j);
	    }
	}
	return in;
    }



} // end namespace LocARNA

#endif // LOCARNA_MATRIX_HH
