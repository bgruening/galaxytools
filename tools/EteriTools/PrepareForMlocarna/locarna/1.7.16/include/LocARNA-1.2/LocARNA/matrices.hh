#ifndef LOCARNA_MATRICES_HH
#define LOCARNA_MATRICES_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* @file Define various generic matrix classes (with templated element
   type): simple matrix, matrix with range restriction, matrix with
   offset, rotatable matrix.
 */

#include <iostream>
#include <vector>
#include <assert.h>

#include <algorithm>

#include "matrix.hh"

namespace LocARNA {

    // ----------------------------------------
    //! \brief Simple matrix class with restriction to a range.
    //!
    //! The matrix features a fix maximal range [0..n]x[0..m].
    //! It can be restricted to [xl..xr]x[yl..yr]. After such
    //! a restriction, the matrix is invalidated and can
    //! only be used with indices (i,j): xl<=i<=xr and yl<=j<=yr  
    //!
    //! @note I planned to use this for the M matrices in LocARNA to optimize locality.
    //! However, I didn't see a performance improvement (maybe for very large instances?)
    //! 
    template <class elem_t>
    class RMatrix : public Matrix<elem_t> {
	typedef typename Matrix<elem_t>::size_type size_type; //!< size type
    protected:
	size_type xl_; //!< left end of restriction in first dimension
	size_type xr_; //!< right end of restriction in first dimension
	size_type yl_; //!< left end of restriction in second dimension
	size_type yr_; //!< right end of restriction in second dimension

	size_type xfactor_; //!< factor for index calculation
	size_type offset_;  //!< offset for index calculation 
    
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
	    assert(xl_<=i && i<=xr_);
	    assert(yl_<=j && j<=yr_);
	
	    return i*xfactor_ + j - offset_;
	}
    
    public:
	
	/** 
	 * Construct as 0x0-matrix
	 * 
	 */
	RMatrix() 
	    : Matrix<elem_t>()
	{}
    
	
	/**
	 * Resize matrix 
	 *
	 * Reserves memory in sufficient size
	 * and removes any existing restrictions
	 * 
	 * @param xdim 
	 * @param ydim 
	 */
	void
	resize(size_type xdim, size_type ydim) {
	    this->xdim_=xdim;
	    this->ydim_=ydim;
	    this->mat_.resize(xdim*ydim);
	
	    restrict(0,xdim-1,0,ydim-1);
	}

	/** 
	 * \brief Set new restricted range
 	 * 
	 * @param xl 
	 * @param xr 
	 * @param yl 
	 * @param yr 
	 */
	void restrict(size_type xl,size_type xr,size_type yl,size_type yr) {
	    assert(xl>=0);
	    assert(yl>=0);
	    assert(xr<this->xdim_);
	    assert(yr<this->ydim_);
	    assert(xl<=xr);
	    assert(yl<=yr);
	
	    this->xl_=xl;
	    this->xr_=xr;
	    this->yl_=yl;
	    this->yr_=yr;

	    this->xfactor_ = (yr-yl)+1;
	    this->offset_  = 0;
	    this->offset_  = addr(xl,yl);
	}
    
	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t & operator() (size_type i,size_type j) const {
	    return this->mat_[addr(i,j)];
	}
    
	/** 
	 * Read/write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return reference to entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	elem_t & operator() (size_type i,size_type j) {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t get(size_type i,size_type j) const {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * @param x element value
	 * @note: redefine since we don't want to use polymorphism
	 */
	void
	set(size_type i,size_type j, const elem_t &x) {
	    this->mat_[addr(i,j)]=x;
	}
	
	/** 
	 * \brief Fill the whole matrix with the given value 
	 * 
	 * @param val value assigned to each entry 
	 * @post all matrix entries are set to val
	 */
	void 
	fill(const elem_t &val) {
	    for (size_type i=xl_; i<xr_; ++i)
		for (size_type j=yl_; j<yr_; ++j)
		    this->mat_(i,j)=val;
	}
    };


    // ----------------------------------------
    //! @brief Simple matrix class with offset
    template <class elem_t>
    class OMatrix : public Matrix<elem_t> {
    protected:
	size_t off_; //!< combined offset for vector access
	size_t xoff_; //!< offset in first dimension
	size_t yoff_; //!< offset in second dimension

	/** 
	 * \brief Computes address/index in 1D vector from 2D matrix indices
	 * 
	 * Shifts indices by offset
	 *
	 * @param i first index
	 * @param j second index
	 * 
	 * @return index in vector
	 * @note this method is used for all internal access to the vector mat_
	 */
	size_t addr(size_t i, size_t j) const {
	    assert(xoff_<=i && i<xoff_+this->xdim_);
	    assert(yoff_<=j && j<yoff_+this->ydim_);
	    return i*this->xdim_ + j - off_;
	}
    
    public:
	
	/** 
	 * Construct as 0x0-matrix
	 * 
	 * @return 
	 */
	OMatrix() : Matrix<elem_t>(0) {
	}
        
	/** 
	 * Resize matrix
	 * 
	 * @param xdim new first dimension
	 * @param ydim new second dimension
	 * @param xoff new offset in first dimension
	 * @param yoff new offset in second dimension
	 */
	void
	resize(size_t xdim, size_t ydim, size_t xoff=0, size_t yoff=0) {
	    xoff_=xoff;
	    yoff_=yoff;
	    off_=xoff*xdim+yoff;
	    this->xdim_=xdim;
	    this->ydim_=ydim;
	    this->mat_.resize(xdim*ydim);
	}

	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t & 
	operator() (size_t i,size_t j) const {
	    return this->mat_[addr(i,j)];
	}
    
	/** 
	 * Read/write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return reference to entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	elem_t & 
	operator() (size_t i,size_t j) {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t get(size_t i,size_t j) const {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * @param x element value
	 * @note: redefine since we don't want to use polymorphism
	 */
	void
	set(size_t i,size_t j, const elem_t &x) {
	    this->mat_[addr(i,j)]=x;
	}
	
    };


    // ----------------------------------------
    //! @brief A matrix class with rotation
    //! 
    template <class elem_t>
    class RotMatrix: public Matrix<elem_t> {

    protected:
	size_t xrot_; //!< rotation in dimension 1
	size_t yrot_; //!< rotation in dimension 2

	/** 
	 * Compute index rotation for one dimension
	 * 
	 * @param x index
	 * @param r amount of rotation
	 * @param d dimension size
	 * 
	 * @return rotated index
	 */
	size_t rot(size_t x, size_t r, size_t d) {
	    assert(r<d);
	    return (x+d-r)%d;
	}
    
	/** 
	 * Computes address/index in 1D vector from 2D matrix indices
	 *
	 * Does the rotation.
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return index in vector
	 * @note this method is used for all internal access to the vector mat_
	 */
	size_t addr(size_t i, size_t j) const {
	    assert(xrot_<=i && i<xrot_+this->xdim_);
	    assert(yrot_<=j && j<yrot_+this->ydim_);
	    return rot(i,xrot_,this->xdim_)*this->xdim_ + rot(j,yrot_,this->ydim_);
	}

    public:    

	/** 
	 * Construct as empty 0x0-matrix. 
	 *  
	 * @return 
	 */
	RotMatrix() : Matrix<elem_t>(0) {
	}
        
	/** 
	 * Resize matrix
	 * 
	 * @param xdim first dimension
	 * @param ydim second dimension
	 * @param xrot rotation in first dimension
	 * @param yrot rotation in second dimension
	 */
	void
	resize(size_t xdim, size_t ydim, size_t xrot=0, size_t yrot=0) {
	    xrot_=xrot;
	    yrot_=yrot;
	    this->xdim_=xdim;
	    this->ydim_=ydim;
	    this->mat_.resize(xdim*ydim);
	}

	/** 
	 * Change rotation
	 * 
	 * @param xrot new first rotation
	 * @param yrot new second rotation
	 * @post matrix is rotated
	 */
	void move(size_t xrot, size_t yrot) {
	    xrot=xrot_;
	    yrot=yrot_;
	}
    
	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t & operator() (size_t i,size_t j) const {
	    return this->mat_[addr(i,j)];
	}
    
	/** 
	 * Read/write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return reference to entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	elem_t & operator() (size_t i,size_t j) {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Read access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * 
	 * @return entry (i,j)
	 * @note: redefine since we don't want to use polymorphism
	 */
	const elem_t get(size_t i,size_t j) const {
	    return this->mat_[addr(i,j)];
	}

	/** 
	 * Write access to matrix element
	 * 
	 * @param i 
	 * @param j 
	 * @param x element value
	 *
	 * @note: redefine since we don't want to use polymorphism
	 */
	void
	set(size_t i,size_t j, const elem_t &x) {
	    this->mat_[addr(i,j)]=x;
	}
	
    };

} // end namespace LocARNA

#endif // LOCARNA_MATRICES_HH
