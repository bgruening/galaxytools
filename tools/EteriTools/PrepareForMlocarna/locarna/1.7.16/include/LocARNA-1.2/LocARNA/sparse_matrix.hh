#ifndef SPARSE_MATRIX_HH
#define SPARSE_MATRIX_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iostream>
#include <tr1/unordered_map>

#include "aux.hh"

namespace LocARNA {

    /**
     * \brief Represents a sparse 2D matrix
     *
     * Sparse matrix of entries val_t implements the matrix by a hash
     * map. The class is designed to be largely exchangable with the
     * non-sparse Matrix class. (A proxy class is used to provide the
     * same syntax for the interface.)
     *
     * @see Matrix
     */

    
    template <typename T>
    class SparseMatrix {
    public:
	
  	typedef T value_t; //!< type of matrix entries

	typedef size_t size_type; //!< usual definition of size_type

	typedef std::pair<size_type,size_type> key_t; //!< type of matrix index pair

    protected:
		
	typedef std::tr1::unordered_map<key_t,value_t,pair_of_size_t_hash > map_t; //!<map type 
	map_t the_map_; //!< internal representation of sparse matrix
	value_t def_; //!< default value of matrix entries
    
    public:

	/**
	 * \brief Stl-compatible constant iterator over matrix elements.
	 *
	 * Behaves like a const iterator of the hash map.
	 */
	typedef typename map_t::const_iterator const_iterator;
	
	
	/**
	 * \brief Element of sparse matrix 
	 * 
	 * Proxy for sparse matrix entries. This is required for
	 * non-const access to matrix elements in order to
	 * provide a very similar syntax for the sparse data
	 * structure and the corresponding non-sparse matrix.
	 */
	class element {
	private:
	    SparseMatrix<T> *m_;
	    key_t k_;
	public:
	    /** 
	     * @brief Construct as proxy for specified element in given sparse matrix
	     * 
	     * @param m pointer to sparse matrix
	     * @param k key/index of entry in given sparse matrix
	     *
	     */
	    element(SparseMatrix<T> *m,key_t k): m_(m),k_(k) {}

	    /** 
	     * @brief Access entry for which the class acts as proxy
	     * 
	     * @return value of matrix entry.
	     *
	     * If entry does not exist, return the default value
	     */
	    operator value_t() {
		typename map_t::const_iterator it = m_->the_map_.find(k_);
		if ( it == m_->the_map_.end() ) 
		    return m_->def_;
		else 
		    return it->second;
	    }
	    
	    /** 
	     * @brief Operator for in place addition 
	     * 
	     * @param x value
	     * 
	     * @post x is added to matrix entry for that the class is proxy
	     *
	     * @return *this after adding x
	     *
	     * @note If entry does not exist, x is added to the default value
	     */
	    element
	    operator +=(const value_t &x) {
		const_iterator it = m_->the_map_.find(k_);
		if ( it == m_->the_map_.end() ) 
		    m_->the_map_[k_] = m_->def_ + x;
		else 
		    m_->the_map_[k_] += x;
	    
		return *this;
	    }
	
	    /** 
	     * @brief Assignment operator
	     * 
	     * @param x value
	     *
	     * @post x is assigned to matrix entry for that the class is proxy
	     * 
	     * @return *this after assigning x
	     *
	     * @note If x equals the default value and the entry exists, it is erased
	     */
	    element &
	    operator =(const value_t &x) {
		if (x==m_->def_) {
		    m_->the_map_.erase(k_);
		} else {
		    // the following replaces m_->the_map_[k_] = x;
		    // but never calls the default constructor for value_t

		    typename map_t::iterator it = m_->the_map_.find(k_);
		    if ( it != m_->the_map_.end() ) { 
			it->second = x;
		    } else { 
			m_->the_map_.insert(typename map_t::value_type(k_,x));
		    }
		}
		return *this;
	    }
	};
    
    
	/** 
	 * @brief Construct with default value
	 * 
	 * @param def default value of entries
	 */
	SparseMatrix(const value_t &def) : the_map_(),def_(def) {}

	/** 
	 * \brief Access to matrix element
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 *
	 * @return proxy to matrix entry (i,j)
	 */
	element operator() (size_type i, size_type j) {
	    return element(this,key_t(i,j));
	}
    
	/** 
	 * \brief Read-only access to matrix element of const matrix
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 *
	 * @return matrix entry (i,j)
	 */
	const value_t & operator() (size_type i, size_type j) const {
	    const_iterator it = the_map_.find(key_t(i,j));
	    if ( it == the_map_.end() ) 
		return def_;
	    else 
		return it->second;
	}
	
	/** 
	 * @brief Write access to matrix entry
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 * @param val value to be written to entry (i,j)
	 *
	 * @note Unlike the assignment operator (via element), there is no
	 * test whether the default value is assigned.
	 * Use reset(i,j) if you want to reset matrix entries to the default. 
	 *
	 * @post writes entry. If entry didn't exist already it is created.
	 */
	void
	set(size_type i, size_type j, const value_t &val) {
	    typename map_t::iterator it = the_map_.find(key_t(i,j));
	    if ( it != the_map_.end() ) { 
		it->second = val;
	    } else { 
		the_map_.insert(typename map_t::value_type(key_t(i,j),val));
	    }
	}

	/** 
	 * \brief Write access to matrix element of const matrix
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 *
	 * @note Creates the entry (i,j) if it is not represented yet.
	 *
	 * @return reference to matrix entry (i,j)
	 */
	value_t & 
	ref(size_type i, size_type j) {
	    typename map_t::iterator it = the_map_.find(key_t(i,j));
	    if ( it == the_map_.end() ) { 
		the_map_.insert(typename map_t::value_type(key_t(i,j),def_));
		it = the_map_.find(key_t(i,j));
	    }
	    return it->second;
	}
    
	/** 
	 * @brief Set matrix entry to default value
	 * 
	 * @param i index first dimension
	 * @param j index second dimension
	 */
	void
	reset(size_type i, size_type j) {
	    typename map_t::iterator it = the_map_.find(key_t(i,j));
	    if ( it != the_map_.end() ) { 
		the_map_.erase(key_t(i,j));
	    } 
	}

	/**
	 * @brief Size of sparse matrix
	 * @return number of non-empty entries
	 */
	size_type
	size() const {
	    return the_map_.size();
	}

	/**
	 * @brief Check for emptiness
	 * @return true, if sparse matrix contains 
	 * only implicite default entries.
	 */
	bool
	empty() const {
	    return the_map_.empty();
	}
	
	/** 
	 * @brief Clear the matrix
	 */
	void
	clear() {
	    the_map_.clear();
	}

	/** 
	 * \brief Begin const iterator over matrix entries
	 * 
	 * @return const iterator pointing to begin of entry hash
	 *
	 * @see end()
	 */
	const_iterator begin() const {
	    return the_map_.begin();
	}
	
	/** 
	 * \brief End const iterator over matrix entries
	 * 
	 * @return const iterator pointing after end of entry hash
	 * @see begin()
	 */
	const_iterator end() const {
	    return the_map_.end();
	}
	
    };

    /** 
     * @brief Output operator
     * 
     * @param out output stream
     * @param m sparse matrix to be writing to stream
     * 
     * @return output stream after writing
     */
    template<class T>
    inline
    std::ostream &
    operator <<(std::ostream &out, const SparseMatrix<T> &m) {
	for (typename SparseMatrix<T>::const_iterator it=m.begin();
	     m.end()!=it; 
	     ++it) {
	    out << "("<<it->first.first<<","<<it->first.second << ") " << it->second << std::endl;
	}
	return out;
    }

} //end namespace LocARNA

#endif // SPARSE_MATRIX_HH
