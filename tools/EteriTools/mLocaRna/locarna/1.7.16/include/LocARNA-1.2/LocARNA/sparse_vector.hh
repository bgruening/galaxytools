#ifndef SPARSE_VECTOR_HH
#define SPARSE_VECTOR_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>

#include <tr1/unordered_map>

namespace LocARNA {
    
    /**
     * \brief Represents a sparse vector
     *
     * Sparse vector of entries val_t implements the vector by a hash
     * map. The class is designed to be largely exchangable with
     * non-sparse vectors.
     *
     * @todo this code is basically a stripped down version of the SparseMatrix code;
     * likely one could reduce redundancy
     */
    template <typename T>
    class SparseVector {
    public:
	
  	typedef T value_t; //!< type of vector entries

	typedef size_t size_type; //!< usual definition of size_type

	typedef size_type key_t; //!< type of vector index pair
	
    protected:
		
	typedef std::tr1::unordered_map<key_t,value_t> map_t; //!< map type  
	map_t the_map_; //!< internal representation of sparse vector
	value_t def_; //!< default value of vector entries
    
    public:

	/**
	 * \brief Stl-compatible constant iterator over vector elements.
	 *
	 * Behaves like a const iterator of the hash map.
	 */
	typedef typename map_t::const_iterator const_iterator;
	
	
	/**
	 * \brief Element of sparse vector 
	 * 
	 * Proxy for sparse vector entries. This is required for
	 * non-const access to vector elements in order to
	 * provide a very similar syntax for the sparse data
	 * structure and the corresponding non-sparse vector.
	 */
	class element {
	private:
	    SparseVector<T> *m_;
	    key_t k_;
	public:
	    /** 
	     * @brief Construct as proxy for specified element in given sparse vector
	     * 
	     * @param m pointer to sparse vector
	     * @param k key/index of entry in given sparse vector
	     *
	     */
	    element(SparseVector<T> *m,key_t k): m_(m),k_(k) {}

	    /** 
	     * @brief Access entry for which the class acts as proxy
	     * 
	     * @return value of vector entry.
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
	     * @post x is added to vector entry for that the class is proxy
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
	     * @post x is assigned to vector entry for that the class is proxy
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
	SparseVector(const value_t &def) : the_map_(),def_(def) {}

	/** 
	 * \brief Access to vector element
	 * 
	 * @param i index first dimension
	 *
	 * @return proxy to vector entry i
	 */
	element operator[] (size_type i) {
	    return element(this,key_t(i));
	}
    
	/** 
	 * \brief Read-only access to vector element of const vector
	 * 
	 * @param i index first dimension
	 *
	 * @return vector entry i
	 */
	const value_t & operator[] (size_type i) const {
	    const_iterator it = the_map_.find(key_t(i));
	    if ( it == the_map_.end() ) 
		return def_;
	    else 
		return it->second;
	}
	
	/** 
	 * @brief Write access to vector entry
	 * 
	 * @param i index first dimension
	 * @param val value to be written to entry i
	 *
	 * @note Unlike the assignment operator (via element), there is no
	 * test whether the default value is assigned.
	 * Use reset(i) if you want to reset vector entries to the default. 
	 *
	 * @post writes entry. If entry didn't exist already it is created.
	 */
	void
	set(size_type i, const value_t &val) {
	    typename map_t::iterator it = the_map_.find(key_t(i));
	    if ( it != the_map_.end() ) { 
		it->second = val;
	    } else { 
		the_map_.insert(typename map_t::value_type(key_t(i),val));
	    }
	}
    
	/** 
	 * @brief Set vector entry to default value
	 * 
	 * @param i index first dimension
	 */
	void
	reset(size_type i) {
	    typename map_t::iterator it = the_map_.find(key_t(i));
	    if ( it != the_map_.end() ) { 
		the_map_.erase(key_t(i));
	    } 
	}

	/**
	 * @brief Size of sparse vector
	 * @return number of non-empty entries
	 */
	size_type
	size() const {
	    return the_map_.size();
	}

	/**
	 * @brief Check for emptiness
	 * @return true, if sparse vector contains 
	 * only implicite default entries.
	 */
	bool
	empty() const {
	    return the_map_.empty();
	}
	
	/** 
	 * @brief Clear the vector
	 */
	void
	clear() {
	    the_map_.clear();
	}

	/** 
	 * \brief Begin const iterator over vector entries
	 * 
	 * @return const iterator pointing to begin of entry hash
	 *
	 * @see end()
	 */
	const_iterator begin() const {
	    return the_map_.begin();
	}
	
	/** 
	 * \brief End const iterator over vector entries
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
     * @param m sparse vector to be writing to stream
     * 
     * @return output stream after writing
     */
    template<class T>
    inline
    std::ostream &
    operator <<(std::ostream &out, const SparseVector<T> &v) {
	for (typename SparseVector<T>::const_iterator it=v.begin();
	     v.end()!=it; 
	     ++it) {
	    out <<it->first <<":" << it->second << " ";
	}
	out << std::endl;
	return out;
    }

} //end namespace LocARNA

#endif // SPARSE_VECTOR_HH
