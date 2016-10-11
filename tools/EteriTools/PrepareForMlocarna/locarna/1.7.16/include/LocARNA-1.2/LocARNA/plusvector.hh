#ifndef LOCARNA_PLUSVECTOR_HH
#define LOCARNA_PLUSVECTOR_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <vector>

namespace LocARNA {
    /**
     * \brief Implements a vector with += operator
     *
     * The class provides an alternate syntax for push_back.
     * 
     * @see Alignment::write_pp() for usage example
     */
    template<class T>
    class plusvector: public std::vector<T> {

    public:
	
	/** 
	 * @brief Add an element to end of vector
	 * 
	 * @param x vector element
	 * 
	 * @return *this
	 * @post x is added at end of *this (like push_back)
	 */
	plusvector& operator += (const T &x) {
	    this->push_back(x);
	    return *this;
	}
    };

}

#endif // LOCARNA_PLUSVECTOR_HH
