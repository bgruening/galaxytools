#ifndef LOCARNA_TYPE_WRAPPER_HH
#define LOCARNA_TYPE_WRAPPER_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

namespace LocARNA {
    /**
     * @brief generic type_wrapper class
     * 
     * wrap a comparable type with +/- operations
     *
     */
    template<class T>
    class type_wrapper {
    	T val_;
    public:

	/**
	 * @brief default constructor
	 */
	type_wrapper(): val_() {}
	
	/**
	 * @brief conversion constructor
	 */
	explicit type_wrapper(const T &i): val_(i) {}	
	
	/**
	 * @brief copy constructor
	 */
	type_wrapper(const type_wrapper &idx): val_(idx.val()) {}

	/**
	 * @brief casting to T
	 *
	 * @note We don't use explicit cast for compatibility, since this is
	 * available with c++0x only:
	 * explicit operator T() const {return val_;}
	 * instead we define this "named cast"
	 */
	const T &val() const {return val_;}

	/** 
	 * @brief equal
	 * @param x 
	 * @return *this == x
	 */
	bool operator ==(const type_wrapper &x) const {return val_ == x.val_;}

	
	/** 
	 * @brief inequal 
	 * @param x 
	 * @return *this != x
	 */
	bool operator !=(const type_wrapper &x) const {return val_ != x.val_;}

	/** 
	 * @brief less equal
	 * @param x 
	 * @return *this <= x
	 */
	bool operator <=(const type_wrapper &x) const {return val_ <= x.val_;}

	/** 
	 * @brief less
	 * @param x 
	 * @return *this < x
	 */
	bool operator  <(const type_wrapper &x) const {return val_ <  x.val_;}

	/** 
	 * @brief greater equal
	 * @param x 
	 * @return *this >= x
	 */
	bool operator >=(const type_wrapper &x) const {return val_ >= x.val_;}

	/** 
	 * @brief greater
	 * @param x 
	 * @return *this > x
	 */
	bool operator  >(const type_wrapper &x) const {return val_ >  x.val_;}
	
	/** 
	 * @brief add
	 * @param x 
	 * @return *this + x
	 */
	type_wrapper operator +(const type_wrapper &x) const {return type_wrapper(val_ + x.val_);}
	
	/** 
	 * @brief subtract
	 * @param x 
	 * @return *this - x
	 */
	type_wrapper operator -(const type_wrapper &x) const {return type_wrapper(val_ - x.val_);}
	
	/** 
	 * @brief add
	 * @param x 
	 * @return *this + x
	 */
	type_wrapper operator +(const T &x) const {return type_wrapper(val_ + x);}
	
	/** 
	 * @brief subtract
	 * @param x 
	 * @return *this - x
	 */
	type_wrapper operator -(const T &x) const {return type_wrapper(val_ - x);}
	
	/** 
	 * @brief prefix increment
	 * @return *this after increment
	 */
	const type_wrapper &operator ++() {++val_; return *this;}
	
	/** 
	 * @brief prefix decrement
	 * @return *this after decrement
	 */
	const type_wrapper &operator --() {--val_; return *this;}

	/** 
	 * @brief postfix increment
	 * @return *this before increment
	 */
	const type_wrapper &operator ++(int) {type_wrapper<T> tmp(*this); ++val_; return tmp;}
	
	/** 
	 * @brief postfix decrement
	 * @return *this before decrement
	 */
	const type_wrapper &operator --(int) {type_wrapper<T> tmp(*this); --val_; return tmp;}

    }; // end class type_wrapper

    template <class T>
    std::ostream & operator << (std::ostream &out,const type_wrapper<T> &x) {
	out<<x.val();
	return out;
    }

} // end namespace LocARNA

#endif // LOCARNA_TYPE_WRAPPER_HH
