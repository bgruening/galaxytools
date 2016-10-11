#ifndef LOCARNA_INFTY_INT_HH
#define LOCARNA_INFTY_INT_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <algorithm>
#include <iosfwd>
#include <assert.h>

namespace LocARNA {

    class FiniteInt;
    class InftyInt;
    
    /**
       Potentially infinite value where finite values are integer.  In
       contrast to class NormInfyInt, the representation of infinity
       can be denormalized, this means that it is not safe (and
       therefore not supported) to subtract or add a second
       potentially infinite value.  In general, such non-normalized
       potentially infinite integers are produced as results of
       adding/subtracting normal potentially infinite integers.

       The infinity integer classes TaintedInftyInt, InftyInt and
       FiniteInt support efficient addition and
       minimization/maximization of potentially infinite integer
       values. For this purpose, the range of the base type (long
       int) is restricted. With this encoding, two normal infinite
       values (InftyInt) of the same sign can be added without
       resulting in overflow. This yields a non-normal
       TaintedInftyInt, which still allows to add finite integers
       (FiniteInt) without overflow (at least as long as the added
       finite integers do not exceed the remaining range of
       [-m..m-1]), where m=2^(s-1) and s is the width of the base
       type, e.g. 64.

       The range of the base type is split into subranges, where s is
       the number of bits in the base type: 
       
       * negative infinity        [ -m..-m/5 [
       * normal negative infinity [ -3m/5..-m/5 [
       * normalized negative infinity -2m/5
       * finite range             [ -m/5 .. m/5 [
       * normal positive infinity [ m/5 .. 3m/5 [
       * normalized positive infinity 2m/5
       * positive infinity        [ m/5 .. m [
       
       Normalizing resets negative infinite values to -2m/5 and
       positive infinite values to 2m/5. In this way, adding two
       normalized infinite values of the same sign and some finite
       value using fast standard integer addition always works without
       overflow and results in a defined value in the infinite range.
    */
    class TaintedInftyInt {
    public:
	//! the base type
	typedef long int base_type;

    protected:
	base_type val; //!< value
	
	//! minimum finite value
	static const base_type min_finity;

	//! maximum finite value
	static const base_type max_finity;
	
	//! minimum normal infinite value
	static const base_type min_normal_neg_infty;
	
	//! maximum normal infinite value
	static const base_type max_normal_pos_infty;
    public:
	
	/** 
	 * @brief Construct empty
	 */
	TaintedInftyInt(): val(0) {
	}

	/** 
	 * @brief Construct from base type
	 * @param x base type value
	 */
	explicit
	TaintedInftyInt(const base_type &x)
	    : val(x) {
	}
	
	/** 
	 * @brief minimum finite value
	 * 
	 * @return minimum finite value
	 */
	static
	base_type
	min_finite() {
	    return min_finity;
	}

	/** 
	 * @brief maximum finite value
	 * 
	 * @return maximum finite value
	 */
	static
	base_type
	max_finite() {
	    return max_finity;
	}
	
	/** 
	 * Test for negative infinity 
	 * 
	 * @return whether object represents negative infinity 
	 */
	bool
	is_neg_infty() const {
	    return val < min_finity;
	}

	/** 
	 * Test for positive infinity 
	 * 
	 * @return whether object represents positive infinity 
	 */
	bool
	is_pos_infty() const {
	    return val > max_finity;
	}

	/** 
	 * Test for finite 
	 * 
	 * @return whether object has a finite value
	 */
	bool
	is_finite() const {
	    return min_finity <= val &&  val <= max_finity;
	}
	
	/** 
	 * Test for finite or normal infinity 
	 * 
	 * @return whether object is in the range of normal infinity (or finite)
	 */
	bool 
	is_normal() const {
	    return min_normal_neg_infty <= val &&  val <= max_normal_pos_infty;
	}
	
	/** 
	 * @brief Convert finite value to base type
	 * 
	 * @return value
	 * @pre is finite
	 */
	base_type 
	finite_value() const {
	    assert(is_finite());
	    return val;
	}

	/** 
	 * @brief Assignment
	 * 
	 * @param x finite int to be assigned
	 * 
	 * @return *this
	 */
	TaintedInftyInt &
	operator =(const FiniteInt &x);

	/** 
	 * @brief Equality test
	 * 
	 * @param x operand 1 (tainted)
	 * @param y operand 2 (tainted)
	 * 
	 * @return whether x equals y
	 */
	friend
	bool
	operator ==(const TaintedInftyInt &x, const TaintedInftyInt &y);
	
	/** 
	 * @brief Add
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x plus y
	 */
	friend
	TaintedInftyInt
	operator +(const TaintedInftyInt &x, const FiniteInt &y);
	
	/** 
	 * @brief Subtract
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x minus y
	 */
	friend
	TaintedInftyInt
	operator -(const TaintedInftyInt &x, const FiniteInt &y);

	/** 
	 * @brief Add
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x plus y
	 */
	friend
	TaintedInftyInt
	operator +(const InftyInt &x, const InftyInt &y);
	
	/** 
	 * @brief Subtract
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x minus y
	 */
	friend 
	TaintedInftyInt
	operator -(const InftyInt &x, const InftyInt &y);
	
	/** 
	 * @brief Minimum
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return min of x and y
	 */
	friend
	TaintedInftyInt
	min(const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * @brief Maximum
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return max of x and y
	 */
	friend
	TaintedInftyInt
	max(const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Greater than operator
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return whether x is greater than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator > (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Less than operator
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return whether x is less than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator < (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Greater or equal than operator
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return whether x is greater or equal than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator >= (const TaintedInftyInt &x, const TaintedInftyInt &y);

	/** 
	 * Less or equal than operator
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return whether x is less or equal than y
	 *
	 * @note The result is undefined when comparing infinte values
	 * of the same sign.
	 */ 
	friend
	bool
	operator <= (const TaintedInftyInt &x, const TaintedInftyInt &y);


	friend class InftyInt;

	/** 
	 * Write TaintedInftyInt object to stream
	 * 
	 * @param out output stream
	 * @param x object
	 * 
	 * @return output stream after writing
	 */
	friend
	std::ostream &
	operator <<(std::ostream &out, const TaintedInftyInt &x);

    };
    
    /**
       Potentially infinite value where finite values are integer.
       The representation of infinite values is normalized.  Due to
       the normalization it is safe to add or subtract a second
       normal potentially infinite integer (thereby generating an
       non-normalized potentially infinite integer).
    */
    class InftyInt : public TaintedInftyInt {
	
    public:
	
	//! normalized negative infinity
	static const InftyInt neg_infty;
	//! normalized positive infinity
	static const InftyInt pos_infty;
	
    private:
	void normalize() {
	    //std::cout << "NORMALIZE" <<std::endl;
	    if (is_neg_infty()) {
		val = neg_infty.val;
	    } else if (is_pos_infty()) {
		val = pos_infty.val;
	    }
	}
    public:
	
	/** 
	 * @brief Construct empty
	 * 
	 */
	InftyInt(): TaintedInftyInt() {
	}

	/** 
	 * @brief Construct from base type
	 *
	 * @param x value of base type 
	 */
	explicit
	InftyInt(const base_type &x):TaintedInftyInt(x) {
	    assert(is_normal());
	}
	
	/** 
	 * @brief Construct from finite int
	 *
	 * @param x value
	 */
	InftyInt(const FiniteInt &x);
	
	/** 
	 * @brief Construct from potentially tainted
	 * 
	 * @param x value
	 */
	InftyInt(const TaintedInftyInt &x) :TaintedInftyInt(x) {
	    normalize();
	}

	/** 
	 * @brief Assignment from potentially tainted infty int
	 * 
	 * @param x value
	 * 
	 * @return *this
	 */
	InftyInt &
	operator =(TaintedInftyInt &x) {
	    val = x.val;
	    normalize();
	    return *this;
	}

	/** 
	 * Add in place
	 * 
	 * @param x operand
	 * 
	 * @return *this after operation
	 */
	InftyInt &
	operator +=(const FiniteInt &x);

	/** 
	 * Subtract in place
	 * 
	 * @param x operand
	 * 
	 * @return *this after operation
	 */
	InftyInt &
	operator -=(const FiniteInt &x);


	/** 
	 * Add
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x plus y
	 */
	friend
	InftyInt
	operator +(const InftyInt &x, const FiniteInt &y);
	
	/** 
	 * Subtract
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * 
	 * @return x minus y
	 */
	friend
	InftyInt
	operator -(const InftyInt &x, const FiniteInt &y);
	
    };
    
    /**
       Finite integer value compatible with potentially infinite
       integers.
    */
    class FiniteInt : public InftyInt {
    public:

	/** 
	 * @brief Construct empty
	 */
	FiniteInt(): InftyInt() {
	}
	
	/** 
	 * @brief Construct from base type value
	 * @param x value
	 */
	FiniteInt(base_type x): InftyInt(x) {
	    assert(is_finite());
	}
	
	/** 
	 * @brief Access finite value 
	 * 
	 * @return value
	 */
	const base_type &
	finite_value() const {
	    return val;
	}
	
	/** 
	 * @brief Add
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * @return x plus y
	 */
	friend
	FiniteInt
	operator +(const FiniteInt &x, const FiniteInt &y);

	/** 
	 * @brief Subtract
	 * 
	 * @param x operand 1
	 * @param y operand 2
	 * @return x minus y
	 */
	friend
	FiniteInt
	operator -(const FiniteInt &x, const FiniteInt &y);
	
    };
    
    inline
    TaintedInftyInt
    operator +(const TaintedInftyInt &x, const FiniteInt &y) {
	TaintedInftyInt res(x);
	res.val+=y.val;
	return res;
    }
    
    inline
    TaintedInftyInt
    operator -(const TaintedInftyInt &x, const FiniteInt &y) {
	TaintedInftyInt res(x);
	res.val-=y.val;
	return res;
    }
    
    inline
    TaintedInftyInt
    operator +(const InftyInt &x, const InftyInt &y) {
	TaintedInftyInt res(x);
	res.val+=y.val;
	return res;    
    }
    
    inline
    TaintedInftyInt
    operator -(const InftyInt &x, const InftyInt &y) {
	TaintedInftyInt res(x);
	res.val-=y.val;
	return res;
    }


    inline
    InftyInt &
    InftyInt::operator +=(const FiniteInt &x) {
	val += x.val;
	return *this;
    }
    
    inline
    InftyInt &
    InftyInt::operator -=(const FiniteInt &x) {
	val -= x.val;
	return *this;
    }
    
    
    inline
    InftyInt
    operator +(const InftyInt &x, const FiniteInt &y) {
	InftyInt res(x);
	res.val+=y.val;
	return res;    
    }
    
    inline
    InftyInt
    operator -(const InftyInt &x, const FiniteInt &y) {
	InftyInt res(x);
	res.val-=y.val;
	return res;
    }

    inline
    FiniteInt
    operator +(const FiniteInt &x, const FiniteInt &y) {
	FiniteInt res(x);
	res.val+=y.val;
	return res;
    }

    inline
    FiniteInt
    operator -(const FiniteInt &x, const FiniteInt &y) {
	FiniteInt res(x);
	res.val-=y.val;
	return res;
    }

    inline
    bool
    operator ==(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val==y.val;
    }

    inline
    TaintedInftyInt &
    TaintedInftyInt::operator =(const FiniteInt &x) {
	val=x.val;
	return *this;
    }

    inline
    TaintedInftyInt
    min(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return TaintedInftyInt(std::min(x.val,y.val));
    }
    
    inline
    TaintedInftyInt
    max(const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return TaintedInftyInt(std::max(x.val,y.val));
    }

    inline
    InftyInt::InftyInt(const FiniteInt &x): TaintedInftyInt(x) {
    }
    
    inline
    bool
    operator > (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val > y.val;
    }

    inline
    bool
    operator < (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val < y.val;
    }

    inline
    bool
    operator >= (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val >= y.val;
    }

    inline
    bool
    operator <= (const TaintedInftyInt &x, const TaintedInftyInt &y) {
	return x.val <= y.val;
    }




} // end namespace LocARNA


#endif
