#ifndef LOCARNA_TUPLES_HH
#define LOCARNA_TUPLES_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

namespace LocARNA {

    /**
     * @brief Represents a 3-tuple
     *
     * triple stores three values first, second, third.
     * extension of std::pair to 3-tuple
     */
    template<class T1,class T2,class T3>
    class triple: public std::pair<T1,T2> {
    public:
	T3 third; //!< third value
	
	/** 
	 * Construct from three values 
	 * 
	 * @param x1 value 1
	 * @param x2 value 2
	 * @param x3 value 3
	 * 
	 */
	triple(const T1 &x1,const T2 &x2,const T3 &x3): std::pair<T1,T2>(x1,x2),third(x3) {
	}
    };
    
    /**
     * @brief Represents a 4-tuple
     *
     * quadruple stores four values first, second, third, fourth.
     * extension of triple to 4-tuple
     */
    template<class T1,class T2,class T3,class T4>
    class quadruple: public triple<T1,T2,T3> {
    public:
	T4 fourth; //!< fourth value
	
	/** 
	 * \brief Construct from four values 
	 * 
	 * @param x1 value 1
	 * @param x2 value 2
	 * @param x3 value 3
	 * @param x4 value 4
	 * 
	 */
	quadruple(const T1 &x1,const T2 &x2,const T3 &x3,const T4 &x4): triple<T1,T2,T3>(x1,x2,x3),fourth(x4) {
	}
    };

    /**
     * @brief Represents a 5-tuple
     *
     * quintuple stores five values first, second, third, fourth, fifth
     * extension of triple to 4-tuple
     */
    template<class T1,class T2,class T3,class T4,class T5>
    class quintuple: public quadruple<T1,T2,T3,T4> {
    public:
	T5 fifth; //!< fifth value

	/**
	 * \brief Construct from five values
	 *
	 * @param x1 value 1
	 * @param x2 value 2
	 * @param x3 value 3
	 * @param x4 value 4
	 * @param x5 value 5
	 *
	 */
	quintuple(const T1 &x1,
		  const T2 &x2,
		  const T3 &x3,
		  const T4 &x4,
		  const T5 &x5)
	    : quadruple<T1,T2,T3,T4>(x1,x2,x3,x4),
	      fifth(x5) {
	}
	
    };

} // end namespace LocARNA

#endif // LOCARNA_TUPLES_HH
