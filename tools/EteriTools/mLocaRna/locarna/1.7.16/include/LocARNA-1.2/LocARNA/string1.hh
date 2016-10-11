#ifndef LOCARNA_STRING1_HH
#define LOCARNA_STRING1_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string>
#include <algorithm>


namespace LocARNA {
    /** 
     * @brief A simple 1-based string
     *
     * Features:
     * - based on  c++ std::string class, but offers very limited interface
     * - conversion from and to std::string
     * - access via operator []
     * - support of various "string" methods
     */
    class string1 {
	std::string s_;
    public:
	
	/** 
	 * \brief Construct from std::string
	 * 
	 * @param s string
	 */
	string1(const std::string &s): s_(s) {}
	
	/** 
	 * \brief Copy constructor
	 * 
	 * @param s string (of type string1)
	 */
	string1(const string1 &s): s_(s.s_) {}
	
	/** 
	 * Cast to std::string (const)
	 *
	 * @return string
	 */
	const std::string & 
	str() const {
	    return s_;
	}

	/** 
	 * Cast to std::string
	 *
	 * @return string
	 */
	std::string & 
	str() {
	    return s_;
	}
    
	/** 
	 * \brief Read access
	 * 
	 * @param i index F
	 * 
	 * @return ith character of string
	 * @note 1-based
	 */
	const char& 
	operator [](size_t i) const {
	    assert(i>=1);
	    assert(i<=s_.length());
	    return s_[i-1];
	}
	
	/** 
	 * \brief Read/write access
	 * 
	 * @param i index 
	 * 
	 * @return (reference to) ith character of string
	 * @note 1-based
	 */
	char& 
	operator [](size_t i) {
	    assert(i>=1);
	    assert(i<=s_.length());
	    return s_[i-1];
	}
    
	/** 
	 * \brief Provide length
	 * 
	 * @return length of string
	 */
	size_t
	length() const {
	    return s_.length();
	}
    
	/** 
	 * \brief Assignment operator
	 * 
	 * @param s string
	 * 
	 * @return *this
	 * @post *this equals s
	 */
	const string1 &
	operator =(const string1 &s) {s_ = s.s_; return *this;}

	/** 
	 * \brief reverse string
	 * 
	 */
	void
	reverse() {
	    std::reverse(s_.begin(),s_.end());
	}

	/** 
	 * \brief push back character
	 * 
	 * @param c character 
	 */
	void
	push_back(char c) {
	    s_.push_back(c);
	}
	
	/** 
	 * @brief Substring
	 * 
	 * @param pos start position of substring, 1-based
	 * @param len length of substring
	 * 
	 * @return substring at pos of length len 
	 */
	string1 
	substr(size_t pos, size_t len) const {
	    assert(pos>=1);
	    assert(pos<=s_.length());
	    return s_.substr(pos-1,len);
	}
	
    };

} // end namespace LocARNA

#endif // LOCARNA_STRING1_HH
