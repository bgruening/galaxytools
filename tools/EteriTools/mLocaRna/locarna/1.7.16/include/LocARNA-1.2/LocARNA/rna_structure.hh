#ifndef LOCARNA_RNA_STRUCTURE_HH
#define LOCARNA_RNA_STRUCTURE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <assert.h>
#include <set>
#include <string>

namespace LocARNA {
    /** 
     * @brief An RNA secondary structure
     *
     * Represents a structure (for a sequence of given length) as set
     * of base pairs. Supports parsing of dot-bracket strings
     * (potentially including pseudoknots) and traversal of base
     * pairs.
     *
     * Generally, base pairs (i,j) have to be oriented, i.e. i<j;
     * compare private method assert_valid_bp()
     */
    class RnaStructure {
    public:
	//! base pair type
	typedef std::pair<size_t,size_t> bp_t;

	//! base pair set type
	typedef std::set<bp_t,std::less<bp_t> > bps_t;
	
    private:
	size_t length_;
	bps_t bps_; 
	
	static const char unpaired_symbol_='.';
	static const std::string open_symbols_;
	static const std::string close_symbols_;
	
	bool
	parse(const std::string &s,
	      bps_t &bps,
	      char op,
	      char cl);

	bool
	parse(const std::string &s,
	      bps_t &bps, 
	      const std::string &op_syms,
	      const std::string &cl_syms);

	void 
	assert_valid_bp(const bp_t &bp) {
	    assert(1 <= bp.first);
	    assert(bp.first<bp.second);
	    assert(bp.second<=length_);
	}

    public:
	/** 
	 * @brief construct empty
	 */
	RnaStructure(): length_(0), bps_() {}

	/** 
	 * @brief construct from dot-bracket string
	 * 
	 * @param structure dot-bracket string
	 *
	 * We recognize different bracket pairs: (),[],{},<>, and
	 * letter pairs Aa, Bb, etc ; the structure string can encode
	 * crossing base pairs like in 
	 *
	 * .(((..[[...]]..)))..
	 *
	 * All such base pairs are recognized.
	 */
	RnaStructure(const std::string &structure);

	
	/** 
	 * @brief Equality operator
	 * 
	 * @param s rna structure to be compared
	 * @return whether equal
	 */
	bool
	operator ==(const RnaStructure &s) const {
	    return 
		this->length_ == s.length_
		&&
		this->bps_ == s.bps_;
	}

	/** 
	 * @brief Base pair for membership test
	 * 
	 * @param x base pair
	 * 
	 * @return whether structure contains the base pair
	 */
	bool
	contains(const bp_t &x) const {
	    return bps_.find(x) != bps_.end();
	}

	/** @brief sequence length 
	 */
	size_t
	length() const {return length_;}

	/** @brief number of base pairs
	 */
	size_t
	size() const {return bps_.size();}

	/** @brief insert base pair
	 * @param bp base pair
	 */
	void
	insert(const bp_t &bp) {
	    assert_valid_bp(bp);
	    bps_.insert(bp);
	}

	/** @brief remove base pair
	 * @param bp base pair
	 */
	void
	remove(const bp_t &bp) {
	    assert_valid_bp(bp);
	    bps_.erase(bp);
	}
	
	/** @brief clear structure
	    set structure to empty
	 */
	void
	clear() {bps_.clear();}

	
	// This is not well supported by the current structure
	// represenation. Therefore, we don't offer such functionality.
	// /**
	//  * @brief Check whether a loop contains a position
	//  * 
	//  * @param k position
	//  * @param x loop; (0,length+1) means external loop
	//  * 
	//  * @return whether the loop enclosed by x contains k
	//  *
	//  * Definition: the loop enclosed by (i,j) contains k iff i<k<j
	//  * and there is no base pair (i',j') in the structure, where
	//  * i<i'<k<j'<j. Note that we don't require (i,j) to be element
	//  * of the structure.
	//  *
	//  * k can be member of more than one loop unless nested().
	//  */
	// bp_t
	// in_loop_of(size_t k, bp_t x) const;
	

	//support contant iteration over base pair set
	
	//! constant iterator over base pairs
	typedef bps_t::const_iterator const_iterator;
	
	/** 
	 * @brief begin of base pair set 
	 * 
	 * @return constant iterator at begin  
	 */
	const_iterator begin() const {return bps_.begin();}
	
	/** 
	 * @brief end of base pair set 
	 * 
	 * @return constant iterator at end  
	 */
	const_iterator end() const {return bps_.end();}
	
	/** 
	 * @brief convert to dot-bracket string 
	 *
	 * @return dot-bracket string
	 *
	 * If the structure contains crossing base pairs, such base
	 * pairs are encoded using more than one pair of bracket
	 * symbols.  ( see constructor from dot bracket string ).  The
	 * use of bracket symbols is greedy from left to right,
	 * following the order defined in the class (by constants
	 * open_symbols_ and close_symbols_).
	 *
	 * @pre structure is in crossing class!
	 */
	std::string
	to_string() const;
	
    private:
	std::string
	to_string( const bps_t &bps ) const;

    public:
	
	/**
	 * @brief Check base pair set for empty structure / class PLAIN
	 *
	 * @param bps set of base pairs
	 *
	 * @return whether given base pair set represents empty structure
	 */
	static
	bool
	empty(const bps_t &bps);
	
	/**
	 * @brief Check for empty structure / class PLAIN
	 *
	 * @return whether structure is empty
	 */
	bool
	empty() const {return empty(bps_);}
	
	/**
	 * @brief Check for class NESTED
	 *
	 * @param bps set of base pairs
	 *
	 * @return whether structure is in the nested class (i.e., not
	 * in classes crossing or unlimited)
	 *
	 * A structure is in the class nested if no base pairs cross
	 * or share common ends.
	 *
	 * @see in_crossing() for naming
	 */
	static
	bool
	nested(const bps_t &bps);
	
	/**
	 * @brief Check for nested structure / class NESTED
	 *
	 * @return whether structure is nested class
	 */
	bool
	nested() const {return nested(bps_);}

	/**
	 * @brief Check for class CROSSING
	 *
	 * @param bps set of base pairs
	 *
	 * @return whether structure is in the crossing class (i.e.,
	 * not in unlimited)
	 *
	 * A structure is in the crossing class if no base pairs share
	 * common ends.
	 *
	 * @note this does *not* test for the presence of crossing
	 * base pairs. For the latter, use !nested().
	 */
	static
	bool
	crossing(const bps_t &bps);
	
	/**
	 * @brief Check for crossing structure / class CROSSING
	 *
	 * @return whether structure is in crossing class
	 */
	bool
	crossing() const {return crossing(bps_);}
	
    }; // end class RnaStructure


} // end namespace LocARNA

#endif

