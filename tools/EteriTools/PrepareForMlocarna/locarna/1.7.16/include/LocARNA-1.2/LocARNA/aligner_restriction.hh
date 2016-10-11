#ifndef LOCARNA_ALIGNER_RESTRICTION
#define LOCARNA_ALIGNER_RESTRICTION

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iostream>

namespace LocARNA {

    /**
       @brief Restricts range of an alignment in Aligner
       
       Contains information for restricting Aligner to sub-sequences
       startA..endA amd startB..endB.

       Take care when using aligner restrictions for multiple
       Alignments with the same aligner object.
       The D-matrix is only computed once, so this works
       as long as the first Aligner::align() is called with
       the most general restriction (e.g. no restriction at all)!
       
       @see Aligner
    */
    class AlignerRestriction {
    private:
	int startA; //!< start position in A
	int startB; //!< start position in B
	int endA; //!< end position in A
	int endB; //! end position in B
    public:
	
	/** 
	 * Constructs with start and end positions of subsequences
	 * 
	 * @param startA_ start position in A
	 * @param startB_ start position in B
	 * @param endA_ end position in A
	 * @param endB_ end position in B
	 */
	AlignerRestriction(int startA_, int startB_, int endA_, int endB_)
	    : startA(startA_), startB(startB_), endA(endA_), endB(endB_)
	{}
	
	/** 
	 * Read access to member
	 * 
	 * @return start position in A
	 */
	size_t get_startA() const {return startA;}
	
	/** 
	 * Read access to member
	 * 
	 * @return end position in A
	 */
	size_t get_endA() const {return endA;}
	
	/** 
	 * Read access to member
	 * 
	 * @return start position in B
	 */
	size_t get_startB() const {return startB;}

	/** 
	 * Read access to member
	 * 
	 * @return end position in B
	 */
	size_t get_endB() const {return endB;}
	
	/** 
	 * Write access to member
	 * 
	 * @param p start position in A
	 */
	void set_startA(size_t p) {startA=p;}

	/** 
	 * Write access to member
	 * 
	 * @param p end position in A
	 */
	void set_endA(size_t p)  {endA=p;}
	
		/** 
	 * Write access to member
	 * 
	 * @param p start position in B
	 */
	void set_startB(size_t p)  {startB=p;}
	
	/** 
	 * Write access to member
	 * 
	 * @param p end position in B
	 */
	void set_endB(size_t p)  {endB=p;}
    };

    /** 
     * Output operator for objects of AlignerRestrictions
     * 
     * @param out output stream
     * @param r   object of AlignerRestriction to be written to stream
     * 
     * @return output stream
     * @note Writes r to out
     */
    inline
    std::ostream & operator<<(std::ostream &out, AlignerRestriction r) {
	return
	    out << r.get_startA() << " "
		<< r.get_startB() << " "
		<< r.get_endA() << " "
		<< r.get_endB();
    }

}

#endif
