#ifndef LOCARNA_SEQUENCE_HH
#define LOCARNA_SEQUENCE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <vector>
#include "multiple_alignment.hh"


namespace LocARNA {
    
    /**
     * @brief "Sequence View" of multiple alignment as array of column
     * vectors
     *
     * @note use MultipleAlignment::as_sequence() to 'convert' from
     * MultipleAlignment to Sequence.
     *
     * @note the conceptual relation between Sequence and
     * MultipleAlignent is that the two classes are equivalent views
     * of a multiple alignment. Therefore, casting between both should
     * be free of cost, in particular for references! operator [] is
     * supported, depending on the static type, only for
     * Sequence. There should be a way to expresse this relation in a
     * nicer way.
     */
    class Sequence:
	public MultipleAlignment {
	
	// BEWARE: don't define attributes in sequence! In
	// MultipleAlignment.as_sequence(), we rely on a simple upcast
	// to convert MultipleAlignment objects to Sequence objects
	
    public:
	
	/** 
	 * @brief Construct empty
	 */
 	Sequence(): MultipleAlignment() {}
	
	/** 
	 * @brief Construct as single sequence
	 * @param name name of sequence
	 * @param sequence sequence string
	 */
	Sequence(const std::string &name,
		 const std::string &sequence)
	    : MultipleAlignment(name,sequence) {
	}
	
	/** 
	 * @brief Access to columns
	 * 
	 * @param col_index column index
	 * 
	 * @return alignment column (proxy class)
	 *
	 * @note allows array notation via [] operator; this is the
	 * main difference to MultipleAlignment class
	 */
	AliColumn
	operator [](size_type col_index) const {
	    return column(col_index);
	}
	
	/** 
	 * \brief names vector (legacy, deprecated)
	 * 
	 * @return vector of sequence names
	 * @note deprecated: in place of names()[i], rather use seqentry(i).name()
	 */
	std::vector<std::string> 
	names() const;
	
    };
    
} // end namespace LocARNA

#endif
