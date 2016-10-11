#ifndef LOCARNA_ALIGNER_HH
#define LOCARNA_ALIGNER_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "aux.hh"
#include "scoring_fwd.hh"
#include "rna_structure.hh"


namespace LocARNA {

    class AlignerImpl;

    class Sequence;
    class AlignerParams;
    class ArcMatches;
    class Alignment;
    
    class AlignerRestriction;
    
    /**
     * \brief Implements comparison by member second
     *
     * Templated function class implementing a comparison operator for
     * member second of class T.
     * @note used for priority queue in Aligner::suboptimal
     */
    template <class T>
    class greater_second {
    public:
	/** 
	 * Compare members second of class T
	 * 
	 * @return whether a.second smaller b.second
	 */
	bool operator() (const T& a, const T&b) {
	    return a.second < b.second;
	}
    };
    
    
    /**
     * \brief Implements locarna alignment algorithm
    
     * Performs the alignment of two sequences and their associated
     * sets of weighted basepairs
    
     * An object always knows about the two sequences and the two
     * weighted base pair sets
    
     * usage: construct, align, trace, get_alignment
    
     * @note Idea "NICE TO HAVE": score matrices should be as small as
     * possible and have offsets.  The M-Matrix can be smaller if the
     * maximal arc len is limited. The D-matrix may be smaller if the
     * number of simultaneously needed arc-pairs is limited, due to
     * limit on the local sub-sequence lengths. A special cyclically
     * rotatable matrix would be required for alignment on the
     * top-level. Such matrices are implemented in matrix.hh but
     * currently not used.
     */
    class Aligner {
	AlignerImpl *pimpl_;

    public:
    

	/** 
	 * @brief copy constructor
	 * @param aligner object to be copied
	 * Copies implementation object (not only pointer) 
	 */
	Aligner(const Aligner &aligner);
	
	/** 
	 * @brief assignment operator
	 * @param aligner object to be assigned
	 * Assigns implementation object (not only pointer) 
	 */
	Aligner &operator =(const Aligner &aligner);
		
	/** 
	 * @brief Construct from parameters
	 * 
	 * @param seqA sequence A
	 * @param seqB sequence B
	 * @param arc_matches arc matches
	 * @param ap parameter for aligner
	 * @param s scoring object
	 */
	Aligner(const Sequence &seqA,
		const Sequence &seqB,
		const ArcMatches &arc_matches,
		const AlignerParams *ap,
		const Scoring *s
		);
    
	//! destructor
	~Aligner();

	//! return the alignment that was computed by trace()
	Alignment const & 
	get_alignment() const; 

	/** 
	 * @brief set the alignment
	 * @param alignment the alignment to be stored in the Aligner object
	 *
	 * This sets the result of the aligner object / trace; it is
	 * useful to evaluate given alignments (where it is used in
	 * place of a computation by align()/trace()).
	 */
	void
	set_alignment(const Alignment &alignment);
	
	//! compute the alignment score
	infty_score_t
	align();
    
	//! offer trace as public method. Calls trace(def_scoring_view).
	void
	trace();
    
	/**
	 * set the restriction on the alignment,
	 * mainly used for the k-best algorithm
	 */
	void 
	set_restriction(const AlignerRestriction &r);
	
	/**
	 * return the current restriction,
	 * mainly used for the k-best algorithm
	 */
	const AlignerRestriction &
	get_restriction() const;
	
  
	/**
	 * @brief Enumerate suboptimal local alignments
	 * 
	 * special operation mode for computing the k best alignments.
	 * Used in place of calls to align() and trace()
	 *
	 * @param k 
	 * @param threshold 
	 * @param opt_normalized 
	 * @param normalized_L 
	 * @param output_width 
	 * @param opt_verbose 
	 * @param opt_local_output 
	 * @param opt_pos_output 
	 * @param opt_write_structure 
	 */
	void 
	suboptimal(size_t k,
		   score_t threshold,
		   bool opt_normalized,
		   score_t normalized_L,
		   size_t output_width,
		   bool opt_verbose,
		   bool opt_local_out,
		   bool opt_pos_output,
		   bool opt_write_structure
		   );
    
    
	//! perform normalized local alignment with parameter L
	infty_score_t
	normalized_align(score_t L, bool opt_verbose);

	/** 
	 * \brief evaluate the alignment according to scoring and scoring parameters
	 * 
	 * @return score of the alignment 
	 */
	score_t 
	evaluate();

	/** 
	 * \brief find optimum consensus structure
	 * 
	 * @return optimum consensus structure for alignment
	 */
	RnaStructure
	optimize_consensus_structure();
	

    };

} //end namespace LocARNA

#endif // LOCARNA_ALIGNER_HH
