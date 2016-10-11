#ifndef LOCARNA_TRACE_CONTROLLER_HH
#define LOCARNA_TRACE_CONTROLLER_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <vector>
#include <assert.h>

#include "aux.hh"
#include "multiple_alignment.hh"

namespace LocARNA {

    class Sequence;
    
    /**
     * \brief Represents a range of traces
     *
     * Represents a range of possible traces in a dynamic
     * programming matrix for aligning two sequences
     * @see TraceController
     */
    class TraceRange {
    public:
	//! alias for MultipleAlignment::SeqEntry
	typedef MultipleAlignment::SeqEntry SeqEntry;
	
	//! pair of sequence entries
	typedef std::pair<SeqEntry, SeqEntry> seqentry_pair_t;
	
	/** 
	 * Remove common gaps in pairwise alignment
	 * 
	 * @param aliA alignment string A with name (SeqEntry)
	 * @param aliB alignment string B with name (SeqEntry)
	 * 
	 * @return pair of "sequence entries" that represent the
	 * pairwise alignment given by aliA and aliB where gap-only
	 * columns are removed
	 */
	static
	seqentry_pair_t
	remove_common_gaps(const SeqEntry &aliA,
			   const SeqEntry &aliB);

    protected:		
	std::vector<size_t> min_col_vector; //!< minimal column in row
	std::vector<size_t> max_col_vector; //!< maximal column in row
    
    public:
	/**
	 * \brief Construct from pair of alignment strings
	 *
	 * construct trace range of two sequences given two alignment strings
	 * of the sequences and the allowed deviation delta
	 * the sequences can contain gaps themselves (which happens,
	 * when sequences orignate from a sequence profile).
	 *
	 * @param pseqA SeqEntry of sequence A
	 * @param pseqB SeqEntry of sequence B
	 * @param aliA alignment SeqEntry for sequence A
	 * @param aliB alignment SeqEntry for sequence B
	 * @param delta the allowed deviation
	 *
	 * side conditions:
	 * remove_gaps(seqA) == remove_gaps(aliA)
	 * && remove_gaps(seqB) == remove_gaps(aliB)
	 * where remove_gaps is a function that removes all gap symbols
	 * length(aliA)==length(aliB)
	 *
	 */
	TraceRange(const SeqEntry &pseqA,
		   const SeqEntry &pseqB,
		   const SeqEntry &aliA,
		   const SeqEntry &aliB,
		   size_type delta);

	
	/** 
	 * @brief Construct as consensus trace range from a set of traces
	 * 
	 * Constructs object as consensus trace ranges of the traces
	 * trs.  The construction follows an idea of relaxing the
	 * deviation constraints by computing a trace with minimal
	 * accumulated distance to all traces and then determining its
	 * delta environment.
	 *
	 * @param lenA length of sequence A
	 * @param lenB length of sequence B
	 * @param trs set of traces
	 * @param delta deviation
	 */
	TraceRange(size_type lenA, size_type lenB, const std::vector<TraceRange> &trs, size_type delta);

	//! \brief Construct empty
	TraceRange() {
	}
    
	/**
	 * \brief Computes cost of a cut in the consensus trace of a trace range set
	 * @param i cut.first  
	 * @param j cut.second
	 * @param trs set of trace ranges
	 * @return cost of cut (i,j) in consensus of trs
	 */
	size_type
	consensus_cost(size_type i,
		       size_type j,
		       const std::vector<TraceRange> &trs) const;
	
	/**
	 * \brief Read number of rows
	 * @return length of seqA, i.e. the maximal row of the trace
	 */
	size_t
	rows() const {return min_col_vector.size()-1;}
	
	/**
	 * \brief Minimal column of trace in a row
	 * @param i: row of matrix, 0<=i<=rows()
	 * @returns minimal valid trace cell in the row i
	 */
	size_t
	min_col(size_t i) const {return min_col_vector[i];}
	
	/**
	 * \brief Maximal column of trace in a row
	 * @param i: row of matrix, 0<=i<=rows()
	 * @returns maximal valid trace cell in the row i
	 */
	size_t 
	max_col(size_t i) const {return max_col_vector[i];}

	/** 
	 * Print object to ouptut stream for debugging
	 * 
	 * @param out output stream
	 */
	void
	print_debug(std::ostream & out) const;
	
    };


    //! abstract class that declares the method is_valid_match()
    class MatchController {
    
    public:
	/**
	 * test for allowed matches due to valid traces
	 * @param i position in sequence A in 1..lenA
	 * @param j position in sequence B in 1..lenB
	 * @returns whether i~j is an allowed match due to valid traces
	 */
	virtual
	bool
	is_valid_match(size_t i, size_t j) const=0;

	virtual 
	~MatchController();
    
    };


    /**
     * @brief Controls the matrix cells valid for traces
     *
     * Controls the matrix cells that need to be filled
     * in a dynamic programming algorithm because they occur on valid
     * traces due to the max-diff heuristic
     *
     * The valid traces can be defined either unrestricted, or by a
     * maximal difference of i and j for matrix cells (i,j) or due to
     * a maximal difference to a given alignment (trace). 
     */
    class TraceController : public TraceRange, public MatchController {

    private:
	// The delimiter character separating the two sequences in the alignment string
	static const char delimiter = '&';
	    
	TraceRange trace_range;
    
	/**
	 * merge in the given trace with delta into current trace range
	 * @param trace the new trace
	 */
	void
	merge_in_trace_range(const TraceRange &tr);

	//! The allowed distance in computing the min and max positions.
	const size_type delta;

	/**
	 * switch between strict and relaxed merging of pairwise trace
	 * ranges
	 */
	const bool relaxed_merging; 
    
    public:
    
	/**
	 * \brief Constructs for the general case of alignment of alignments
	 * @param seqA sequence A
	 * @param seqB sequence B
	 * @param ma multiple reference alignment
	 * @param delta the allowed difference
	 * @param relaxed_merging whether to use relaxed merging of trace ranges 
	 * @note If delta == -1 then min_col is 1 and max_col is lenB
	 * @note If delta != -1 and ma==NULL, then define min j, max j by deviation |i-(lenA/lenB)*j|<=delta
	 * @note These values are chosen such that for all j between
	 *  min and max, the delta constraint holds for all pairs of
	 *  sequences in seqA and seqB.
	 */
	TraceController(const Sequence &seqA, 
			const Sequence &seqB, 
			const MultipleAlignment *ma, 
			int delta,
			bool relaxed_merging=false);

	//! \brief Virtual destructor
	virtual
	~TraceController();
        
	/**
	 * test for matrix entries on valid trace
	 * @param i position in sequence A in 1..lenA or 0
	 * @param j position in sequence B in 1..lenB or 0
	 * @returns whether matrix cell (i.j) is valid 
	 */
	bool
	is_valid(size_type i, size_type j) const;

	/**
	 * test for allowed matches due to valid traces
	 * @param i position in sequence A in 1..lenA
	 * @param j position in sequence B in 1..lenB
	 * @returns whether i~j is an allowed match due to valid traces
	 */
	virtual
	bool
	is_valid_match(size_type i, size_type j) const;

	/**
	 * \brief Read deviation
	 * @return deviation Delta
	 */
	size_type get_delta() const {return delta;}

	/**
	 * print debugging information to stream
	 * @param out output stream
	 */
	void
	print_debug(std::ostream & out) const;

    private:
	/**
	 * constrain the min/max j without reference alignment by
	 * delta only such that match i~j is allowed iff | i/n - j/m |
	 * <= delta/((n+m)/2), unless delta is too small to connect
	 * the matrix entries; in the latter case delta is chosen
	 * minimally.
	 */
	void
	constrain_wo_ref(size_type lenA, size_type lenB, size_type delta);

    };

    inline
    bool
    TraceController::is_valid(size_type i, size_type j) const {
	return min_col(i)<=j && j<=max_col(i);
    }

    inline
    bool
    TraceController::is_valid_match(size_type i, size_type j) const {
	return is_valid(i,j) && is_valid(i-1,j-1);
    }

} //end namespace


#endif /* LOCARNA_TRACE_CONTROLLER_HH */
