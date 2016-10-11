#ifndef LOCARNA_ALIGNER_IMPL_HH
#define LOCARNA_ALIGNER_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "aligner.hh"

#include "aligner_restriction.hh"
#include "scoring.hh"
#include "alignment.hh"
#include "arc_matches.hh"

namespace LocARNA {

    class Sequence;
    template <class T> class Matrix;

    /**
     * @brief Implementation of Aligner
     */
    struct AlignerImpl {

	/**
	 * type of matrix M
	 * @note 'typedef RMtrix<infty_score_t> M_matrix_t;' didn't improve performance
	 */
	typedef ScoreMatrix M_matrix_t;

	//! an arc
	typedef BasePairs__Arc Arc;
	
	
	const Scoring *scoring_; //!< the scores
	Scoring *mod_scoring_; //!< used in normalized scoring, when we need to modify the scoring
    
	const AlignerParams *params_; //!< the parameter for the alignment
	
	const Sequence &seqA_; //!< sequence A
	const Sequence &seqB_; //!< sequence B
    
	const ArcMatches &arc_matches_; //!< the potential arc matches between A and B
    
	const BasePairs &bpsA_; //!< base pairs of A
	const BasePairs &bpsB_; //!< base pairs of B
    
	/**
	   \brief restriction of alignment for k-best
       
	   The AlignerRestriction r is used for the k-best local
	   alignments.
       
	   Aligner always works on sub-sequences/structures of seqA and
	   seqB. These are specified by r.
       
	   In the standard case, i.e. align the whole sequences, the
	   values are set to startA=1 and endA=seqA.length() and
	   analogously for B
       
	   The constructor initializes to the above standard values
	   methods provide possibility to restriction.
	*/
	AlignerRestriction r_;

	//! matrix indexed by the arc indices of rnas A and B
	ScoreMatrix Dmat_;
    
	/**
	 * M matrices
	 * @note in the case of structure local alignment, 
	 * the algo uses eight M matrices 
	 */
	std::vector<M_matrix_t> Ms_;
    
	/**
	 * for cool affine gap cost, we need two additional matrices E
	 * and F.  However we only need to store one row for E and one
	 * scalar for F.  
	 * @note for structure local, we need one such
	 * matrix per state 0..3
	 * @see Fs
	 */
	std::vector<ScoreVector> Es_;
    
	/**
	 * for affine gap cost.
	 * @see Es
	*/
	std::vector<infty_score_t> Fs_;
    
	int min_i_; //!< subsequence of A left end, computed by trace back
	int min_j_; //!< subsequence of B left end, computed by trace back

	int max_i_; //!< subsequence of A right end, computed by align_top_level
	int max_j_; //!< subsequence of B right end, computed by align_top_level
    
	bool D_created_; //!< flag, is D already created?
    
	Alignment alignment_; //!< resulting alignment
    
	/**
	 * \brief different states for computation of structure-local alignment.
	 *
	 * \note The idea of the names is E=exclusion, NO=no exclusion, X=one exclusion,
	 * OP=open exclusion. In E_1_2, 1 refers to sequence A and 2 to sequence B.
	 */
	enum {E_NO_NO, E_X_NO, E_NO_X, E_X_X,
	      E_OP_NO, E_NO_OP, E_OP_X, E_X_OP};


	// ============================================================
	/**
	 * @brief Provides the standard view on the scoring
	 * @see ModifiedScoringView
	 *
	 * @note We use a template-based scheme to switch between use of
	 * the unmodified score and the modified score without run-time
	 * penalty for the standard case the mechanism is used for methods
	 * align_noex and trace_noex
	 */
	class UnmodifiedScoringView {
	private:
	    const AlignerImpl *aligner_impl_; //!< aligner object for that the view is provided
	public:
	
	    /** 
	     * Construct for Aligner object
	     * 
	     * @param aligner_impl the aligner implementation object
	     */
	    UnmodifiedScoringView(const AlignerImpl *aligner_impl): aligner_impl_(aligner_impl) {};
	    
	    /** 
	     * Get scoring object
	     * 
	     * @return (unmodified) scoring object of aligner
	     */
	    const Scoring *scoring() const {return aligner_impl_->scoring_;}
	
	    /** 
	     * View on matrix D
	     * 
	     * @param a arc in A
	     * @param b arc in B
	     * 
	     * @return D matrix entry for match of a and b
	     */
	    infty_score_t D(const Arc &a, const Arc &b) const {
		return aligner_impl_->Dmat_(a.idx(),b.idx());
	    }
	
	    /** 
	     * View on matrix D
	     * 
	     * @param am arc match
	     * 
	     * @return D matrix entry for arc match am
	     */
	    infty_score_t D(const ArcMatch &am) const {
		return D(am.arcA(),am.arcB());
	    }
	};
    
    
	/**
	 * @brief Provides a modified view on the scoring
	 *
	 * This view is used when
	 * computing length normalized local alignment.  
	 * @see UnmodifiedScoringView
	 */
	class ModifiedScoringView {
	private:
	    const AlignerImpl *aligner_impl_; //!< aligner object for that the view is provided
    	
	    score_t lambda_; //!< factor for modifying scoring
	
	    /** 
	     * Computes length of an arc
	     * 
	     * @param a the arc
	     * 
	     * @return length of arc a
	     */
	    size_t
	    arc_length(const Arc &a) const {
		return a.right()-a.left()+1;
	    }
	public:

	    /** 
	     * Construct for Aligner object
	     * 
	     * @param aligner_impl The aligner implementation object
	     *
	     * @note scoring object in aligner has to be modified by lambda already
	     */
	    ModifiedScoringView(const AlignerImpl *aligner_impl)
		: aligner_impl_(aligner_impl),lambda_(0) {}
	
	    /** 
	     * Change modification factor lambda
	     * 
	     * @param lambda modification factor
	     */
	    void
	    set_lambda(score_t lambda) {
		lambda_=lambda;
	    }
	
	    /** 
	     * Get scoring object
	     * 
	     * @return modified scoring object of aligner
	     */
	    const Scoring *scoring() const {return aligner_impl_->mod_scoring_;}

	    /** 
	     * View on matrix D
	     * 
	     * @param a arc in A
	     * @param b arc in B
	     * 
	     * @return modified D matrix entry for match of a and b
	     */
	    infty_score_t D(const Arc &a,const Arc &b) const {
		return aligner_impl_->Dmat_(a.idx(),b.idx())
		    -FiniteInt(lambda_*(arc_length(a)+arc_length(b)));
	    }

	    /** 
	     * View on matrix D
	     * 
	     * @param am arc match
	     * 
	     * @return modified D matrix entry for arc match am
	     */
	    infty_score_t D(const ArcMatch &am) const {
		return aligner_impl_->Dmat_(am.arcA().idx(),am.arcB().idx())
		    -FiniteInt(lambda_*(arc_length(am.arcA())+arc_length(am.arcB())));
	    }
	};
    

	const UnmodifiedScoringView def_scoring_view_; //!< Default scoring view
	ModifiedScoringView mod_scoring_view_; //!< Modified scoring view for normalized alignment
    
	// ============================================================
	
	/** 
	 * @brief copy constructor
	 * 
	 * @param a Aligner implementation
	 */
	AlignerImpl(const AlignerImpl &a);
	
	/** 
	 * @brief Construct from parameters
	 * 
	 * @param seqA sequence A
	 * @param seqB sequence B
	 * @param arc_matches arc matches
	 * @param ap parameter for aligner
	 * @param s scoring object
	 */
	AlignerImpl(const Sequence &seqA, 
		    const Sequence &seqB,
		    const ArcMatches &arc_matches,
		    const AlignerParams *ap,
		    const Scoring *s
		    );

	/** 
	 * Destructor
	 */
	~AlignerImpl();
	
	// ============================================================
	

	/**
	 * \brief initialize matrices M and E
	 *
	 * initialize first column and row of matrices M and E for
	 * the alignment below of arc match (a,b).
	 * The initialization depends on the state.
	 * First row/column means the row al and column bl.
	 * For correct initialization (in particular in local modes),
	 * globalA/B and exclA/B need to be given correctly
	 * for the state!
	 *
	 * @param state the state, selects the matrices M,E
	 * @param al left end of arc a
	 * @param ar right end of arc a
	 * @param bl left end of arc b
	 * @param br right end of arc b
	 * @param globalA allow no free deletion of prefix of sequence A
	 * @param exclA allow deletion of prefix of sequence A with cost exclusion()
	 * @param globalB analogous for sequence B
	 * @param exclB analogous for sequence B
	 * @param sv Scoring view
	 * 
	*/
	template <class ScoringView>
	void init_state(int state, pos_type al, pos_type ar, 
			pos_type bl, pos_type br, 
			bool globalA, bool exclA,
			bool globalB, bool exclB, 
			ScoringView sv);
    
    
	/**
	 * \brief standard cases for alignment (without exlusion handling).
	 *
	 * recursion cases that handle everything but exclusions
	 * (in the LSSA-paper this function was called NoEx
	 *    
	 * @param state necessary for structure local, there state refers to a set of matrices M,E,F
	 * @param al position in sequence A: left end of current arc match
	 * @param bl position in sequence B: left end of current arc match
	 * @param i position in sequence A, for which score is computed
	 * @param j position in sequence B, for which score is computed
	 * @param sv the scoring view to be used
	 * @returns score of i,j in matrix set state that results from standard cases
	 * 
	 * @pre state in 0..4, in non-structure local alignment state has to be 0;
	 * @pre i,j is allowed by edge controller
	 */
	template<class ScoringView>
	infty_score_t align_noex(int state, pos_type al, pos_type bl, pos_type i, pos_type j, ScoringView sv);
     
	/**
	 * align the loops closed by arcs (al,ar) and (bl,br).
	 * in structure local alignment, this allows to introduce exclusions
	 *
	 * @param al left end of arc a
	 * @param ar right end of arc a
	 * @param bl left end of arc b
	 * @param br right end of arc b
	 * @param allow_exclusion whether to allow exclusions
	 * 
	 * @pre arc-match (al,ar)~(bl,br) valid due to constraints and heuristics
	 */
	void align_in_arcmatch(pos_type al,pos_type ar,pos_type bl,pos_type br,
			       bool allow_exclusion);
  

	/**
	 * align the top-level with potential free end gaps
	 * and return the maximal score
	 */
	infty_score_t
	align_top_level_free_endgaps();

    
	/**
	 * align the top-level in a sequence local alignment
	 * and return the maximal score
	 */
	template<class ScoringView>
	infty_score_t align_top_level_locally(ScoringView sv);
    
	//! align top level in the scanning version
	infty_score_t align_top_level_localB();
  
	/** 
	 * \brief trace back within an match of arcs
	 * 
	 * @param state the state selects M/E/F matrices (used in structure local alig)
	 * @param al left end of arc in A
	 * @param i  right end of subsequence in A
	 * @param bl left end of arc in B
	 * @param j right end of subsequence in B
	 * @param top_level whether alignment is on top level
	 * @param sv scoring view 
	 */
	template<class ScoringView>
	void trace_in_arcmatch(int state,int al,int i,int bl,int j,bool top_level,ScoringView sv);
    
	/** 
	 * \brief standard cases in trace back (without handling of exclusions)
	 * 
	 * @param state the state selects M/E/F matrices (used in structure local alig)
	 * @param al left end of arc in A
	 * @param i  right end of subsequence in A
	 * @param bl left end of arc in B
	 * @param j right end of subsequence in B
	 * @param top_level whether alignment is on top level
	 * @param sv scoring view 
	 */
	template <class ScoringView>
	void trace_noex(int state,
			pos_type al, pos_type i,
			pos_type bl,pos_type j,
			bool top_level,
			ScoringView sv);
    
	/**
	 * trace an arc match
	 * @param am the arc match
	 */
	void trace_arcmatch(const ArcMatch &am);

	/**
	 * trace an arc match in case of forbidden lonely pairs
	 * @param am the arc match
	 */
	void trace_arcmatch_noLP(const ArcMatch &am);

	//! compute the alignment score
	infty_score_t
	align();
	
	/**
	   create the entries in the D matrix
	   This function is called by align() (unless D_created)
	*/
	void align_D();

	/**
	   fill in D the entries with left ends al,bl
	*/
	void 
	fill_D_entries(pos_type al, pos_type bl);
    
	/**
	   fill D entries when no-lonely-pairs option given
	   after computation of M-matrices with left ends al,bl
       
	   this fills D entries for arc matches with left ends al-1,bl-1,
	   since the positions refer to the inner arc
	   of a stacked arc pair
	*/
	void 
	fill_D_entries_noLP(pos_type al, pos_type bl);
    
	/** 
	 * Read/Write access to D matrix
	 * 
	 * @param am Arc match
	 * 
	 * @return entry of D matrix for am
	 */
	infty_score_t &D(const ArcMatch &am) {
	    return Dmat_(am.arcA().idx(),am.arcB().idx());
	}

	/** 
	 * Read/Write access to D matrix
	 * 
	 * @param arcA arc in sequence A
	 * @param arcB arc in sequence B
	 * 
	 * @return entry of D matrix for match of arcA and arcB
	 */
	infty_score_t &D(const Arc &arcA,const Arc &arcB) {
	    return Dmat_(arcA.idx(),arcB.idx());
	}

	/**
	 * do the trace back through the alignment matrix
	 * with partial recomputation
	 * pre: call align() to fill the top-level matrix
	 */
	template <class ScoringView>
	void trace(ScoringView sv);


    };
	
} // end namespace LocARNA

#endif // LOCARNA_ALIGNER_IMPL_HH
