#ifndef LOCARNA_ARC_MATCHES_HH
#define LOCARNA_ARC_MATCHES_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <algorithm>
#include <vector>
#include <tr1/unordered_map>

#include "scoring_fwd.hh"
#include "aux.hh"
#include "matrix.hh"
#include "basepairs.hh"

#include <assert.h>


namespace LocARNA {
    class Scoring;
    class Sequence;
    class RnaData;
    class AnchorConstraints;
    class TraceController;
    class MatchController;
    
    /**
     * @brief Represents a match of two base pairs (arc match) 
     *
     * Maintains pointers to the single arcs (in the base pairs
     * strucures) and an arc match index.
     *
     * @see ArcMatches
     */
    class ArcMatch {
    public:
	typedef std::vector<int>::size_type size_type; //!< size type
	typedef size_type idx_type; //!< arc match index
	typedef BasePairs__Arc Arc; //!< arc 
    private:
	const Arc *arcA_; //!< the arc in A
	const Arc *arcB_; //!< the arc in B
	idx_type idx_;    //!< (unique) index of arc match

    public:
	
	/** 
	 * Construct with member values
	 * 
	 * @param arcA arc in B
	 * @param arcB arc in A
	 * @param idx  index
	 */
	ArcMatch(const Arc *arcA,const Arc *arcB, idx_type idx)
	    : arcA_(arcA),
	      arcB_(arcB),
	      idx_(idx)
	{}
    
	/** 
	 * Read access
	 * 
	 * @return the matched arc in A
	 */
	const Arc &
	arcA() const {return *arcA_;}
    
	/** 
	 * Read access
	 * 
	 * @return the matched arc in B
	 */
	const Arc &
	arcB() const {return *arcB_;}
    
	/** 
	 * Read access
	 * 
	 * @return (unique) index of arc match
	 */
	idx_type
	idx()  const {return idx_;}
    };

    //! Vector of arc matches
    typedef std::vector<ArcMatch> ArcMatchVec;

    //! Vector of arc match indices
    typedef std::vector<ArcMatch::idx_type> ArcMatchIdxVec;

    /**
       @brief Maintains the relevant arc matches and their scores
   
       It works as an interface between the source of arc match scoring
       information, i.e.  dot plots of the single sequences or explicit
       listing of all arc matches, and the use of this information in the
       class Scoring and in the alignment algorithm. For the latter
       purpose, it allows to generate the necessary objects of class
       BasePairs for the iteration of base pairs in the single structures.

       ArcMatches knows about all arc matches that shall be considered for the alignment.
       Each arc match gets an index between 0..(number of arc_matches-1).
       The index is used to access the arc match score, the single arcs of the match,
       the matrix D, ...
   
       The object offers iteration over 1.) all arc matches 2.) all arc
       matches that share left ends/right ends (i,j) During iteration the
       index of the current arc match is known, thus one never needs to
       infer the arc match index from the arc ends or arc indices!
       (This could be done in constant time using an O(n^2) lookup table.)
    */
    class ArcMatches {
    public:
	typedef std::vector<int>::size_type size_type; //!< size
	typedef BasePairs__Arc Arc; //!< arc
    protected:
	
	size_type lenA; //!< length of sequence A
	size_type lenB; //!< length of sequence B
        
	BasePairs *bpsA; //!< base pairs of RNA A
	BasePairs *bpsB; //!< base pairs of RNA B
    
	/* Constraints and Heuristics */
    
	size_type max_length_diff; //!< for max-diff-am heuristics
    
	const MatchController &match_controller; //!< allowed alignment traces by max-diff heuristics
    
	const AnchorConstraints &constraints; //!< for constraints
    
    
	/**
	 * decide according to constraints
	 * and heuristics whether an arc match is valid.
	 * @param arcA arc (i,j) in first sequence
	 * @param arcB arc (k.l) in second sequence
	 * @returns whether match of arcA and arcB is valid
	 * An arc match is valid, if and only if:
	 * 1.) matches i~k and j~l are valid due to trace_controller (max-diff-match heuristic) and constraints (anchor constraints)
	 * 2.) length difference of arcs <= max_length_diff
	 */
	bool is_valid_arcmatch(const Arc &arcA,const Arc &arcB) const;
    
	/* END constraints and heuristics */
    
	bool maintain_explicit_scores; //!< whether scores are maintained explicitely or computed from pair probabilities
    
	//! vector of all maintained arc matches
	ArcMatchVec arc_matches_vec;

	/**
	 * the number of valid arc matches
	 * @note this number can differ from the size of arc_matches_vec 
	 */
	size_type number_of_arcmatches;

	//! vector of scores (of arc matches with the same index)
	std::vector<score_t> scores;
    

	//! for each (i,j) maintain vector of the indices of the arc matchs that share the common right end (i,j) 
	Matrix<ArcMatchIdxVec> common_right_end_lists;
       
    
	//! for each (i,j) maintain vector of the indices of the arc matchs that share the common left end (i,j) 
	Matrix<ArcMatchIdxVec> common_left_end_lists;
    
    
	//! vector of indices of inner arc matches 
	ArcMatchIdxVec inner_arcmatch_idxs;
    
	//! initialize the vector of inner arc match indices
	void
	init_inner_arc_matchs();
    
	//! Compare two arc match indices by lexicographically comparing their left ends
	class lex_greater_left_ends {
	    const ArcMatches &arc_matches;
	public:
	    
	    /** 
	     * Construct with access to arc matches
	     * 
	     * @param arc_matches_ 
	     */
	    lex_greater_left_ends(const ArcMatches &arc_matches_)
		: arc_matches(arc_matches_)
	    {}
	    
	    /** 
	     * @brief Compare to arc matches
	     *
	     * Compares two arc matches by their left ends lexicographically 
	     *
	     * @param i index of first arc match
	     * @param j index of second arc match
	     *
	     * @return whether the first arc match is larger than the second
	     */
	    bool
	    operator () (const ArcMatch::idx_type &i, const ArcMatch::idx_type &j) const {
		size_type ali = arc_matches.arcmatch(i).arcA().left();
		size_type bli = arc_matches.arcmatch(i).arcB().left();
		size_type alj = arc_matches.arcmatch(j).arcA().left();
		size_type blj = arc_matches.arcmatch(j).arcB().left();
	    
		return (ali>alj) || (ali==alj && bli>blj);
	    }
	};

	/**
	 * A simple 5-tuple of 4 positions and a score
	 * 
	 */
	class tuple5 {
	public:
	    typedef std::vector<int>::size_type size_type; //!< size type

	    size_type i; //!< position i
	    size_type j; //!< position j
	    size_type k; //!< position k
	    size_type l; //!< position l
	    score_t score; //!< the score (as used below: score of arc match (i,j)~(k,l))
	
	    /** 
	     * Construct with member values
	     * 
	     * @param i_ position i
	     * @param j_ position j
	     * @param k_ position k
	     * @param l_ position l
	     * @param score_ the score
	     */
	    tuple5(size_type i_,size_type j_,size_type k_,size_type l_,score_t score_)
		:i(i_), j(j_), k(k_), l(l_), score(score_)
	    {}
	};
	    
    
    public:

	
	/** 
	 *  \brief construct with explicit arc match score list
	 * 
	 * construct from seqnames and explicit list of all scored arc
	 * matches together with their score.
	 *
	 *
	 * @note registers constraints and heuristics and then calls read_arcmatch_scores.
	 * The constructed object explicitely represents/maintains the scores of arc matchs. 
	 *
	 * @param seqA_ sequence A
	 * @param seqB_ sequence B
	 * @param arcmatch_scores_file file containing arc match scores
	 * @param probability_scale if >=0 read probabilities and multiply them by probability_scale 
	 * @param max_length_diff accept arc matches only up to maximal length difference
	 * @param trace_controller accept only due to trace controller
	 * @param constraints accept only due to constraints 
	 */
	ArcMatches(const Sequence &seqA_, 
		   const Sequence &seqB_, 
		   const std::string &arcmatch_scores_file,
		   int probability_scale,
		   size_type max_length_diff,
		   const MatchController &trace_controller,
		   const AnchorConstraints &constraints);
	
	
	
	/** 
	 *  \brief construct from single base pair probabilities. 
	 * 
	 * In this case, the object filters for relevant base
	 * pairs/arcs by min_prob.  Registers constraints and
	 * heuristics and then calls read_arcmatch_scores.  Constructs
	 * BasePairs objects for each single object and registers
	 * them.  Generates adjacency lists of arc matches for
	 * internal use and sorts them. Lists contain only valid arc
	 * matches according to constraints and heuristics (see
	 * is_valid_arcmatch()). The constructed object does not
	 * explicitely represent/maintain the scores of arc matchs.
	 
	 * @param rnadataA data for RNA A
	 * @param rnadataB data for RNA B
	 * @param min_prob consider only arcs with this minimal probability
	 * @param max_length_diff consider arc matches only up to maximal length difference
	 * @param trace_controller arc matches only due to trace controller
	 * @param constraints arc matches only due to constraints 
	 */
	ArcMatches(const RnaData &rnadataA, 
		   const RnaData &rnadataB,
		   double min_prob,
		   size_type max_length_diff, 
		   const MatchController &trace_controller,
		   const AnchorConstraints &constraints);
    
	//! clean up base pair objects
	~ArcMatches();
    
	// for the mea probabilistic consistency transformation, support to read and write the arcmatch scores
	// this allows in general to have user defined arc-match scores
    
	/**
	 * \brief Reads scores for arc matches 
	 *
	 * Reads from a file
	 * reads a list i j k l score, where
	 * score is the score for matching arcs (i,j) and (k,l).
	 * @param arcmatch_scores_file file containing arc match scores
	 * @param probability_scale if >=0 read probabilities and multiply them by probability_scale 
	 *
	 * @note All registered arc matches are valid (is_valid_arcmatch()).
	 */
	void read_arcmatch_scores(const std::string &arcmatch_scores_file, int probability_scale);
    
    
	/**
	 * write arc match scores to a file
	 * (this is useful after the scores are generated from base pair probabilities)
	 */
	void write_arcmatch_scores(const std::string &arcmatch_scores_file, const Scoring &scoring) const;
    
    
	//! returns the base pairs object for RNA A
	const BasePairs &
	get_base_pairsA() const {
	    return *bpsA;
	}

	//! returns the base pairs object for RNA B
	const BasePairs &
	get_base_pairsB() const {
	    return *bpsB;
	}

	//! true, if arc match scores are explicit (because they are read in from a list)
	bool
	explicit_scores() const {
	    return maintain_explicit_scores;
	}
    
	/**
	 * @brief Make arcmatch scores explicit
	 * 
	 * @param scoring An scoring object
	 */
	void
	make_scores_explicit(const Scoring &scoring);
	
	
	/**
	 * get the score of an arc match
	 * @param am arc match
	 * @returns score of arc match
	 * @pre object represents arc match scores explicitely
	 */
	score_t
	get_score(const ArcMatch &am) const {
	    assert(maintain_explicit_scores);
	    return scores[am.idx()];
	}
        
	//! total number of arc matches
	size_type num_arc_matches() const {
	    return number_of_arcmatches;
	}
    
	//! get arc match by its index
	const ArcMatch &arcmatch(size_type idx) const {
	    assert(idx<number_of_arcmatches);
	    return arc_matches_vec[idx];
	}
        
	// ============================================================
	// Iteration over arc matches
	//
    
	//! list of all arc matches that share the common right end (i,j)
	const ArcMatchIdxVec &
	common_right_end_list(size_type i, size_type j) const {
	    return common_right_end_lists(i,j);
	}
    
	//! list of all arc matches that share the common left end (i,j)
	const ArcMatchIdxVec &
	common_left_end_list(size_type i, size_type j) const {
	    return common_left_end_lists(i,j);
	}
    
    
	// ============================================================
    

	/** 
	 * @brief get the maximal right ends of any arc match with left ends (al,bl).
	 *
	 * @pre max_ar, max_br are initialized with smallest possible values
	 *
	 * @param al left and in sequence A
	 * @param bl left and in sequence B
	 * @param[in,out] max_ar maximal right end in sequence A
	 * @param[in,out] max_br maximal right end in sequence B
	 * @param no_lonely_pairs whether in lonely pair mode 
	 *
	 * @return out parameters max_ar, max_br are set to the maximal right ends
	 *
	 * Determine minimal max_ar and max_br such that there is no arc match with left ends
	 * al, bl and larger right ends ar>max_ar or al>max_br.
	 * In lonely pair mode, consider only arc matchs that have an immediately enclosing arc match.
	 */
	void get_max_right_ends(size_type al,size_type bl,size_type *max_ar,size_type *max_br, bool no_lonely_pairs) const; 

	
	/**
	 * get the minimal right ends of any arc match with left ends (al,bl).
	 * pre: min_ar, min_br are initialized with largest possible values
	 * returns result in out parameters min_ar, min_br
	 * @note Does not support noLP like get_max_right_ends() yet 
	 */
	void get_min_right_ends(size_type al,size_type bl,size_type *min_ar,size_type *min_br) const; 
    
	// ------------------------------------------------------------
	// inner arc matches

    
	/**
	 * whether there is an inner (valid) arc match for the
	 * arc with the given index
	 */
	bool
	exists_inner_arc_match(const ArcMatch &am) const {
	    return inner_arcmatch_idxs[am.idx()] < num_arc_matches();
	}

	/**
	 * tests, whether there is an inner arc match (then its arcs have non-zero probability)
	 */
	bool
	is_stackable(const ArcMatch &am) const {
	    return 
		exists_inner_arc_match(am);
	    
	    // previously there was an additional condition: 
	    // * and the arcs of the match have non-zero joint probability to occur simultaneously with the inner arc!
	    // WHY?
	    // 	&&
	    // 	bpsA->get_arc_2_prob(am.arcA().left(),am.arcA().right())>0 
	    // 	&& 
	    // 	bpsB->get_arc_2_prob(am.arcB().left(),am.arcB().right())>0;
	}

        
	/**
	 * index of inner arc match for the
	 * arc with the given index.
	 * Call only if there is such an arc match
	 */
	const ArcMatch &
	inner_arc_match(const ArcMatch &am) const {
	    return arcmatch(inner_arcmatch_idxs[am.idx()]);
	}
    
	/**
	 * sort the lists of arc matches with common right ends in "common_right_end_list"
	 * by their left ends in lexicographically descending order.
	 * Additionally, generate data structure for optimized traversal.
	 */
	void
	sort_right_adjacency_lists();

	// ------------------------------------------------------------
	// iteration (in no specific order)

	//! const iterator over arc matches
	typedef ArcMatchVec::const_iterator const_iterator;

	//! begin of arc matches vector
	const_iterator begin() const {return arc_matches_vec.begin();}

	//! end of arc matches vector
	const_iterator end() const {return arc_matches_vec.begin()+number_of_arcmatches;}
    
    };

    
    /**
     * @brief class ArcMatches with additional mapping
     * 
     * Like ArcMatches, maintain the relevant arc matches and their
     * scores. Additionally, build an index to support mapping from
     * pairs of arcs to arc matches.
     */
    class ArcMatchesIndexed : public ArcMatches {
    public:
	/**
	 *  \brief construct with explicit arc match score list
	 * 
	 * construct from seqnames and explicit list of all scored arc
	 * matches together with their score.
	 *
	 *
	 * @note registers constraints and heuristics and then calls read_arcmatch_scores.
	 * The constructed object explicitely represents/maintains the scores of arc matchs. 
	 *
	 * @param seqA_ sequence A
	 * @param seqB_ sequence B
	 * @param arcmatch_scores_file file containing arc match scores
	 * @param probability_scale if >=0 read probabilities and multiply them by probability_scale 
	 * @param max_length_diff accept arc matches only up to maximal length difference
	 * @param trace_controller accept only due to trace controller
	 * @param constraints accept only due to constraints 
	 */
	ArcMatchesIndexed(const Sequence &seqA_, 
			  const Sequence &seqB_, 
			  const std::string &arcmatch_scores_file,
			  int probability_scale,
			  size_type max_length_diff,
			  const MatchController &trace_controller,
			  const AnchorConstraints &constraints)
	    :ArcMatches(seqA_,seqB_,arcmatch_scores_file,probability_scale,max_length_diff,trace_controller,constraints),
	     am_index_()
	{
	    build_arcmatch_index();
	}
	
	/** 
	 *  \brief construct from single base pair probabilities. 
	 * 
	 * In this case, the object filters for relevant base
	 * pairs/arcs by min_prob.  Registers constraints and
	 * heuristics and then calls read_arcmatch_scores.  Constructs
	 * BasePairs objects for each single object and registers
	 * them.  Generates adjacency lists of arc matches for
	 * internal use and sorts them. Lists contain only valid arc
	 * matches according to constraints and heuristics (see
	 * is_valid_arcmatch()). The constructed object does not
	 * explicitely represent/maintain the scores of arc matchs.
	 
	 * @param rnadataA data for RNA A
	 * @param rnadataB data for RNA B
	 * @param min_prob consider only arcs with this minimal probability
	 * @param max_length_diff consider arc matches only up to maximal length difference
	 * @param trace_controller arc matches only due to trace controller
	 * @param constraints arc matches only due to constraints 
	 */
	ArcMatchesIndexed(const RnaData &rnadataA,
			  const RnaData &rnadataB,
			  double min_prob,
			  size_type max_length_diff,
			  const MatchController &trace_controller,
			  const AnchorConstraints &constraints)
	    :ArcMatches(rnadataA,rnadataB,min_prob,max_length_diff,trace_controller,constraints),
	     am_index_()
	{
	    build_arcmatch_index();
	}
	
    private:
	//! type of a pair of indices
	typedef std::pair<size_type,size_type> idx_pair_t;
	
	//! type of index for mapping arc pairs to arc matchs
	typedef std::tr1::unordered_map<idx_pair_t,ArcMatch::idx_type,pair_of_size_t_hash> am_index_type;
	
	//! index for mapping arc pairs to arc matchs
	am_index_type am_index_;
	
	/** 
	 * @brief construct index of arcmatch indices by arc pair indices
	 * 
	 */
	void
	build_arcmatch_index();

	
    public:
	
	/** 
	 * @brief the invalid arc match index 
	 * 
	 * @return invalid arc match index 
	 *
	 * @note this index is returned by am_index if the queried
	 * pair of arcs does not correspond to a valid arcmatch
	 */
	const ArcMatch::idx_type invalid_am_index() const {
	    // The invalid arc match index is implemented to be the
	    // maximum valid index + 1.
	    // This allows an efficient implementation, where we push 
	    // an invalid arc match to the end of vector arc_matches_vec.
	    // Thus, we return size-1!

	    return number_of_arcmatches;
	}
	
	/** 
	 * @brief Lookup arc match index by pair of arc indices
	 * 
	 * @param arcAIdx index of arc in A
	 * @param arcBIdx index of arc in B
	 * 
	 * @return index of arcmatch of arcA and arcB or invalid_am_index()
	 */
	const ArcMatch::idx_type
	am_index(const size_type &arcAIdx,const size_type &arcBIdx) const {
	    am_index_type::const_iterator it = am_index_.find(idx_pair_t(arcAIdx,arcBIdx));
	    if (am_index_.end() != it) {
		return it->second;
	    } else {
		return invalid_am_index();
	    }
	}

	/** 
	 * @brief Lookup arc match by pair of arcs
	 * 
	 * @param arcA arc in A
	 * @param arcB arc in B
	 * 
	 * @return arcmatch of arcA and arcB (or an invalid arc match with index invalid_am_index())
	 */
	const ArcMatch &
	am_index(const Arc &arcA,const Arc &arcB) const {
	    return arc_matches_vec[am_index(arcA.idx(),arcB.idx())];
	}
    };

} // end namespace LocARNA

#endif // LOCARNA_ARC_MATCHES_HH
