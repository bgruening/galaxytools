#ifndef LOCARNA_BASEPAIRS_HH
#define LOCARNA_BASEPAIRS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include<iosfwd>

#include <vector>
#include <set>
#include <assert.h>

#include "params.hh"
#include "sparse_matrix.hh"

namespace LocARNA {

    class RnaData;
    class Sequence;

    /**
     * \brief Represents a base pair
     * 
     * stores base pair left end, right end
     * and an arc index
     * 
     * @note index uniqueness is not guaranteed by the class
     * itself but can be used to locate the arc in a vector
     * by the caller
     *
     * @note To indicate the relation to the BasePair class, we define
     * with prefix BasePairs__.  Using a nested class does not work
     * properly (C++ does not support forward references to nested
     * classes.) Using a namespace instead of prefix should work as well,
     * but caused problems at first when including the adjecency entry classes in the namespace!
     *
     */
    class BasePairs__Arc {
    private:
	size_t idx_;
	size_t left_;
	size_t right_;
    public:
    
	/** 
	 * Construct from member values
	 * 
	 * @param idx  Index
	 * @param left  Left position of arc
	 * @param right  Right position of arc
	 */
	BasePairs__Arc(size_t idx, size_t left, size_t right):
	    idx_(idx),
	    left_(left),
	    right_(right)
	{}
	
	/*
	 * @brief Virtual destructor
	 */
	virtual
	~BasePairs__Arc();

	/** 
	 * Read access
	 *  
	 * @return left arc end
	 */
	size_t
	left() const {return left_;}

	/** 
	 * Read access
	 *  
	 * @return right arc end
	 */
	size_t
	right() const {return right_;}
	
	/** 
	 * Read access
	 *  
	 * @return index of arc
	 */
	size_t
	idx() const {return idx_;}
    };

    // ============================================================
    /**
     * @brief Describes sequence and structure ensemble of an RNA
     *
     * Stores and maintains the list of potential base pairs together
     * with their score contributions in arc matches.
     *
     * In contrast to RnaData, which stores the raw data for an RNA,
     * a BasePairs object knows about sparsification of base pairs by
     * a probability threshold and provides traversal of base pairs
     * suited for alignment algorithms.
     *
     * If the base pairs object is constructed from an RnaData object,
     * it knows about its corresponding RnaData object.
     *
     */
    class BasePairs
    {
    private:
	const RnaData *rna_data_;
	double min_prob_;
	double len_;

    public:
	typedef size_t size_type; //!< size

	typedef BasePairs__Arc Arc; //!< arc
	
	/**
	 * @brief Entry in a left adjacency list
	 *
	 * @see RightAdjEntry
	 *
	 * @note Deriving the class(es) from Arc is not equivalent to a
	 * type definition, like typedef Arc LeftAdjEntry, since it
	 * supports overloading of the operator < for left and right
	 * adjacency list entries (in contrast to typedef!).
	 * Damn you C++ :).
	 */
	class LeftAdjEntry : public Arc {
	public:
	    /** 
	     * Construct from arc 
	     * 
	     * @param a arc
	     */
	    LeftAdjEntry(const Arc &a): Arc(a) {}
	};
	
	/**
	 * @brief Entry in a right adjacency list
	 * @see LeftAdjEntry
	 */
	class RightAdjEntry : public Arc {
	public:
	    /** 
	     * Construct from arc 
	     * 
	     * @param a arc
	     */
	    RightAdjEntry(const Arc &a): Arc(a) {}
	};

	//! Vector of arcs
	typedef std::vector<Arc> arc_vec_t;
	
	/* types for data structures for the access of an arc in the structure,
	   by its right end, its left end, or left and right end 
	
	   the access structure is implemented as a hash map,
	*/
	
	//! type of left adjacency list
	typedef std::vector<LeftAdjEntry> LeftAdjList; 
	
	//! type of right adjacency list
	typedef std::vector<RightAdjEntry> RightAdjList; 
	
	//! type for matrix of arcs (actually arc indices)
	typedef SparseMatrix<int> arc_matrix_t;

	//! type for pair of positions (base pairs)
	typedef std::pair<size_type,size_type> bpair_t;
	
	//! type for set of position pairs
	typedef std::set<bpair_t> bpair_set_t;
    
    private:
	std::vector<LeftAdjList> left_;
	std::vector<RightAdjList> right_;

	arc_vec_t arc_vec_;
	arc_matrix_t arcs_;
    
	/**
	 * @brief resize data structures
	 * @param seq_len length of sequence
	 *
	 * Resizes the data structures left_ and right_ that allow
	 * fast acces to an arc arcs by its left end, its right end,
	 * or its left and right end
	 */
	void
	resize(size_type seq_len);
	
	//! generate the datastructures that allow fast access to arcs
	void
	generateBPLists(const RnaData &rna_data);
    
	//! sort the adjacency lists as expected by the alignment algorithm
	void
	sortAdjLists();
    
    public:
	
	/** 
	 * Construct from rna data
	 * 
	 * @param rna_data rna data
	 * @param min_prob minimal probability for filtering base pairs
	 *
	 * @note while rna data maintains base pairs regardless of
	 * their probability, an object of BasePairs represents the
	 * sparsified set of base pairs (due to min_prob_)
	 */
	BasePairs(const RnaData *rna_data,double min_prob):
	    rna_data_(rna_data),
	    min_prob_(min_prob),
	    len_(get_length_from_rna_data()),
	    left_(),
	    right_(),
	    arc_vec_(),
	    arcs_(-1)
	{
	    generateBPLists(*rna_data_);
	}

	/** 
	 * \brief Construct from a set of base pairs
	 * 
	 * @param len length of sequence
	 * @param bps set of base pairs
	 */
	BasePairs(size_type len, const bpair_set_t &bps ):
	    rna_data_(0),
	    min_prob_(1.0),
	    len_(len),
	    left_(),
	    right_(),
	    arc_vec_(),
	    arcs_(-1)
	{
	    resize(seqlen());
	    for (bpair_set_t::const_iterator it=bps.begin(); bps.end()!=it; ++it) {
		register_arc(it->first,it->second);
	    }
	    sortAdjLists();
	}
	
	// /**
	//  * @brief Copy constructor
	//  */
	// BasePairs(const BasePairs &bps);

	// /**
	//  * @brief Assignment operator
	//  */
	// BasePairs &
	// operator =(const BasePairs &bps);

	/**
	 * registers a basepair (i,j),
	 * maintains the basepair access data structures
	 */
	void
	register_arc(int i, int j);

	//! returns the list of arcs with right end i
	const LeftAdjList &
	left_adjlist(int i) const {
	    //std::cout<<"size of left adjlist of "<<i<<"is "<<left_[i].size()<<std::endl; 
	    return left_[i];
	}
	
	//! returns the list of arcs with left end i
	const RightAdjList &
	right_adjlist(int i) const { return right_[i];}

	//! accesses basepair by (i,j)
	const Arc &
	arc(int i,int j) const {return arc_vec_[arcs_(i,j)];}
    
	/**
	 * \param idx an arc index
	 * \returns arc with index idx
	 */
	const Arc &
	arc(size_type idx) const {
	    assert(idx<arc_vec_.size());
	    return arc_vec_[idx];
	}
        
	//! returns whether basepair (i,j) exists
	bool 
	exists_arc(int i,int j) const {return -1 != arcs_(i,j);}
    
	//! returns number of basepairs in the object
	size_type
	num_bps() const {return arc_vec_.size();}
    
	//! returns length of sequence
	size_type
	seqlen() const;
	
	double
	prob_min() const; //!< return minimal probability
    
	// /** 
	//  * @brief Access to corresponding RnaData object
	//  * 
	//  * @return reference to RnaData object
	//  */
	// const RnaData &
	// get_rna_data() const {
	//     assert(rna_data_!=NULL);
	//     return *rna_data_;
	// }
	

    private:
	
	// return length from rna data
	// pre: rna data available
	size_type
	get_length_from_rna_data() const;

    };

    std::ostream &
    operator <<(std::ostream &out, const BasePairs::Arc &arc);



}
#endif
