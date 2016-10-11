#ifndef LOCARNA_RNA_DATA_HH
#define LOCARNA_RNA_DATA_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include "aux.hh"
#include "sparse_matrix.hh"

namespace LocARNA {

    class MultipleAlignment;
    class Sequence;
    class Alignment;
    class RnaEnsemble;
    class RnaDataImpl;
    class PFoldParams;
    class SequenceAnnotation;
    
    /**
     * @brief represent sparsified data of RNA ensemble
     * 
     * knows sequence, cutoff probability and base pair probabilities
     * greater than the cutoff probability; potentially knows stacking
     * probabilities
     *
     * @note This class knows the cutoff probability *of its data*. This
     * cutoff can be different from the cutoff in classes like
     * BasePairs which defines the structure elements that are
     * considered by algorithms.
     *
     * @note the class guarantees that sequences are normalized
     * (uppercase, T->U) even when read in unnormalized form,
     * e.g. from file or stream
     */
    class RnaData {
    private:
	friend class RnaDataImpl;
	RnaDataImpl *pimpl_;  //!<- pointer to corresponding implementation object

    public:
	
	//! arc probability matrix
	typedef SparseMatrix<double> arc_prob_matrix_t;
	
	typedef size_t size_type; //!< usual size type
	
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rna_ensemble RNA ensemble data
	 * @param p_bpcut cutoff probability
	 * @param pfoldparams folding parameters (controls stacking)
	 *
	 * @note RnaData copies all required data from rna_ensemble
	 * and does not keep a reference; if pfoldparams.stacking is true,
	 * copy stacking terms
	 */
	RnaData(const RnaEnsemble &rna_ensemble,
		double p_bpcut,
		const PFoldParams &pfoldparams);

	/** 
	 * @brief Construct from file
	 * 
	 * @param filename input file name
	 * @param p_bpcut cutoff probability
	 * @param pfoldparams folding parameters
	 *
	 * @note autodetect format of input; 
	 * for fa or aln input formats, predict base pair probabilities 
	 *
	 * @todo consider to allow reading from istream; use
	 * istream::seekg(0) to reset stream to beginning (needed for
	 * format autodetect.) Is there a problem due to fail on
	 * eofbit in C++98?
	 */
	RnaData(const std::string &filename,
		double p_bpcut,
		const PFoldParams &pfoldparams);
	
	/** 
	 * @brief Construct as consensus of two aligned RNAs
	 * 
	 * @param rna_dataA RNA ensemble data A
	 * @param rna_dataB RNA ensemble data B
	 * @param alignment Alignment of A and B
	 * @param p_expA background probability for A
	 * @param p_expB background probability for B
	 * @param only_local if true, construct only local alignment
	 *
	 * The object uses mean cutoff probability of the given
	 * objects; The background probability is used in computing
	 * consensus probabilities.  If both input rna data objects
	 * have stacking probabilities, stacking consensus
	 * probabilities are computed as well. If the object contain
	 * sequence anchors, we construct the new object with a
	 * consensus anchor string. (The latter is done as part of the
	 * consensus sequence computation.)
	 */
	RnaData(const RnaData &rna_dataA,
		const RnaData &rna_dataB,
		const Alignment &alignment,
		double p_expA,
		double p_expB, 
		bool only_local=false
		);
	
    protected:
    	/** 
	 * @brief Almost empty constructor
	 * 
	 * @param p_bpcut cutoff probability
	 */
	RnaData(double p_bpcut);
	
    private:
	/**
	 * @brief copy constructor
	 */
	RnaData(const RnaData &);
    public:
	
	/**
	 * @brief destructor 
	 */
	virtual 
	~RnaData();

    private:
	/**
	 * @brief assignment operator
	 */
	RnaData &
	operator =(const RnaData &);

    public:

	/**
	 * @brief Get the multiple alignment as sequence
	 * @return sequence
	 */
	const Sequence &
	sequence() const;

	/**
	 * @brief Get the multiple alignment
	 * @return multiple alignment
	 */
	const MultipleAlignment &
	multiple_alignment() const;

	/**
	 * @brief Get the sequence length
	 * @return length of RNA sequence
	 */
	size_type
	length() const;

	/**
	 * @brief Get base pair cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	arc_cutoff_prob() const;
	
	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return probability p_ij of basepair (i,j) if
	 * p_ij>p_bpcut; otherwise, 0
	 */
	double 
	arc_prob(pos_type i, pos_type j) const;

    protected:
	//! type of constant iterator over arcs with probability above cutoff
	typedef arc_prob_matrix_t::const_iterator arc_probs_const_iterator;
	
	/**
	 * @brief begin of arcs with probability above cutoff
	 * Supports iteration over arcs
	 * @returns constant iterator
	 */
	arc_probs_const_iterator
	arc_probs_begin() const;

	/**
	 * @brief begin of arcs with probability above cutoff
	 * Supports iteration over arcs
	 * @returns constant iterator
	 */
	arc_probs_const_iterator
	arc_probs_end() const;
    public:
	

	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return joint probability p^(2)_ij of basepair (i,j) and
	 * (i+1,j-1) if p_i+1j-1>p_bpcut and p^(2)_ij > p_bpcut; otherwise, 0
	 */
	double 
	joint_arc_prob(pos_type i, pos_type j) const;
	
	/**
	 * @brief Get arc probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return conditional probability p_ij|i+1j-1 of basepair
	 * (i,j) under condition of base pair (i+1,j-1) if
	 * p_i+1j-1>p_bpcut and p^(2)_ij > p_bpcut; throw exception if
	 * p_i+1j-1<=p_bpcut; otherwise, 0
	 */
	double 
	stacked_arc_prob(pos_type i, pos_type j) const;
	
	// some computed probabilities (for convenience)
	/**
	 * \brief Probability that a position is paired upstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j>i (upstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_downstream
	 */
	double
	prob_paired_upstream(pos_type i) const;
        
	/**
	 * \brief Probability that a position is paired downstream
	 * 
	 * \param i sequence position
	 * \return probability that a position i is paired with a position j<i (downstream)
	 * @note O(sequence.length()) implementation
	 * @see prob_paired_upstream
	 */
	double
	prob_paired_downstream(pos_type i) const;
    
	/**
	 * \brief Unpaired probability 
	 * \param i sequence position
	 * \return probability that a position i is unpaired
	 * @note O(sequence.length()) implementation
	 */
	double
	prob_unpaired(pos_type i) const;
	
	// IO
	/** 
	 * Write data in pp format
	 * 
	 * @param out output stream
	 * @param p_outbpcut cutoff probability
	 *
	 * @return stream
	 *
	 * Writes only base pairs with probabilities greater than
	 * p_outbpcut
	 */
	std::ostream &
	write_pp(std::ostream &out, double p_outbpcut=0) const;

	/**
	 * @brief Write object size information
	 *
	 * @param out output stream
	 *
	 * Writes numbers of stored probabilities to stream
	 */
	std::ostream &
	write_size_info(std::ostream &out) const;

	/** 
	 * @brief Availability of stacking terms 
	 * 
	 * @return whether stacking terms are available
	 */
	bool
	has_stacking() const;

	/**
	 * @brief Write access to alignment anchors
	 * @param anchors alignment anchors
	 * @see MultipleAlignment::set_annotation()
	 */
	void
	set_anchors(const SequenceAnnotation &anchors);

    protected:
	
	/** 
	 * @brief initialize from fixed structure
	 * 
	 * @param structure fixed structure
	 * @param stacking whether to initialize stacking terms
	 *
	 * @note can be overloaded to initialize with additional
	 * information (in loop probabilities)
	 */
	virtual
	void
	init_from_fixed_structure(const SequenceAnnotation &structure,
				  bool stacking);

	/** 
	 * @brief initialize from rna ensemble 
	 * 
	 * @param rna_ensemble rna ensemble
	 * @param stacking whether to initialize stacking terms
	 * 
	 * @note can be overloaded to initialize with additional
	 * information (in loop probabilities)
	 */
	virtual
	void
	init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
			       bool stacking);
	

	/** 
	 * @brief read and initialize from file, autodetect format
	 * 
	 * @param filename name of input file
	 * @param stacking whether to read stacking terms
	 *
	 * @return whether probabilities were read completely
	 *
	 * @note: this method is designed such that it can be used for
	 * RnaData and ExtRnaData
	 */
	bool
	read_autodetect(const std::string &filename,
			bool stacking);
		
	/**
	 * @brief check in loop probabilities
	 *
	 * @return true iff loop probabilities are available or not
	 * required 
	 * @note use to indicate the need for recomputation
	 * in read_autodetect(); always true in RnaData
	 */
	virtual
	bool
	inloopprobs_ok() const {return true;}
		
	
	/** 
	 * Read data in pp format 2.0
	 * 
	 * @param filename name of input file
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_bpcut_; reads stacking probabilities only if
	 * has_stacking_ is true
	 *
	 * @note can be overloaded to read extension sections
	 *
	 * @note pp is a proprietary format of LocARNA. In its
	 * simplest version, it starts with the sequence/alignment and
	 * then simply lists the arcs (i,j) with their probabilities p
	 * and optionally stacking probabilities p2.  pp-files contain
	 * entries i j p [p2]. p denotes the probabilitiy of base pair
	 * (i,j). The optional stacking probability p2 is the joint
	 * probability of base pairs (i,j) and (i+1,j+1).
	 *
	 * @note handling of stacking: after the call, has_stacking_
	 * is true only if the file specified the STACK keyword and
	 * has_stacking_ was true before.
	 */
	virtual
	void
	read_pp(const std::string &filename);
	
	/**
	 * Read data in pp format 2.0
	 * 
	 * @param in input stream
	 *
	 * @see read_pp(std::string)
	 */
	virtual
	std::istream &
	read_pp(std::istream &in);

	/** 
	 * Read data in the old pp format
	 * 
	 * @param filename name of input file 
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_bpcut_; reads stacking probabilities only if
	 * has_stacking_ is true
	 *
	 * @note the old pp format starts with the sequence/alignment
	 * and then simply lists the arcs (i,j) with their
	 * probabilities p and optionally stacking probabilities p2.
	 * pp-files contain entries i j p [p2]. p denotes the
	 * probabilitiy of base pair (i,j). The optional stacking
	 * probability p2 is the joint probability of base pairs (i,j)
	 * and (i+1,j+1).
	 *
	 * @note handling of stacking: after the call, has_stacking_
	 * is true only if the file specified at least one stacking
	 * probability and has_stacking_ was true before.
	 */
	void
	read_old_pp(const std::string &filename);

	/** 
	 * Read data in Vienna's dot plot ps format
	 * 
	 * @param filename name of input file
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_bpcut_; reads stacking probabilities only if
	 * has_stacking_ is true
	 *
	 * @note Recently changed behavior: reads sequence name from
	 * file (instead of guessing from filename!); stacking
	 * probabilities are read if available (then, sets
	 * has_stacking_ to true)
	 *
	 * @note throws wrong_format exception if not in ps format 
	 *
	 * @note reading dot plot ps format is extremely fragile since
	 * the dot plot ps format was designed for viewing not for
	 * data input 
	 */
	void
	read_ps(const std::string &filename);

	
    }; // end class RnaData
  
}


#endif // LOCARNA_RNA_DATA_HH
