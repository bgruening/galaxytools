#ifndef LOCARNA_EXT_RNA_DATA_HH
#define LOCARNA_EXT_RNA_DATA_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include "aux.hh"
#include "sparse_matrix.hh"
#include "rna_data.hh"

namespace LocARNA {

    class RnaEnsemble;
    class ExtRnaDataImpl;
    class PFoldParams;

    /**
     * @brief represent sparsified data of RNA ensemble extended by in loop probabilities
     * 
     * knows sequence, cutoff probability and base pair
     * probabilities greater than the cutoff probability
     *
     * @note This class knows the cutoff probabilities of its
     * data. These cutoffs can be different from the cutoffs in classes
     * like BasePairs which define the structure elements that are
     * considered by algorithms.
     * 
     * @note put the extension of extended pp format as additional
     * record after the regular pp-data, such that extended pp files
     * can be read as non-extended ones
     */
    class ExtRnaData: public RnaData {
    private:
	friend class ExtRnaDataImpl;
	ExtRnaDataImpl *pimpl_;  //!<- pointer to corresponding implementation object
    public:
	/** 
	 * @brief Construct from RnaEnsemble with cutoff probability
	 * 
	 * @param rna_ensemble RNA ensemble data
	 * @param p_bpcut cutoff probability for base pairs
	 * @param p_bpilcut cutoff probability for base pairs in loops
	 * @param p_uilcut cutoff probability for unpaired bases in loops
	 * @param pfoldparams parameters for partition folding
	 *
	 *
	 * @note requires that rnaensemble has in loop probabilities
	 */
	ExtRnaData(const RnaEnsemble &rna_ensemble,
		   double p_bpcut,
		   double p_bpilcut,
		   double p_uilcut,
		   const PFoldParams &pfoldparams);
	
	/** 
	 * @brief Construct from input file
	 * 
	 * @param filename input file name
	 * @param p_bpcut cutoff probability for base pairs
	 * @param p_bpilcut cutoff probability for base pairs in loops
	 * @param p_uilcut cutoff probability for unpaired bases in loops
	 * @param pfoldparams parameters for partition folding
	 *	 
	 * @note autodetect format of input; 
	 * for fa or aln input formats, predict base pair probabilities 
	 */
	ExtRnaData(const std::string &filename,
		   double p_bpcut,
		   double p_bpilcut,
		   double p_uilcut,
		   const PFoldParams &pfoldparams);

    private:
	/**
	 * @brief copy constructor
	 */
	ExtRnaData(const ExtRnaData &);
    public:

	/** 
	 * @brief destructor 
	 */
	~ExtRnaData();
	
    private:
	/**
	 * @brief assignment operator
	 */
	ExtRnaData &
	operator =(const ExtRnaData &);
    public:

	/**
	 * @brief Get base pair in loop cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	arc_in_loop_cutoff_prob() const;
	
	/**
	 * @brief Get base pair in loop probability
	 * @param i left sequence position of inner base pair
	 * @param j right sequence position of inner base pair
	 * @param p left sequence position of closing base pair
	 * @param q right sequence position of closing base pair
	 *
	 * @return joint probability of basepair (i,j) in loop of (p,q) if
	 * above threshold; otherwise, 0
	 */
	double 
	arc_in_loop_prob(pos_type i, pos_type j,pos_type p, pos_type q) const;
	
	/**
	 * @brief Get base pair in loop probability
	 * @param i left sequence position  
	 * @param j right sequence position
	 *
	 * @return probability of basepair (i,j) in external loop if
	 * above threshold; otherwise, 0
	 */
	double 
	arc_external_prob(pos_type i, pos_type j) const;
	

	/**
	 * @brief Get unpaired base in loop cutoff probability
	 * @return cutoff probability p_bpcut
	 */
	double
	unpaired_in_loop_cutoff_prob() const;
	
	/**
	 * @brief Get base pair in loop probability
	 * @param k sequence position (unpaired base)  
	 * @param p left sequence position of closing base pair  
	 * @param q right sequence position of closing base pair
	 *
	 * @return joint probability of unpaired base k in loop of (p,q) if
	 * above threshold; otherwise, 0
	 */
	double 
	unpaired_in_loop_prob(pos_type k,pos_type p, pos_type q) const;

	/**
	 * @brief Get base pair in loop probability
	 * @param k sequence position (unpaired base)
	 *
	 * @return probability of base k unpaired in external loop if
	 * above threshold; otherwise, 0
	 */
	double 
	unpaired_external_prob(pos_type k) const;

	
	/**
	 * @brief Write object size information
	 *
	 * @param out output stream
	 *
	 * Writes numbers of stored probabilities to stream
	 */
	std::ostream &
	write_size_info(std::ostream &out) const;
	
	
	// IO
	/** 
	 * Write data in extended pp format
	 * 
	 * @param out output stream
	 * @param p_outbpcut cutoff probability
	 * @param p_outbpilcut cutoff probability base pairs in loop
	 * @param p_outuilcut cutoff probability unpaired in loop
	 *
	 * @return stream
	 *
	 * Writes only base pairs with probabilities greater than
	 * p_outbpcut; base pairs in loops, p_outbpilcut; unpaired
	 * bases in loops, p_outuilcut
	 */
	std::ostream &
	write_pp(std::ostream &out,
		 double p_outbpcut=0,
		 double p_outbpilcut=0,
		 double p_outuilcut=0) const;

    protected:

	/** 
	 * Read data in pp format 2.0 with in-loop probabilities
	 * 
	 * @param in input stream
	 *
	 * @see RnaData::read_pp()
	 */
	virtual
	std::istream &
	read_pp(std::istream &in);

	/** 
	 * @brief initialize from fixed structure
	 * 
	 * @param structure fixed structure
	 * @param stacking whether to initialize stacking terms
	 *
	 * @note overloaded to initialize with additional information
	 * (in loop probabilities)
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
	 * @note overloaded to initialize with additional information
	 * (in loop probabilities)
	 *
	 * @pre bp cutoff prob is initialized 
	 */
	virtual
	void
	init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
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
	inloopprobs_ok() const;
		
    }; // end class ExtRnaData
  
}

#endif // LOCARNA_EXT_RNA_DATA_HH
