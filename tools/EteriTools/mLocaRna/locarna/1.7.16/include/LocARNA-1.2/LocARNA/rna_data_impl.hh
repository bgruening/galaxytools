#ifndef LOCARNA_RNA_DATA_IMPL_HH
#define LOCARNA_RNA_DATA_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include "rna_data.hh"
#include "sequence.hh"

namespace LocARNA {

    class MultipleAlignment;
    class RnaEnsemble;
    class PFoldParams;
    //    template<class T> class SparseVector<T>;

    /**
     * @brief Implementation of RnaData
     */
    struct RnaDataImpl {
    
	//! type for matrix of arc probabilities
	typedef RnaData::arc_prob_matrix_t arc_prob_matrix_t;
	
	RnaData *self_; //!<- pointer to corresponding non-impl object

	//! the sequence
	MultipleAlignment sequence_; 

	//! cutoff probabilitiy for base pair
	double p_bpcut_;
	
	/**
	 * sparse array for all arc probabilities above threshold; the
	 * array is used when reading in the probabilities and for
	 * merging probs during pp-output
	 */
	arc_prob_matrix_t arc_probs_; 
	
	/**
	 * sparse array for all probabilities that a pair (i,j) and
	 * its immediately inner pair (i+1,j-1) are formed
	 * simultaneously above threshold; analogous to arc_probs_
	 *
	 * @note arc_2_probs_ has entry (i,j) implies arc_probs_ has entry (i,j)
	 */
	arc_prob_matrix_t arc_2_probs_; 
	
	//! whether stacking probabilities are available
	bool has_stacking_; 
	
	/** 
	 * @brief Construct as consensus of two aligned RNAs
	 * 
	 * @param self pointer to corresponding RnaData object
	 * @param rna_dataA data of RNA A 
	 * @param rna_dataB data of RNA B
	 * @param alignment pairwise alignment of A and B
	 * @param p_expA background probability A
	 * @param p_expB background probability B
	 */
	RnaDataImpl(RnaData *self,
		    const RnaData &rna_dataA,
		    const RnaData &rna_dataB,
		    const Alignment::edges_t &alignment,
		    double p_expA,
		    double p_expB);

    	/** 
	 * @brief Almost empty constructor
	 * 
	 * @param self pointer to corresponding RnaData object
	 * @param p_bpcut cutoff probability
	 */
	RnaDataImpl(RnaData *self,
		    double p_bpcut);
	
	// ----------------------------------------
	// METHODS


	/** 
	 * @brief initialize from fixed structure
	 * 
	 * @param structure fixed structure
	 * @param stacking whether to initialize stacking terms
	 */
	void
	init_from_fixed_structure(const SequenceAnnotation &structure,
				  bool stacking);

	/** 
	 * @brief initialize from rna ensemble 
	 * 
	 * @param rna_ensemble rna ensemble
	 * @param stacking whether to initialize stacking terms
	 */
	void
	init_from_rna_ensemble(const RnaEnsemble &rna_ensemble,
			       bool stacking);

	/**
	 * @brief read sequence section of pp-format
	 *
	 * @param in input stream
	 * @return stream
	 *
	 * this section comprises sequence/multiple alignment
	 * (possibly including sequence anchor annotation)
	 */
	std::istream &
	read_pp_sequence(std::istream &in);

	/**
	 * @brief read section of base pair probabilities of pp-format
	 *
	 * @param in input stream
	 *
	 * Reads only base pairs with probabilities greater than
	 * p_bpcut_; reads stacking only if has_stacking_
	 */
	std::istream &
	read_pp_arc_probabilities(std::istream &in);
	
	/**
	 * @brief write section of base pair probabilities of pp-format
	 *
	 * @param out ouput stream
	 * @return stream
	 *
	 */
	std::ostream &
	write_pp_sequence(std::ostream &out) const;

	/**
	 * @brief write section of base pair probabilities of pp-format
	 *
	 * @param out ouput stream
	 * @param p_outbpcut cutoff probabilitiy
	 * @param stacking whether to write stacking probabilities; if
	 *   stacking but !has_stacking_, no stacking terms are
	 *   written but flag #STACKS is written to output
	 *
	 * @return stream
	 *
	 * Write only base pairs with probabilities greater than
	 * p_outbpcut
	 */
	std::ostream &
	write_pp_arc_probabilities(std::ostream &out,
				   double p_outbpcut,
				   bool stacking) const;


	/** 
	 * @brief Initialize as consensus of two aligned RNAs
	 * 
	 * @param edges alignment edges
	 * @param rna_dataA rna data A
	 * @param rna_dataB rna data B
	 * @param p_expA background probability A
	 * @param p_expB background probability B
	 * @param stacking if true, stacking consensus is computed
	 */
	void
	init_as_consensus_dot_plot(const Alignment::edges_t &edges,
				   const RnaData &rna_dataA,
				   const RnaData &rna_dataB,
				   double p_expA,
				   double p_expB,
				   bool stacking
				   );
	
	/** 
	 * @brief Consensus probability
	 * 
	 * @param pA probability A
	 * @param pB probability B
	 * @param sizeA number of rows in sequence A
	 * @param sizeB number of rows in sequence B
	 * @param p_expA background probability A
	 * @param p_expB background probability B
	 * 
	 * @pre p_bpcut_ is initialized
	 *
	 * Essentially computes a weighted geometric mean; some care
	 * is taken, to avoid total extinction and to restrict the
	 * accumulation of (small) probabilities.
	 * 
	 * @return consensus probability
	 */
	double
	consensus_probability(double pA,
			      double pB, 
			      size_t sizeA,
			      size_t sizeB,
			      double p_expA,
			      double p_expB) const;

    }; // end class RnaDataImpl
    

} //end namespace LocARNA


#endif // LOCARNA_RNA_DATA_IMPL_HH
