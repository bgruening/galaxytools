#ifndef LOCARNA_STRAL_SCORE_HH
#define LOCARNA_STRAL_SCORE_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <math.h>


#include "aux.hh"
#include "sequence.hh"

namespace LocARNA {

    template <class T> class Matrix;
    template <class T> class Alphabet;
    class RnaData;

    //! \brief Implements the stral-like scoring function
    class StralScore {
    
	typedef std::vector<double> p_vec_t;
    
	Sequence seqA;
	Sequence seqB;
    
	p_vec_t p_upA;   //!< probability paired upstream seq A
	p_vec_t p_downA; //!< probability paired downstream seq A
	p_vec_t p_unA;   //!< probability unpaired seq A
    
	p_vec_t p_upB;   //!< probability paired upstream seq B
	p_vec_t p_downB; //!< probability paired downstream seq B
	p_vec_t p_unB;   //!< probability unpaired seq B
    
	const Matrix<double> &sim_mat;
	const Alphabet<char> &alphabet;
	double pf_struct_weight;
	double gap_opening;
	double gap_extension;
    
    private:
	void init_prob_vecs(const RnaData &rna,
			    p_vec_t &p_up,
			    p_vec_t &p_down,
			    p_vec_t &p_un);
    public:
    
	/** 
	 * Construct for pair of RNAs with parameters for alignment
	 * 
	 * @param rnaA data of first RNA
	 * @param rnaB data of second RNA
	 * @param sim_mat_ similarity matrix for bases
	 * @param alphabet_ alphabet
	 * @param pf_struct_weight_ structure weight 
	 * @param gap_opening_ gap opening cost
	 * @param gap_extension_ gap extension cost
	 */
	StralScore(const RnaData &rnaA,
		   const RnaData &rnaB, 
		   const Matrix<double> &sim_mat_,
		   const Alphabet<char> &alphabet_,
		   double pf_struct_weight_,
		   double gap_opening_,
		   double gap_extension_
		   );

	/** 
	 * \brief Compute STRAL-like similarity of two residues in the two RNAs
	 *
	 * @param i position in sequence A
	 * @param j position in sequence B
	 *

	 * @note Computes the average similarity over all pairs of
	 * alignment rows in the RNA sequence, which are alignments in
	 * general.
	 *
	 * @note The treatment of gaps and unknown nucleotide symbols
	 * in the aligned alignments is quite ad hoc.
	 * 
	 * @return similarity of residues i in A and j in B. 
	 */
	double sigma(size_type i, size_type j) const;
    
	/** 
	 * \brief Read gap opening cost
	 * 
	 * @return gap opening cost
	 */
	double alpha() const {
	    return gap_opening;
	}

	/** 
	 * \brief Read gap extension cost
	 * 
	 * @return gap extension cost
	 */
	double beta() const {
	    return gap_extension;
	}
    
	//! \brief Reverse the scoring
	//!
	//! @post the object scores the reverted RNAs
	void 
	reverse();
    
    };

}

#endif //LOCARNA_STRAL_SCORE_HH
