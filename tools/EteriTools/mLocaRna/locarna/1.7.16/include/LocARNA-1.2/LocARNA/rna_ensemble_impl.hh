#ifndef LOCARNA_RNA_ENSEMBLE_IMPL_HH
#define LOCARNA_RNA_ENSEMBLE_IMPL_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "rna_ensemble.hh"
#include "multiple_alignment.hh"
#include "sparse_matrix.hh"

#ifdef HAVE_LIBRNA	
#  include "mcc_matrices.hh"
#endif

namespace LocARNA {

    /**
     * @brief Implementation of RnaEnsemble
     */
    struct RnaEnsembleImpl {

	// we don't use the entire rna ensemble implementation if lib rna is not available
#ifdef HAVE_LIBRNA	
	
	//RnaEnsemble *self_; //!<- pointer to corresponding RnaEnsemble object
    
	MultipleAlignment sequence_; //!< the sequence
	
	//! whether pair probabilities are availabe
	bool pair_probs_available_; 
	
	//! whether stacking probabilities are available
	bool stacking_probs_available_; 
	
	//! whether "in loop" probabilities are availabe
	bool in_loop_probs_available_; 
		
	// std::vector<FLT_OR_DBL> qm1; // store qm1 for debugging
	std::vector<FLT_OR_DBL> qm2_;     //!< matrix qm2_ (stored VRNA-style in a vector)
	std::vector<FLT_OR_DBL> scale_;   //!< table for precomputed scaling of pf values
	std::vector<FLT_OR_DBL> expMLbase_; //!< table for precomputed multi loop terms
	
	McC_matrices_base *McCmat_; //!< DP matrix data structures of VRNA's McCaskill algorithm
	
	//! whether alifold was used to compute the McCaskill matrices
	bool used_alifold_;

	double min_free_energy_; //!< minimum free energy (if computed anyway)
	std::string min_free_energy_structure_; //!< minimum free energy structure (if computed)

	
	// /** 
	//  * @brief Construct from file
	//  * 
	//  * @param file file name
	//  * @param readPairProbs whether to read pair probabilities
	//  * @param readStackingProbs whether to read stacking probabilities
	//  * @param readInLoopProbs whether to read in loop probabilities
	//  */
	// RnaEnsembleImpl(const std::string &file,
	// 		bool readPairProbs,
	// 		bool readStackingProbs,
	// 		bool readInLoopProbs);
	
	/** 
	 * @brief Construct from sequence or multiple alignment
	 * 
	 * @param sequence  sequence or multiple alignment
	 * @param pfparams partition folding parameters
	 * @param inLoopProbs whether to compute in loop probabilities
	 * @param use_alifold whether to use alifold (required unless
	 * sequence is a single sequence)
	 */
	RnaEnsembleImpl(const MultipleAlignment &sequence,
			const PFoldParams &pfparams,
			bool inLoopProbs, 
			bool use_alifold);
	
	~RnaEnsembleImpl();

	////////////////////////////////////////////////////////////
	
	/** 
	 * @brief Pair type of an admissible basepair.
	 * 
	 * @param i left end of base pair
	 * @param j right end of base pair
	 * 
	 * @return pair type unless the base pair is not admissible,
	 * i.e. it is not complementary or has probability 0.0. Then
	 * return 0.
	 */
	int ptype_of_admissible_basepair(size_type i,size_type j) const;	

	/** 
	 * \brief (re)compute the pair probabilities
	 * 
	 * @param params pfolding parameters
	 * @param inLoopProbs whether in loop probabilities should be made available
	 * @param use_alifold whether alifold should be used
	 *
	 * @pre unless use_alifold, sequence row number has to be 1
	 */
	void
	compute_ensemble_probs(const PFoldParams &params,bool inLoopProbs, bool use_alifold);
	

	/**
	 * \brief Get joint probability of stacked arcs
	 * @param i left sequence position  
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	*/
	double
	arc_2_prob_noali(size_type i, size_type j) const;

	/**
	 * \brief Get joint probability of stacked arcs (alifold version)
	 * @param i left sequence position 
	 * @param j right sequence position
	 * \return probability of basepairs (i,j) and (i+1,j-1) occuring simultaneously
	 */
	double
	arc_2_prob_ali(size_type i, size_type j) const;
	
	/** 
	 * \brief Unpaired probabilty of base in a specified loop (alifold) 
	 *
	 * alifold-specific version of prob_unpaired_in_loop()
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * @see prob_unpaired_in_loop()
	 *
	 * @note pre: loop probs available, alifold used
	 */
	double
	unpaired_in_loop_prob_ali(size_type k,
				  size_type i,
				  size_type j) const;
	
	/** 
	 * \brief Unpaired probabilty of base in a specified loop (no alifold) 
	 *
	 * single sequence folding-specific version of prob_unpaired_in_loop()
	 *
	 * @param k unpaired sequence position
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that k is unpaired in the loop closed by i and j
	 *
	 * @see prob_unpaired_in_loop()
	 *
	 * @note pre: in loop probs are available, alifold not used
	 */
	double
	unpaired_in_loop_prob_noali(size_type k,size_type i,size_type j) const;
	
	/** 
	 * \brief Probabilty of base pair in a specified loop (alifold)
	 *
	 *`alifold-specific code
	 * 
	 * @param ip left end of inner base pair
	 * @param jp right end of inner base pair
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that (ip,jp) is inner base pair in the loop closed by i and j
	 *
	 * @see prob_basepair_in_loop()
	 *
	 * @note pre: loop probs available, alifold used
	 */
	double
	arc_in_loop_prob_ali(size_type ip,
			     size_type jp,
			     size_type i,
			     size_type j) const;

	/** 
	 * \brief Probabilty of base pair in a specified loop
	 *
	 *`single sequence folding-specific code
	 * 
	 * @param ip left end of inner base pair
	 * @param jp right end of inner base pair
	 * @param i left end of loop enclosing base pair
	 * @param j right end of loop enclosing base pair
	 * 
	 * @return probability that (ip,jp) is inner base pair in the loop closed by i and j
	 *
	 * @see prob_basepair_in_loop()
	 *
	 * @note pre: loop probs available, alifold not used
	 */
	double
	arc_in_loop_prob_noali(size_type ip,
			       size_type jp,
			       size_type i,
			       size_type j) const;

	/** 
	 * \brief Computes the Qm2 matrix
	 *
	 * The method creates and fills the Qm2 matrix needed for
	 * prob_unpaired_in_loop().
	 * 
	 * @pre McCaskill matrices are computed and accessible.
	 */
	void
	compute_Qm2();

	/** 
	 * \brief Computes the Qm2 matrix (alifold)
	 *
	 * The method creates and fills the Qm2 matrix needed for
	 * prob_unpaired_in_loop() if alifold is used.
	 * 
	 * @pre McCaskill alifold matrices are computed and accessible.
	 */
	void
	compute_Qm2_ali();

	/** 
	 * \brief Computes the McCaskill matrices and keeps them accessible
	 * 
	 * Allocates and fills the McCaskill matrices.
	 *
	 * @pre sequence_ has exactly one row
	 *
	 * @param params parameters for partition folding
	 * @param inLoopProbs whether to compute information for in loop probablities
	 * @param local_copy whether to make a local copy of ViennaRNA's data structures
	 * 
	 * @note requires linking to librna
	 */
	void
	compute_McCaskill_matrices(const PFoldParams &params, bool inLoopProbs, bool local_copy=true);
	
	
	
	/** 
	 * \brief Computes the McCaskill matrices and keeps them accessible (alifold)
	 * 
	 * Allocates and fills the McCaskill alifold matrices.
	 *
	 * @param params parameters for partition folding
	 * @param inLoopProbs whether to compute and keep information for in loop probablities
	 * @param local_copy whether to make a local copy of ViennaRNA's data structures
	 * 
	 * @note requires linking to librna
	 */
	void
	compute_McCaskill_alifold_matrices(const PFoldParams &params, bool inLoopProbs, bool local_copy=true);



#endif // HAVE_LIBRNA
	
    };

} // end namespace LocARNA

#endif // LOCARNA_RNA_ENSEMBLE_IMPL_HH

