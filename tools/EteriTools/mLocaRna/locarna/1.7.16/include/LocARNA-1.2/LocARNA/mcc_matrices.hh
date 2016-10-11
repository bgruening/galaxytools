#ifndef LOCARNA_MCC_MATRICES_HH
#define LOCARNA_MCC_MATRICES_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <assert.h>

#define PUBLIC // for Vienna
extern "C" {
#include <ViennaRNA/params.h> // import pf_paramT definition
}

namespace LocARNA {

    class McC_matrices_base {
    protected:
    	size_t length_;     //!< sequence length
	bool local_copy_; //!< whether pointers point to local copies of data structures
	
	FLT_OR_DBL *qb_;  //!< Q<sup>B</sup> matrix					
	FLT_OR_DBL *qm_;  //!< Q<sup>M</sup> matrix					
	
	FLT_OR_DBL *bppm_;  //!< base pair probability matrix
	
    	int* iindx_;        //!< iindx from librna's get_iindx()
    	/** 
	 * @brief construct empty
	 */
	McC_matrices_base();

	/** 
	 * @brief initialize
	 * 
	 * @param length length of sequence 
	 */
	void
	init(size_t length);

    public:
	FLT_OR_DBL *q1k_; //!< 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)	
	FLT_OR_DBL *qln_; //!< 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)      
	
	pf_paramT *pf_params_; //!< parameters for pf folding
	
		
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual
	~McC_matrices_base();

	
	//! \brief index in triagonal matrix
	size_t iidx(size_t i,size_t j) const {
	    assert(1<=i);
	    assert(i<=j);
	    assert(j<=length_);

	    return iindx_[i]-j;
	}

	/** 
	 * @brief Read access matrix bppm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL bppm(size_t i, size_t j) const { return bppm_[iidx(i,j)]; }

	/** 
	 * @brief Read access matrix qb
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL qb(size_t i, size_t j) const { return qb_[iidx(i,j)]; }

	/** 
	 * @brief Read access matrix qm
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	FLT_OR_DBL qm(size_t i, size_t j) const { return qm_[iidx(i,j)]; }
	    
    protected:
	
	/** 
	 * @brief free all local copies of data structures
	 */
	void free_all_local();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_base &McCmat);
    };
    
    //! @brief  structure for McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_pf_arrays() and get_bppm() of Vienna librna
    class McC_matrices_t : public McC_matrices_base {
	char *ptype_;	   //!< pair type matrix					
	
    public:

	char *sequence_;  //!< 0-terminated sequence string
	short *S_;        //!< 'S' array (integer representation of nucleotides)	
	short *S1_;	   //!< 'S1' array (2nd integer representation of nucleotides)	
	
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param sequence the sequence as 0-terminated C-string 
	 * @param local_copy  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_matrices_t(char *sequence, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual 
	~McC_matrices_t();


	/** 
	 * @brief Access matrix ptype
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	char ptype(size_t i, size_t j) const { return ptype_[iidx(i,j)]; }

	/** 
	 * @brief Reverse ptype
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	char 
	rev_ptype(size_t i, size_t j) const;

    protected:

	/** 
	 * Free all data structures of the Vienna package
	 */
	void free_all();

	//! \brief deep copy all data structures 
	void
	deep_copy(const McC_matrices_t &McCmat);
    };

     //! @brief  structure for Alifold-McCaskill matrices pointers
    //!
    //! Contains pointers to matrices made accessible through
    //! get_alipf_arrays() and alipf_export_bppm() of Vienna librna
    class McC_ali_matrices_t : public McC_matrices_base {
    protected:
	size_t n_seq_;     //!< sequence length
	
    public:

	short **S_;       //!< 'S' array (integer representation of nucleotides)	
	short **S5_;	   //!< 'S5' array
	short **S3_;	   //!< 'S3' array
	unsigned short  **a2s_;  //!< 'a2s' array
	char **Ss_;	   //!< 'Ss' array
	
    protected:
	short *pscore_; //!< alifold covariance/conservation scores
    public:
	/** 
	 * @brief construct by call to VRNA lib functions and optionally make local copy
	 * 
	 * @param n_seq number of sequenes in alignment
	 * @param length length of sequences in alignment 
	 * @param local_copy  if TRUE, copy the data structures; otherwise, only store pointers
	 */
	McC_ali_matrices_t(size_t n_seq, size_t length, bool local_copy);
	
	/** 
	 * @brief destruct, optionally free local copy
	 */
	virtual
	~McC_ali_matrices_t();


	/** 
	 * @brief Access matrix pscore
	 * 
	 * @param i first index
	 * @param j second index
	 * 
	 * @return matrix entry 
	 */
	short pscore(size_t i, size_t j) const { return pscore_[iidx(i,j)]; }


    protected:

	/** 
	 * @brief Free McCaskill/VRNA data structures
	 */
	void free_all();

	/**
	 * @brief deep copy all data structures
	 * @param McCmat object to copy
	 */
	void
	deep_copy(const McC_ali_matrices_t &McCmat);
    };
    
} // end namespace LocARNA


#endif // LOCARNA_MCC_MATRICES_HH
