#ifndef LOCARNA_CONFUSION_MATRIX
#define LOCARNA_CONFUSION_MATRIX

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstddef>
#include <cassert>
#include <string>

namespace LocARNA {
    class RnaStructure;

    /**
     * @brief Compare RNA secondary structure by their confusion matrix
     *
     * Computes confusion matrix and Matthews' correlation coefficient
     * @todo implement slide rule
     *
     * Base pair filters provide a generic way to count only a subset
     * of potential base pairs, e.g. only canonical base pairs, or
     * base pairs in a certain range l..u
     */
    class ConfusionMatrix {
    public:
	typedef std::pair<size_t,size_t> bp_t;
	
	/**
	 * @brief basic class for base pair filters (no filtering)
	 */
	class BPFilter {
	public:
	    BPFilter()
	    {};

	    virtual
	    ~BPFilter()
	    {};
	    
	    virtual
	    bool
	    operator () (size_t i, size_t j) const {
		assert(1<=i); // sequences are 1-based
		assert(i<=j);
		return true;
	    }
	    
	    virtual 
	    bool operator () (const bp_t &bp) const {
		return (*this)(bp.first,bp.second);
	    }
	};

	/**
	 * @brief loop size base pair filter
	 */
	class BPMinLoopSizeFilter : public BPFilter {
	    size_t mls_;
	public:
	    BPMinLoopSizeFilter(size_t mls)
	    : BPFilter(),
	      mls_(mls)
	    {}

	    ~BPMinLoopSizeFilter() {};
	    
	    bool
	    operator () (size_t i, size_t j) const {
		assert(i>=1);
		return i+mls_<j;
	    }
	};
	
	/**
	 * @brief base pair range filter 
	 */
	class SpanRangeBPFilter: public BPFilter {
	    size_t lo_;
	    size_t up_;
	public:
	    /** 
	     * @brief Construct with range
	     * 
	     * @param lo lower bound
	     * @param up upper bound, if u>0; else if up==0, no upper bound
	     */
	    SpanRangeBPFilter(size_t lo, size_t up)
		: BPFilter(), lo_(lo), up_(up)
	    {
		assert(1<=lo);
	    }

	    ~SpanRangeBPFilter() {};
	    
	    bool
	    operator () (size_t i, size_t j) const {
		assert(i<=j);
		return lo_<=(j-i) && (up_==0 || (j-i)<=up_);
	    }
	};

	/**
	 * @brief base pair filter to allow only canonical base pairs 
	 */
	class CanonicalBPFilter: public BPFilter {
	    const std::string &sequence_;
	public:
	    CanonicalBPFilter(const std::string &sequence)
		: BPFilter(), 
		  sequence_(sequence)
	    {}
	    
	    ~CanonicalBPFilter() {};
	    
	    bool
	    operator () (size_t i, size_t j) const {
		return BPFilter::operator ()(i,j) && canonical(sequence_[i],sequence_[j]);
	    }
	private:
	    static
	    bool
	    canonical(char x, char y) {
		static std::string cpairs="AUCGGUUAGCUG";
		for(size_t i=0; i<cpairs.length(); i+=2) {
		    if (x==cpairs[i] && y==cpairs[i+1]) return true;
		}
		return false;
	    }
	};

        class CombinedBPFilter: public BPFilter {
	    const BPFilter &fa_;
	    const BPFilter &fb_;
	public:
	    CombinedBPFilter(const BPFilter &fa, const BPFilter &fb): 
		BPFilter(),fa_(fa),fb_(fb)
	    {}
	    
	    ~CombinedBPFilter()
	    {}
	    
	    bool
	    operator () (size_t i, size_t j) const {
		return fa_(i,j) && fb_(i,j);
	    }
	};
	
	
	/** 
	 * Construct with reference and predicted structure (given as dot-bracket strings) 
	 * 
	 * @param ref  reference structure
	 * @param pred predicted structure
	 */
	ConfusionMatrix(const std::string &ref,
			const std::string &pred,
			const BPFilter &filter=BPFilter()
			);
	
	/** 
	 * Construct with reference and predicted structure (given as base pair sets) 
	 * 
	 * @param ref  reference structure
	 * @param pred predicted structure
	 */	
	ConfusionMatrix(const RnaStructure &ref,
			const RnaStructure &pred,
			const BPFilter &filter=BPFilter()
			);

	/** 
	 * True positives
	 * 
	 * @return number of true positive base pairs 

	 * @note A true positive is a pair which is basepaired in
	 * both prediction and annotation (slide not implemented yet).
	 */
	size_t
	tp() const { return tp_; }
    
	/** 
	 * True negatives
	 * 
	 * @return number of true negative base pairs 
	 *
	 * @note A true negative is a pair which is neither annotated
	 * nor predicted to basepair
	 */
	size_t
	tn() const { return tn_; }

	/** 
	 * False positives
	 * 
	 * @return number of false positive base pairs
	 *
	 * @note A false positive is a pair predicted to basepair but
	 * not annotated.
	 */
	size_t
	fp() const {	
	    return fp_;
	}

	/** 
	 * False positives in conflict with reference
	 * 
	 * @return number of false positive base pairs, where only conflicting base pairs
	 * are counted as false
	 *
	 * @note A false positive in conflict is a predicted pair that
	 * conflicts with reference
	 */
	size_t
	cfp() const {	
	    return cfp_;
	}

	/** 
	 * False negatives
	 * 
	 * @return number of false negative base pairs 

	 * @note A false negative is pair annotated to basepair but
	 * not predicted.
	 */
	size_t
	fn() const {	
	    return fn_;
	}

	/** 
	 * Positive prediction value
	 *
	 * aka precision
	 * 
	 * @return PPV = TP/(TP+FP)
	 */
	double
	ppv() const;

	/** 
	 * Sensitivity
	 * 
	 * aka recall
	 *
	 * @return SENS = TP/(TP+FN)
	 */
	double
	sens() const;

	/** 
	 * Specificity
	 * 
	 * @return SPEC = TN/(TN+FP)
	 */
	double
	spec() const;

 	/** 
	 * F1 score (aka F-score, F-measure)
	 * 
	 * harmonic mean of PPV and SENS
	 *
	 * @return F1 = PPV*SENS / (PPV+SENS), if PPV+SENS!=0; 0, otherwise 
	 */
	double
	f1_score() const;

 	/** 
	 * F1 score (aka F-score, F-measure), considering only conflicting false positives
	 * 
	 * harmonic mean of CPPV and SENS
	 *
	 * @return F1 = CPPV*SENS / (CPPV+SENS), if CPPV+SENS!=0; 0, otherwise 
	 */
	double
	cf1_score() const;

	/** 
	 * Matthews' correlation coefficient
	 * 
	 * @return MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) 
	 */
	double
	mcc() const;

    protected:
	
	/** 
	 * @brief Count common base pairs
	 * 
	 * @param s1   first structure 
	 * @param s2  second structure
	 * 
	 * @return number of base pairs common to both structures  
	 */
	static
	size_t
	count_common_bps(const RnaStructure &s1,
			 const RnaStructure &s2,
			 const BPFilter &filter);

	/** 
	 * @brief Count conflicting base pairs
	 * 
	 * @param s1   first structure 
	 * @param s2  second structure
	 *
	 * Two base pairs (i,j) and (i',j') conflict if and
	 * only if they share exactly one end.
	 * 
	 * @return number of base pairs in s1 that conflict with s2
	 */
	static
	size_t
	count_conflicting_base_pairs(const RnaStructure &s1,
				     const RnaStructure &s2,
				     const BPFilter &filter);
	/**
	 * @brief Count potential base pairs
	 * 
	 * @param length length of RNA
	 * @return number of potential base pairs
	 */
	static
	size_t
	count_potential_base_pairs(size_t length,
				   const BPFilter &filter);
	

	/**
	 * @brief Count base pairs in a structure
	 * 
	 * @param structure RNA structure
	 * @return number of base pairs
	 */
	static
	size_t
	count_base_pairs(const RnaStructure &s,
			 const BPFilter &filter);
	
    private:
	/** 
	 * @brief Compute confusion matrix
	 * 
	 * @param ref  reference structure 
	 * @param pred predicted structure
	 *
	 * @see count_base_pairs() for sequence argument returns result
	 * in attributes
	 */
	void
	compute_confusion_matrix(const RnaStructure &ref, 
				 const RnaStructure &pred,
				 const BPFilter &filter);
	
	size_t tp_;  //!< num of true positives
	size_t tn_;  //!< num of true negatives
	size_t fp_;  //!< num of false positives
	size_t cfp_; //!< num of conflicting false positives
	size_t fn_;  //!< num of false negatives
    };

} // end namespace LocARNA


#endif // LOCARNA_CONFUSION_MATRIX
