#ifndef LOCARNA_PFOLD_PARAMS_HH
#define LOCARNA_PFOLD_PARAMS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

namespace LocARNA {

    /**
     * \brief Parameters for partition folding
     *
     * Describes certain parameters for the partition folding of 
     * a sequence or alignment.
     *
     * @see RnaEnsemble
     *
    */
    class PFoldParams {
	bool noLP_;
	bool stacking_;
    public:
	/** 
	 * Construct with all parameters
	 * 
	 * @param noLP
	 * @param stacking 
	 */
	PFoldParams(bool noLP,
		    bool stacking
		    )
	    : noLP_(noLP),
	      stacking_(stacking) 
	{}
	
	/** 
	 * @brief Check no LP flag
	 * 
	 * @return value of flag 
	 */
	bool noLP() const {return noLP_;}
	
	/** 
	 * @brief Check stacking flag
	 * 
	 * @return value of flag 
	 */
	bool stacking() const {return stacking_;}
    };


}

#endif // LOCARNA_PFOLD_PARAMS_HH
