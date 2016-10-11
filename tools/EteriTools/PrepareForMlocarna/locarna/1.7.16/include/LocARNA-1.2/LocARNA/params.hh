#ifndef LOCARNA_PARAMS_HH
#define LOCARNA_PARAMS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <vector>
#include <string>

namespace LocARNA {

    class AnchorConstraints;
    class TraceController;

    /**
       \brief Description of free end gaps.
       
       Decodes the description given by a string of 4 characters '+'/'-'
       and provides methods with reasonable names.
    */
    class FreeEndgapsDescription {
	std::vector<bool> desc;
    public:
    
	/** 
	 * @brief Construct from string description
	 * 
	 * @param d description given by a string of 4 characters '+'/'-'
	 *
	 * @note the string description is suited to specify free end gaps in this way on the command line
	 */
	FreeEndgapsDescription(const std::string &d)
	    : desc(4)
	{
	    if (d.length()>=4) {
		for (size_t i=0; i<4; i++) desc[i] = (d[i]=='+');
	    } else {
		for (size_t i=0; i<4; i++) desc[i] = false;
	    }
	}
    
	/** 
	 * Are gaps free at left end of first sequences?
	 * @return whether free end gaps are allowed
	 */
	bool
	allow_left_1() const  {
	    return desc[0];
	}

	/** 
	 * Are gaps free at right end of first sequences?
	 * @return whether free end gaps are allowed
	 */
	bool
	allow_right_1() const  {
	    return desc[1];
	}

	/** 
	 * Are gaps free at left end of second sequences?
	 * @return whether free end gaps are allowed
	 */
	bool
	allow_left_2() const  {
	    return desc[2];
	}

	/** 
	 * Are gaps free at right end of second sequences?
	 * @return whether free end gaps are allowed
	 */
	bool
	allow_right_2() const  {
	    return desc[3];
	}
    };


    /**
       \brief Parameter for alignment by Aligner
             
       Collects the parameters for the aligner object.  These parameters
       controll the kind of alignment (local/global),
       restrictions/constraints on the alignment and certain heuristics.
       Parameters for the score are collected in a different class.
       
       @see Aligner
       @see ScoringParams
    */
    class AlignerParams {  
    public:

	const bool no_lonely_pairs; //!< no lonely pairs option
  
	const bool STRUCT_LOCAL; //!< allow exclusions for maximizing alignment of connected substructures 
	const bool SEQU_LOCAL; //!< sequence local alignment / maximize alignment of subsequences

	const FreeEndgapsDescription free_endgaps; //!< description of potentially allowed free end gaps
    
	const bool DO_TRACE; //!< whether do perfom trace back

	const TraceController &trace_controller; //!< trace controller controlling allowed trace cells
    
	const int max_diff_am; //!< maximal difference of arc lengths in arc match
	const double min_am_prob; //!< minimal probability of an arc match
	const double min_bm_prob;  //!< minimal probability of a base match
    
	const bool stacking; //!< whether to use stacking

	const AnchorConstraints &constraints; //!< anchor constraints
    
    
	/** 
	 * Construct with parameters
	 * 
	 * @param _no_lonely_pairs no lonely pairs option
	 * @param _STRUCT_LOCAL allow exclusions for maximizing alignment of connected substructures 
	 * @param _SEQU_LOCAL  sequence local alignment / maximize alignment of subsequences
	 * @param _free_endgaps  description of potentially allowed free end gaps
	 * @param _trace_controller  trace controller controlling allowed trace cells
	 * @param _max_diff_am maximal difference of arc lengths in arc match
	 * @param _min_am_prob minimal probability of an arc match
	 * @param _min_bm_prob  minimal probability of a base match
	 * @param _stacking whether to use stacking
	 * @param _constraints  anchor constraints
	 */
	AlignerParams(bool _no_lonely_pairs, 
		      bool _STRUCT_LOCAL, 
		      bool _SEQU_LOCAL, 
		      std::string _free_endgaps,
		      const TraceController &_trace_controller,
		      int _max_diff_am,
		      const double _min_am_prob,
		      const double _min_bm_prob,
		      bool _stacking,
		      const AnchorConstraints &_constraints		  
		      ):
	    no_lonely_pairs(_no_lonely_pairs),
	    STRUCT_LOCAL(_STRUCT_LOCAL),
	    SEQU_LOCAL(_SEQU_LOCAL),
	    free_endgaps(_free_endgaps),
	    DO_TRACE(true),
	    trace_controller(_trace_controller),
	    max_diff_am(_max_diff_am),
	    min_am_prob(_min_am_prob), 
	    min_bm_prob(_min_bm_prob),	   
	    stacking(_stacking),
	    constraints(_constraints) 
	{}
  
    };

}

#endif // LOCARNA_PARAMS_HH
