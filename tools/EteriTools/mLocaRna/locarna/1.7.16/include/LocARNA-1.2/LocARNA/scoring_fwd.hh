#ifndef LOCARNA_SCORING_FWD_HH
#define LOCARNA_SCORING_FWD_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "infty_int.hh"

namespace LocARNA {

    //! type of the locarna score as defined by the class Scoring
    typedef long int score_t;
    
    //! an extended score_t that can store and calculate with
    //! infinite values (i.p. we use -infty for invalid matrix entries)
    typedef InftyInt infty_score_t;

    typedef TaintedInftyInt tainted_infty_score_t;


    //! type of partition functions
#ifdef VERY_LARGE_PF
    typedef long double pf_score_t;
#else
    typedef double pf_score_t;
#endif

    class Scoring;
    
} // end namespace LocARNA

#endif // LOCARNA_SCORING_FWD_HH


