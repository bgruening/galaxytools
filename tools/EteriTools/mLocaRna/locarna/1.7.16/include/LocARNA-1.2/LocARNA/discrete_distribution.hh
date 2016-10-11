#ifndef LOCARNA_DISCRETE_DISTRIBUTION
#define LOCARNA_DISCRETE_DISTRIBUTION

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

namespace LocARNA {
    
    // ------------------------------------------------------------
    /** Generate discrete distributions from pseudo random numbers
     */
    class DiscreteDistribution {
    private:
	std::vector<double> distrib_acc_;
    public:
	/** 
	 * \brief Construct with distribution vector
	 *  
	 * @param distvec vector defining discrete distribution 
	 */
	DiscreteDistribution(const std::vector<double> &distvec)
	{
	    if (distvec.size()==0) return;
	
	    distrib_acc_.resize(distvec.size());

	    distrib_acc_[0] = distvec[0];
	    for (size_t i=1; i<distvec.size(); i++) {
		distrib_acc_[i]+=distrib_acc_[i-1]+distvec[i];
	    }
	    double max_acc = distrib_acc_[distvec.size()-1];
	    for (size_t i=0; i<distvec.size(); i++) {
		distrib_acc_[i]/=max_acc;
	    }
	}
    
	/** 
	 * \brief get random number
	 * @param x random number
	 * 
	 * @return random number from discrete distribution
	 */
	size_t
	operator () (size_t x) {
	    double y=x/(double)RAND_MAX;
	
	    size_t i=0;
	    for (; i<distrib_acc_.size() && y>distrib_acc_[i]; ++i)
		;
	
	    assert(i<=distrib_acc_.size());
	    return i;
	}
    };

} // end namespace LocARNA

#endif //LOCARNA_DISCRETE_DISTRIBUTION

