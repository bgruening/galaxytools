#ifndef LOCARNA_MATCH_PROBS
#define LOCARNA_MATCH_PROBS

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string>

#include "matrix.hh"

namespace LocARNA {

    class StralScore;
    class RnaData;
    template <class T> class Alphabet;
    class Sequence;

    /**
     * \brief Provides probabilities for each match.
     *
     * The probabilities are either computed or
     * read in from file
     * 
     * for computing probabilities, we offer two methods:
     *
     * 1)
     * For computing the probabilities the class 
     * uses a pairHMM analogously to PROBCONS.
     * Also, the class reads transition probabilities
     * from a file in the format of Probcons, which
     * allows to use their parameter files
     *
     * 2)
     * similar to proba/probalign,  use statistical-mechanics-like
     * model. Assume alignments are Boltzman distributed,
     * calc match probs via partition function
     *
     * the second approach supports Stral-like scoring (using pf_struct_weight as "alpha")
     *
     */
    class MatchProbs {
    public:    
	typedef size_t size_type; //!< size
	
	//! construct as empty object
	MatchProbs();

	//! construct from file
	MatchProbs(const std::string &filename);
    
	/**
	 * read probcons parameter file
	 * and compute match probabilities
	 * for the two given sequences
	 * @throws failure
	 */
	void
	pairHMM_probs(const Sequence &seqA,
		      const Sequence &seqB,
		      const std::string & file);
    
	/**
	 * calculate edge probabilities via statistical mechanics model (partition function)
	 * get match scores via matrix sim_mat.
	 * The method accepts the matrix sim_mat together with an alphabet.
	 * The alphabet object is necessary to translate sequence symbols
	 * to the indices in this matrix.
	 */
	void
	pf_probs(const RnaData &rnaA,
		 const RnaData &rnaB,
		 const Matrix<double> &sim_mat,
		 const Alphabet<char> &alphabet,
		 double gap_opening,
		 double gap_extension,
		 double pf_struct_weight,
		 double temp,
		 bool flag_local);
    
	/**
	 * read the probabilities from a stream
	 * assumes matrix starting 0,0 whereas sequences start 1,1
	 */
	std::istream &
	read(std::istream &in);

	//! read the probabilities from a file
	void
	read(const std::string &filename);


	/**
	 * read the probabilities from a stream
	 * read 'sparse' format "i j p"
	 */
	std::istream &
	read_sparse(std::istream &in, size_type lenA, size_type lenB);

	/**
	 * read the probabilities from a file
	 * read 'sparse' format "i j p"
	 */
	void
	read_sparse(const std::string &filename, size_type lenA, size_type lenB);
    
	/**
	 * write the probabilities to a stream
	 * writes matrix starting 0,0 whereas sequences start 1,1
	 */
	std::ostream &
	write(std::ostream &out) const;
    
	//! write the probabilities to a file
	void
	write(const std::string &filename) const;

	/**
	 * write the probabilities to a stream, only probs >= threshold
	 * use format "i j p"
	 */
	std::ostream &
	write_sparse(std::ostream &out, double threshold) const;
    
	/**
	 * write the probabilities to a file, only probs >= threshold
	 * use format "i j p"
	 */
	void
	write_sparse(const std::string &filename, double threshold) const;

	//! get the length of the first sequence
	size_type get_lenA() const {return probs.sizes().first;}
    
	//! get the length of the second sequence
	size_type get_lenB() const {return probs.sizes().second;}
    
	//! return the match probability for the two bases
	double prob(size_t i, size_t j) const {
	    assert(1<=i && i<probs.sizes().first);
	    assert(1<=j && j<probs.sizes().second);

	    return probs(i,j);
	}
    
    private:
	//! perform the partition version of Gotoh's algorithm
	void
	pf_gotoh(size_type lenA,
		 size_type lenB,
		 Matrix<double> &zM,
		 Matrix<double> &zA,
		 Matrix<double> &zB,
	     
		 const StralScore &score,

		 double temp,
	     
		 bool local
		 );
    
	Matrix<double> probs; //!< the base match probabilities
    
	/**
	 * @brief Maintains parameter for computing match probabilities 
	 *
	 * The class is used in sequence alignment a la Probcons by MatchProbs::pairHMM_probs()
	 *
	 * Supports reading parameter from file
	 *
	 * @see MatchProbs 
	 */
	class ProbConsParameter {
	public:
	    // ------------------------------------------------------------
	    // transition probabilities
	    // there are three states M, X, Y (plus implicitely start and end state)
	    double initM; //!< transition probability initM
	    double initX; //!< transition probability initX 
	    double initY; //!< transition probability initY 
	    double startX; //!< transition probability startX
	    double startY; //!< transition probability startY
	    double extendM; //!< transition probability extendM
	    double extendX; //!< transition probability extendX
	    double extendY; //!< transition probability extendY
	    double startMFromX; //!< transition probability startMFromX
	    double startMFromY; //!< transition probability startMFromY

	    std::string basenames; //!< base names
    
	    Matrix<double> emmission; //!< matrix of emmission probabilities
	    std::vector<double> background; //!< vector of background probabilities
    
	    /** 
	     * Construct from file
	     * 
	     * @param filename
	     *
	     * @throws failure
	     */
	    ProbConsParameter(const std::string &filename);
	};

    };
    
}

#endif // LOCARNA_MATCH_PROBS
