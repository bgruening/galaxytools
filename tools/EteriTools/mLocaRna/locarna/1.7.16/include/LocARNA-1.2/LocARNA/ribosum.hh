#ifndef LOCARNA_RIBOSUM_HH
#define LOCARNA_RIBOSUM_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>

#include "matrix.hh"
#include "alphabet.hh"

namespace LocARNA {

    //! @brief Represents ribosum similarity matrices
    //!
    //! Reads ribosum data from file and provides access.
    //! 
    class Ribosum {
    public:
	typedef Matrix<double> matrix_t; //!< type of a matrix
    protected:
	typedef Alphabet<std::string> alphabet_type; //!< type of alphabet

	std::string name; //!< name of ribosum
	matrix_t bm;  //!< scores for base matches, 4x4 matrix
	matrix_t am;  //!< scores for basepair/arc matches,
	//!< 16x16 matrix


	alphabet_type basename_alphabet; //!< alphabet of base names
	alphabet_type arcname_alphabet; //!< alphabet of arc names

	Alphabet<char> char_basename_alphabet; //!< alphabet of base names as characters 


    protected:


	/** 
	 * Read matrix from input stream
	 * 
	 * @param in input stream
	 * @param[out] mat matrix to be read 
	 * @param names alphabet names provided as Alphabet (of strings)
	 * 
	 * @return input stream after reading
	 */
	std::istream &
	read_matrix(std::istream &in, 
		    matrix_t &mat,
		    const alphabet_type &names) const;
    
	/** 
	 * Write matrix to output stream
	 * 
	 * @param out output stream
	 * @param mat matrix to be written
	 * @param alph alphabet names provided as Alphabet (of strings)
	 * 
	 * @return output stream after writing
	 */    
	std::ostream &
	write_matrix(std::ostream &out, 
		     const matrix_t &mat,
		     const alphabet_type &alph) const;
    
	//! @brief Construct empty
	Ribosum();
	
	//! reads the standard ribosum file format
	//! @param in input stream
	void
	read_ribosum(std::istream &in);

    
    protected:

	//! transform the basename alphabet to alphabet over characters 
	Alphabet<char> make_char_alphabet() const;

	/** 
	 * Set alphabet of base names
	 * 
	 * @param a array of strings of the base names
	 * @post basename_alphabet and char_basename_alphabet are initialized
	 */
	void set_basename_alphabet(const std::string a[]) {
	    basename_alphabet=alphabet_type(std::vector<std::string>(&a[0],&a[4]));
	    char_basename_alphabet = make_char_alphabet();
	}
    
	/** 
	 * Set alphabet of base names
	 * 
	 * @param a array of strings of the arc names
	 * @post arcname_alphabet is initialized
	 */
	void set_arcname_alphabet(const std::string a[]) {
	    arcname_alphabet=alphabet_type(std::vector<std::string>(&a[0],&a[16]));
	}

    public:
	/** 
	 * Construct from file
	 * 
	 * @param filename name of the input file
	 */
	Ribosum(const std::string &filename);
    
	//! @brief virtual destructor
	virtual
	~Ribosum();
    
	/** 
	 * Get base match scores
	 * @return the matrix of base match scores
	 */
	const matrix_t &get_basematch_scores() const {return bm;}
    
	/** 
	 * Get arc match scores
	 * @return the matrix of arc match scores
	 */
	const matrix_t &get_arcmatch_scores() const {return am;}

	//! Get the basename alphabet as alphabet over strings
	//! @return alphabet of strings
	const alphabet_type &string_alphabet() const {return basename_alphabet;}
    
	//! Get the basename alphabet as alphabet over characters
	//! @return basename alphabet
	const Alphabet<char> &alphabet() const {return char_basename_alphabet;}
    
	//! Get name of ribosum
	//! @return name of ribosum
	const std::string & get_name() const {return name;}
    
	//! \brief Get base match score
	//! 
	//! @param i character of first nucleotide
	//! @param j character of second nucleotide
	//! @return ribosum score for matching nucleotides i and j.
	double basematch_score(char i,char j) const {
	    return bm(alphabet().idx(i),alphabet().idx(j));
	}

	//! \brief Get arc match score
	//! 
	//! @param i left character of first arc
	//! @param j right character of first arc
	//! @param k left character of second arc
	//! @param l right character of second arc
	//! @return ribosum score for matching an arc of nucleotides i and j 
	//! with an arc of nucleotides k and l.
	double arcmatch_score(char i,char j,char k,char l) const {
	    return am(alphabet().idx(i)*4+alphabet().idx(j), alphabet().idx(k)*4+alphabet().idx(l));
	}
    
	friend std::ostream & operator << (std::ostream &out, const Ribosum &ribosum);
	
    };


    //! @brief Represents ribosum similarity matrices including raw frequencies
    //!
    //! Extension of the ribosum class that maintains additional "raw"
    //! information in particular frequencies of bases, basematches,
    //! basepairs, and arcmatches
    //!
    class RibosumFreq : public Ribosum {
    public:

	/** 
	 * Construct from file
	 * 
	 * @param filename Name of input file
	 * 
	 * @note The file has to be in an extended ribosum file format
	 * including frequencies.
	 *
	 */
	RibosumFreq(const std::string &filename);

    protected:
	RibosumFreq();
    
	matrix_t base_probs_; //!< matrix of base probabilities
	matrix_t base_nonstruct_probs_; //!< matrix of base probabilities in non-structural context
	matrix_t basepair_probs_; //!< matrix of base pair probabilities
	matrix_t basematch_probs_; //!< matrix of base match probabilties 
	matrix_t arcmatch_probs_; //!< matrix of arc match probabilities

    public:

	//! Get probability of a base
	//! @param i nucleotide character
	//! @return probability of nucleotide character
	double
	base_prob(char i) const {
	    return base_probs_(alphabet().idx(i),0);
	}

	//! Get probability of a base occuring in a non-structural match
	//! @param i nucleotide character
	//! @return probability of nucleotide 
	//! character in non-structural matches
	double
	base_nonstruct_prob(char i) const {
	    return base_nonstruct_probs_(alphabet().idx(i),0);
	}
    
	//! Get base prob matrix
	//! @return matrix of probabilities of all bases/nucleotides 
	const matrix_t &
	get_base_probs() const {
	    return base_probs_;
	}

	//! Get base nonstruct prob matrix
	//! @return matrix of probabilities of all bases/nucleotides
	//! in nonstructural matches
	const matrix_t &
	get_base_nonstruct_probs() const {
	    return base_nonstruct_probs_;
	}
    
	//! Get probability of a basepair
	//! @param i left nucleotide character
	//! @param j right nucleotide character
	//! @return probability of base pair of nucleotides i and j
	double
	basepair_prob(char i,char j) const {
	    return basepair_probs_(alphabet().idx(i),alphabet().idx(j));
	}

	//! Get base pair prob matrix
	//! @return matrix of probabilities of base pairs
	const matrix_t &
	get_basepair_probs() const {
	    return basepair_probs_;
	}
    
	//! Get the probability of a base match
	//! @param i left nucleotide character
	//! @param j right nucleotide character
	//! @return probability of match of nucleotides i and j
	double
	basematch_prob(char i,char j) const {
	    return basematch_probs_(alphabet().idx(i),alphabet().idx(j));
	}

	//! Get basematch prob matrix
	//! @return matrix of probabilities of base matches
	const matrix_t &
	get_basematch_probs() const {
	    return basematch_probs_;
	}
    
	//! Get the probability of an arcmatch
	//! @param i left nucleotide character of first base pair
	//! @param j right nucleotide character of first base pair
	//! @param k left nucleotide character of second base pair
	//! @param l right nucleotide character of second base pair
	//! @return probability of matching a 
	//! @return probability for matching an arc of nucleotides i and j 
	//! with an arc of nucleotides k and l.
	double
	arcmatch_prob(char i, char j, char k, char l) const {
	    return arcmatch_probs_(alphabet().idx(i)*4+alphabet().idx(j), alphabet().idx(k)*4+alphabet().idx(l));
	}

	//! Get arcmatch prob matrix
	//! @ return matrix of arc match probabilities
	const matrix_t &
	get_arcmatch_probs() const {
	    return arcmatch_probs_;
	}
    
    
	//! \brief probability that a nucleotide/base occurs unpaired
	//! @param i nucleotide character
	//! @return probability that nucleotide is unpaired
	double
	base_unpaired_prob(char i) const;
    

	//! \brief Get corrected score for a base match
	//!
	//! @param i first nucleotide character
	//! @param j second nucleotide character
	//! @return corrected score that i and j are matched
	//! 
	//! @note the score is computed as log odd of frequencies
	//! for seeing i and j matched without incident structure divided by 
	//! the background to see i and j without incident structure.
	//! 
	//! @note The arguments are characters of the alphabet/nucleotides.
	//! @note Currently, not tabellized. Thus, we have some computational overhead.
	double
	basematch_score_corrected(char i,char j) const;

    
	//! \brief Print the corrected score of base matches
	//! @see basematch_score_corrected()
	void
	print_basematch_scores_corrected() const;
    
	/** 
	 * Read one matrix in ribosum input 
	 * 
	 * @param in input stream 
	 * @param header header string
	 * @param[out] mat matrix object to be read into 
	 * @param xdim first dimension of matrix
	 * @param ydim second dimension of matrix
	 */
	void
	read_matrix(std::istream &in, const std::string &header, matrix_t &mat, size_t xdim, size_t ydim);

	//! Write the ribosum matrices as C++ code to cout
	//! @param ribname Name of ribosum
	//! @note this is used to precompile ribosum matrices, 
	//! such that their parameters can be linked to the code 
	//! and we do not have to depend on data files.
	void
	write_CC_code(const std::string &ribname) const;

     
	/** 
	 * Write single matrix to output stream
	 * 
	 * @param out output stream
	 * @param name name of matrix
	 * @param mat matrix
	 * 
	 * @return output stream after writing
	 */
	std::ostream &
	write_matrix(std::ostream &out, const std::string &name, const Matrix<double> &mat) const;
    
	friend std::ostream & operator << (std::ostream &out, const RibosumFreq &ribosum);

    private:
    
	void write_CC_matrix(const std::string &ribname,const std::string &matname,
			     int x, int y, const Ribosum::matrix_t &m) const;

	void
	read_frequencies(std::istream &in);
    };

} // end namespace LocARNA

#endif //LOCARNA_RIBOSUM_HH
