#ifndef LOCARNA_MULTIPLE_ALIGNMENT_HH
#define LOCARNA_MULTIPLE_ALIGNMENT_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include "aux.hh"
#include "string1.hh"
#include "scoring_fwd.hh"
#include "sequence_annotation.hh"

#include <assert.h>

#include <iostream>


namespace LocARNA {

    class Alignment;
    class AlignmentEdges;
    template<class T> class Alphabet;
    class BasePairs;
    class Scoring;
    class Sequence;
    
/**
 * @brief Represents a multiple alignment
 *
 * The multiple alignment is implemented as vector of name/sequence
 * pairs.
 *
 * Supports traversal of name/sequence pairs. The sequence entries support
 * mapping from columns to positions and back.
 *
 * Names are unique in a multiple alignment object.
 *
 * Sequences positions and column indices are 1..len.
 *
 * MultipleAlignment can have anchor and structure annotation and can
 * read and write them.
 *
 * @note this class is agnostic of the type of sequences in the
 * alignment; it does not check for 'allowed characters' nor transform
 * characters.  However, normalize_rna_bases() is provided to perform
 * a normalization in the case of RNAs, which is generally assumed by
 * the alignment engines and Vienna folding routines. 
 *
 * @todo because this class does not know whether it contains RNA, it
 * would be useful to have a derived class RNAMultipleAlignment. This
 * class could guarantee that its sequences are normalized RNA
 * sequences. Consequently, we could enforce by the type system that
 * RnaEnsemble is generated only from RNAMultipleAlignment etc.
 */
class MultipleAlignment {

public:
    typedef size_t size_type; //!< size type

    /**
     * @brief file format type for multiple alignments
     */
    struct FormatType {
	//! inner type
	enum type {
	    CLUSTAL, //!< (extended) clustal file format
	    FASTA  //!< fasta file format
	};
    };
	

    //! @brief type of sequence annotation.
    //! enumerates legal annotation types
    struct AnnoType {
	//! inner type
	enum type {
	    structure,       //!< structure annotation (often, constraint)
	    fixed_structure, //!< structure annotation (a single structure)
	    anchors          //!< anchor annotation (for anchor constraints)
	};
    };
    
private:
    //! prefix strings for annotations (shall be prefix unique)
    //! (no one outside of MultipleAlignment should have to know about this!)
    static const std::vector<std::string> annotation_tags;
public:
    //! @brief number of annotation types
    //! @return number of annotation types
    static
    size_t
    num_of_annotypes() {
	return annotation_tags.size();
    }
    
    /**
     * @brief A row in a multiple alignment
     * 
     * pair of a name string and a sequence string; support
     * projections 
     *
     * @see MultipleAlignment
     */
    class SeqEntry {
    public:
	typedef MultipleAlignment::size_type size_type; //!< size type
	
	typedef std::pair<pos_type,pos_type> pos_pair_t; //!< pair of positions
	
    private:	
	std::string name_; //!< name of the sequence
	std::string description_; //!< optional sequence description
	string1 seq_; //<! alignment string of the sequence 
	
    public:

	/** 
	 * @brief Construct from strings name and seq
	 * 
	 * @param name Sequence name
	 * @param seq  Sequence string
	 * @note empty description
	 */
	SeqEntry(const std::string &name,
		 const std::string &seq)
	    : name_(name), description_(""), seq_((string1)seq)
	{}
	
	/** 
	 * @brief Construct from strings name and 1-based string seq
	 * 
	 * @param name Sequence name
	 * @param seq  Sequence string
	 * @note empty description
	 */
	SeqEntry(const std::string &name, const string1 &seq)
	    : name_(name), description_(""), seq_(seq)
	{}
	
	/** 
	 * @brief Construct from strings name, description and seq
	 * 
	 * @param name Sequence name
	 * @param description Sequence description
	 * @param seq  Sequence string
	 */
	SeqEntry(const std::string &name, 
		 const std::string &description,
		 const std::string &seq)
	    : name_(name), description_(description), seq_((string1)seq)
	{}
	
	/** 
	 * @brief Construct from strings name, description and 1-based string seq
	 * 
	 * @param name Sequence name
	 * @param description Sequence description
	 * @param seq  Sequence string
	 */
	SeqEntry(const std::string &name,
		 const std::string &description,
		 const string1 &seq)
	    : name_(name), description_(description), seq_(seq)
	{}
	
	// access
	
	//! @brief (read-only) access to name
	const std::string &
	name() const {return name_;}
	
	//! @brief (read-only) access to description
	const std::string &
	description() const {return description_;}

	//! @brief (read-only) access to seq
	const string1 &
	seq() const {return seq_;}

	//! @brief length without gaps
	size_type
	length_wogaps() const;
	
	//****************************************
	// projections
	
	/**
	 * @brief map sequence position -> alignment column.
	 * @note time O(len)
	 * @param pos position in sequence (without gaps)
	 * as marginal cases: pos 0 maps to 0 and
	 * a too large position maps to length+1
	*/
	pos_type
	pos_to_col(pos_type pos) const;
	
	/**
	 * map alignment column -> sequence positions
	 * @note time O(len)
	 * @param col column index in aligmnent
	 * @returns pair of positions (pos1,pos2)
	 *   if column col contains a non-gap, then pos1=pos2 is the position of the gap
	 *   if column col contains a gap, then pos1 is the sequence position left of the gap or 0 and pos2 the position right of the gap or sequence length+1
	*/
	pos_pair_t
	col_to_pos(pos_type col) const;

	/** 
	 * @brief reverse sequence
	 * 
	 */
	void
	reverse() {
	    seq_.reverse();
	}
	
	/** 
	 * @brief append character to sequence
	 * @param c character
	 */
	void
	push_back(char c) {
	    seq_.push_back(c);
	}

	//! @brief write access to seq
	void
	set_seq(const string1 &seq) {seq_=seq;}


    };

    /**
     * @brief read only proxy class representing a column of the alignment 
     *
     * Allow read only access to the symbols in the column by their row index
     */
    class AliColumn {
	const MultipleAlignment &ma_;
	size_type col_index_;
    public:
	/** 
	 * @brief Construct from multiple alignment column
	 * 
	 * @param ma multiple alignment
	 * @param col_index column index
	 */
	AliColumn(const MultipleAlignment &ma,size_type col_index): ma_(ma),col_index_(col_index) {
	    assert(1<=col_index);
	    assert(col_index<=ma.length());
	}
	
	/** 
	 * @brief element access
	 * 
	 * @param row_index 0-based index of alignment row
	 * 
	 * @return character at row in the represented column
	 */
	const char &
	operator [](size_type row_index) const {
	    return ma_.seqentry(row_index).seq()[col_index_];
	}

	/** 
	 * @brief Size / Number of rows 
	 * @return number of rows of the multiple alignment
	 */
	size_type 
	size() const {
	    return ma_.num_of_rows();
	}

	/** 
	 * @brief Test equality
	 * 
	 * @param ac second alignment column
	 * 
	 * @return whether columns are equal
	 */
	bool
	operator ==(const AliColumn &ac) const {
	    bool ret = this->size()==ac.size();
	    for (size_type i=0; ret && i<size(); i++) {
		ret = ( this->ma_.seqentry(i).seq()[this->col_index_]
			== 
			ac.ma_.seqentry(i).seq()[ac.col_index_] );
	    }
	    return ret;
	}

	/** 
	 * @brief Test inequality
	 * 
	 * @param ac second alignment column
	 * 
	 * @return whether columns are equal
	 */
	bool
	operator !=(const AliColumn &ac) const {
	    return !(*this == ac);
	}

    };
    
private:
    
    //! map from string to index
    typedef std::map<std::string,size_type> str2idx_map_t;
    
    //! map annotation type to sequence annotation
    typedef std::map<size_t,SequenceAnnotation> annotation_map_t;
    
    //************************************************************
    // attributes of MultipleAlignment
    
    //! vector of alignment rows
    std::vector<SeqEntry> alig_;
    
    //! alignment/sequence annotation
    annotation_map_t annotations_;
    
    /**
     * association between names and indices, use to 
     * locate sequences by name in log time
     */
    str2idx_map_t name2idx_;
    
    // end attributes
    //************************************************************

    //! @brief create the map for translating names to indices
    void
    create_name2idx_map();

    /**
     * @brief Read alignment from input stream, expect clustalw-like format.
     *
     * @param in input stream
     * @note
     * - A header starting with CLUSTAL is ignored, but not required.
     * - Lines can be empty or of the form <name> <seq>.
     * - Names may occur multiple times. in this case seq strings <seq> are appended.
     * - The order of first occurrences of names in the stream is preserved.
     * @note overwrites/clears existing data     
     */
    void
    read_aln_clustalw(std::istream &in);

    /**
     * @brief Read alignment from input stream, expect fasta format.
     * 
     * @param in input stream
     *
     * @note Sequence descriptors have the form '>descriptor'. Any
     * white space between '>' and the name is ignored.  The sequence
     * name is the descriptor until the first blank. The rest of the
     * line is understood as sequence description.
     *
     * @note Sequences can be multiline, white space in sequences is ignored.
     * @note The order of sequences in the stream is preserved.
     * @note overwrites/clears existing data
     *
     * @todo read_aln_fasta() currently does not read anchor
     * constraints and structure. Should it? If yes, likely using
     * special fa headers >#A, >#S.
     */
    void
    read_aln_fasta(std::istream &in);
    
public:
    
    //! @brief const iterator of sequence entries
    typedef std::vector<SeqEntry>::const_iterator const_iterator;

    //! @brief Construct empty
    MultipleAlignment();
    
    /**
     * @brief Construct from file
     *
     * @param file name of input file
     * @param format file format (CLUSTAL or FASTA) 
     * @throw failure on read problems
     * @see MultipleAlignment(std::istream &in)
    */
    MultipleAlignment(const std::string &file, FormatType::type format=FormatType::CLUSTAL);

    /**
     * @brief Construct from stream
     *
     * @param in input stream with alignment in clustalW-like format
     * @param format file format (CLUSTAL or FASTA) 
     * @throw failure on read errors
    */
    MultipleAlignment(std::istream &in, FormatType::type format=FormatType::CLUSTAL);

    /**
     * @brief Construct as degenerate alignment of one sequence
     * @param name name of sequence
     * @param sequence sequence strings
     */
    MultipleAlignment(const std::string &name,
		      const std::string &sequence);
    
    /**
     * @brief Construct as pairwise alignment from names and alignment strings
     * @param nameA name of sequence A
     * @param nameB name of sequence B
     * @param alistringA alignment strings of sequence A
     * @param alistringB alignment strings of sequence B
     * 
     * @note handling of gap-symbols: use same gap symbols as in given alistrings
     */
    MultipleAlignment(const std::string &nameA,
		      const std::string &nameB,
		      const std::string &alistringA,
		      const std::string &alistringB);
    
    /**
     * @brief Construct from Alignment object
     * @param alignment object of type Alignment
     * @param only_local if true, construct only local alignment
     *
     * Automatically computes a consensus anchor string if anchors are
     * available. Consensus anchors containing duplicate names are cleared.
     * Does not compute some kind of consensus structure,
     * even if structure annotation of sequences A and B in Alignment
     * is available.
     */
    MultipleAlignment(const Alignment &alignment, bool only_local=false);
    
    /**
     * @brief Construct from alignment edges and sequences
     * @param edges alignment edges
     * @param seqA sequence A
     * @param seqB sequence B
     *
     * Automatically computes a consensus anchor string if anchors are
     * available. Consensus anchors containing duplicate names are
     * cleared. Does not compute some kind of consensus structure,
     * even if structure annotation of sequences A and B is available.
     */
    MultipleAlignment(const AlignmentEdges &edges,
		      const Sequence &seqA,
		      const Sequence &seqB);

protected:
    /**
     * @brief Initialize from alignment edges and sequences
     * @param edges alignment edges
     * @param seqA sequence A
     * @param seqB sequence B
     */
    void
    init(const AlignmentEdges &edges,
	 const Sequence &seqA,
	 const Sequence &seqB);
public:

    /**
     * @brief virtual destructor
     */
    virtual
    ~MultipleAlignment();
    
    /**
     * @brief "cast" multiple alignment to sequence
     *
     * @note this works like an upcast; this is ok, as long as
     * sequence does not specify attributes
     */
    const Sequence & as_sequence() const;  

    /**
     * @brief normalize rna symbols
     * @see normalize_rna_sequence()
     *
     * Normalize the symbols in all aligned sequences assuming that
     * they code for RNA
     */
    void
    normalize_rna_symbols();
    
    /**
     * @brief Number of rows of multiple aligment
     * @return number of rows
    */
    size_type
    num_of_rows() const { 
	return alig_.size();
    }
    
    /** 
     * @brief Emptiness check
     * 
     * @return whether the object contains no sequences 
     *
     * @note an alignment containing one or more empty sequences is
     * not empty in this sense.
     */
    bool
    empty() const {
	return alig_.empty();
    }

    /**
     * @brief Read access of annotation by prefix
     * @param type of annotation
     * @return sequence annotation
     * @note returns ref to empty annotation if annotation is not available
     */
    const SequenceAnnotation &
    annotation(const AnnoType::type &annotype) const;

    /**
     * @brief Write access to annotation
     * @param prefix annotation prefix
     * @param annotation sequence annotation
     * @todo check that annotation is valid for multiple alignment;
     * throw failure if annotation is not valid
     */
    void
    set_annotation(const AnnoType::type &annotype,
		   const SequenceAnnotation &annotation) {
	assert(0<=annotype && annotype<num_of_annotypes());
	annotations_[(size_t)annotype] = annotation;
    }
    
    /**
     * Annotation availability
     * @param prefix annotation prefix
     * @return wheter annotions with prefix are available
     */
    bool
    has_annotation(const AnnoType::type &annotype) const {
	assert(0<=annotype && annotype<num_of_annotypes());
	return annotations_.find(annotype)!=annotations_.end();
    }

    /**
     * @brief Test whether alignment is proper
     * @return whether all sequences have the same length
    */
    bool
    is_proper() const;
    
    /**
     * @brief Length of multiple aligment
     *
     * @note Assumes proper alignment. Does not check, whether all
     * sequences have the same length!
     * @return length of first sequence in alignment
    */
    pos_type 
    length() const { return alig_.empty() ? 0 : alig_[0].seq().length(); }
    
    /**
     * @brief Begin for read-only traversal of name/sequence pairs
     * @return begin iterator
    */
    const_iterator
    begin() const {
	return alig_.begin();
    }
    
    /**
     * @brief End for read-only traversal of name/sequence pairs
     * @return end iterator
    */
    const_iterator
    end() const {
	return alig_.end();
    }
    
    /**
     * @brief Test whether name exists
     * @param name name of a sequence
     * @return whether sequence with given name exists in multiple alignment
    */
    bool
    contains(std::string name) const;
      
    /* index access saves time over access by sequence name */
    
    /**
     * @brief Access index by name
     *
     * @pre name exists
     * @param name name of a sequence
     * @return index of name/sequence pair with given name
    */
    size_type
    index(const std::string &name) const {
	str2idx_map_t::const_iterator it = name2idx_.find(name);
	assert(it!=name2idx_.end());
	return it->second;
    }
    
    /**
     * @brief Access name/sequence pair by index
     *
     * @pre index in range 0..size()-1
     * @param index index of name/sequence pair (0-based)
     * @return sequence (including gaps) with given index
    */
    const SeqEntry &
    seqentry(size_type index) const {
	return alig_[index];
    }
    
    /**
     * @brief Access name/sequence pair by name
     *
     * @param name name of name/sequence pair
     * @return sequence (including gaps) with given name
    */
    const SeqEntry &
    seqentry(const std::string &name) const {
	return alig_[index(name)];
    }
    

    /**
     * @brief Deviation of a multiple alignment from a reference alignment
     *
     * @param ma multiple alignment
     * @return deviation of ma from reference alignment *this
     * deviation is defined for realignment in limited deviation from a
     * reference alignment as preformed when --max-diff-aln is given with
     * --max-diff to locarna.
     * @pre the sequences of ma have to occur in the alignment *this 
    */
    size_type
    deviation(const MultipleAlignment &ma) const; 
    
    /**
     * @brief Sum-of-pairs score between a multiple alignment and a reference alignment
     *
     * @param ma multiple alignment
     * @param compalign whether to compute score like compalign
     *
     * @return sum-of-pairs score of ma from reference alignment *this
     *
     * @note Whereas the sps score for compalign==FALSE
     * counts common matches only, the compalign score additionally
     * counts common indels.
     *
     * @pre the sequences of ma have to occur in the alignment *this 
    */
    double
    sps(const MultipleAlignment &ma, bool compalign=true) const; 
    
    /**
     * @brief Cmfinder realignment score of a multiple alignment to a reference alignment
     *
     * @param ma multiple alignment
     *
     * @return cmfinder realignment score of ma to reference alignment *this
     *
     * @note this score was defined in Elfar Torarinsson, Zizhen Yao,
     * Eric D. Wiklund, et al. Comparative genomics beyond
     * sequence-based alignments: RNA structures in the ENCODE
     * regions. Genome Res. 2008 (Section Realignment calculation)
     *
     * @pre the sequences of ma have to occur in the alignment *this 
    */
    double
    cmfinder_realignment_score(const MultipleAlignment &ma) const; 

    /** 
     * @brief Average deviation score
     * 
     * @param ma multiple alignment
     * 
     * @return average deviation fo alignment ma to reference alignment *this
     *
     * @pre the sequences of ma have to occur in the alignment *this 

     * @note this is not the same as deviation (and may be even
     * not very similar)!
     */
    double
    avg_deviation_score(const MultipleAlignment &ma) const;

    
    /** 
     * @brief Consensus sequence of multiple alignment
     * 
     * Consensus sequence by simple majority in each column. Assume
     * that only ascii < 127 characters occur
     *
     * @return consensus sequence as string
     */
    std::string
    consensus_sequence() const;
    
    /** 
     * @brief Access alignment column
     * 
     * @param col_index column index 
     * 
     * @return reference to alignment column with index i (1-based)
     */
    AliColumn
    column(size_type col_index) const {
	return AliColumn(*this,col_index);
    }

    /** 
     * @brief Append sequence entry
     * 
     * @param seqentry new sequence entry
     *
     * @pre *this is empty or entry must have same size as *this
     */
    void
    append(const SeqEntry &seqentry);

    /** 
     * @brief Prepend sequence entry
     * 
     * @param seqentry new sequence entry
     *
     * @pre *this is empty or entry must have same size as *this
     *
     * @note prepend is a lot more costly then append; it has cost
     * linearly in the number of rows
     */
    void
    prepend(const SeqEntry &seqentry);
    
    /**
     * @brief Append a column
     *
     * @param c column that is appended
     */
    void
    operator += (const AliColumn &c);
    
    /**
     * @brief Append the same character to each row
     *
     * @param c character that is appended
     */
    void
    operator += (char c);
    
    /**
     * @brief reverse the multiple alignment
     */
    void
    reverse();


    // ------------------------------------------------------------
    // output
    
    /**
     * @brief Write alignment to stream
     *
     * @param out output stream
     * @return output stream
     *
     * Writes one line "<name> <seq>" for each sequence.
     */
    std::ostream &
    write(std::ostream &out) const;

    /**
     * @brief Write alignment to stream
     *
     * @param out output stream
     * @param width output stream
     * @return output stream
     *
     * Writes lines "<name> <seq>" per sequence, wraps lines at width
     */
    std::ostream &
    write(std::ostream &out, size_t width) const;
    
    /**
     * @brief Write formatted line of name and sequence
     *
     * The line is formatted such that it fits the output of the write
     * methods.
     *
     * @param out output stream
     * @param name name string
     * @param sequence sequence string
     * @return output stream
     */
    std::ostream &
    write_name_sequence_line(std::ostream &out,
			     const std::string &name,
			     const std::string &sequence) const;
    
    /**
     * @brief Write sub-alignment to stream 
     *
     * Write from position start to position end to output stream
     * out; write lines "<name> <seq>"
     *
     * @param out output stream
     * @param start start column (1-based)
     * @param end end column (1-based)
     * @return output stream
     */
    std::ostream &
    write(std::ostream &out, size_type start, size_type end) const;
    
    /**
     * @brief check character constraints
     *
     * Check whether the alignment contains characters from the given
     * alphabet only and, if warn, print warnings otherwise.
     *
     * @param alphabet alphabet of admissible characters
     *
     * @return whether all characters are in the alphabet
     */
    bool 
    checkAlphabet(const Alphabet<char> &alphabet) const;
    
private:
        
    /**
     * @brief Deviation of a pairwise alignment from a pairwise reference alignment
     * @param a1 first alignment string of alignment a
     * @param a2 second alignment string of alignment a
     * @param ref1 first alignment string of reference alignment ref
     * @param ref2 second alignment string of reference alignment ref
     *
     * @return deviation of alignment a from reference alignment ref
     */
    static
    size_type
    deviation2(const string1 &a1,
	       const string1 &a2,
	       const string1 &ref1,
	       const string1 &ref2
	       );

    
    /** 
    * @brief Pairwise match score for calculation of match_sps
    * 
    * @param a1 row 1 of test alignment
    * @param a2 row 2 of test alignment
    * @param ref1 row 1 of reference alignment
    * @param ref2 row 2 of reference alignment
    * @param score_common_gaps whehter to score common gaps
    * 
    * @return alignment comparison match score for pairwise alignments (a1,a2) and (ref1,ref2)
    *
    * @see sps()
    */
    static
    double
    pairwise_match_score(const SeqEntry &a1,
			 const SeqEntry &a2,
			 const SeqEntry &ref1,
			 const SeqEntry &ref2,
			 bool score_common_gaps
			 );
    
    /** 
     * @brief Determine matching positions for each string position
     * 
     * @param s string 1
     * @param t string 2
     * 
     * @return vector v of length length(s+1), such that for each
     * position i in s (1<=i<=|s|), v[i] is the matching position in t
     * or -1 if there is no match.
     */
    static
    std::vector<int>
    match_vector(const string1 &s,
		 const string1 &t);
    
    /** 
     * @brief Determine matching positions for each string position
     * 
     * @param s string 1
     * @param t string 2
     * 
     * @return vector v of length length(s+1), such that for each
     * position i in s (1<=i<=|s|), v[i] is the matching position in t
     * or the position after that i is deleted.
     */
    static
    std::vector<int>
    match_vector2(const string1 &s,
		  const string1 &t);

    /** 
     * @brief Count matches in pairwise alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string 2
     * 
     * @return number of matches
     */
    static
    size_t
    count_matches(const SeqEntry &a1,
		  const SeqEntry &a2);
    
    /** 
     * @brief Count matches in pairwise alignment that do not occur in
     * a second alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string1 
     * @param ref1 reference alignment string 1
     * @param ref2 reference alignment string 1
     * 
     * @return number of matches exclusively in alignment a (and not
     * in reference)
     */
    static
    size_t
    count_exclusive_matches(const SeqEntry &a1,
			    const SeqEntry &a2,
			    const SeqEntry &ref1,
			    const SeqEntry &ref2
			    );

    /** 
     * @brief Average deviation score for pairwise alignment
     * 
     * @param a1 alignment string 1
     * @param a2 alignment string1 
     * @param ref1 reference alignment string 1
     * @param ref2 reference alignment string 1
     * 
     * @return avg deviation score for alignment (a1,a2) from
     * reference alignment (ref1,ref2)
     *
     * @note This score averages over the differences of positions j_a
     * and j_ref for all positions i and computes a sum-of-pairs
     * score. In case i is matched to a gap between j^left and
     * j^right, we define j as the average value (for j in {j_a,
     * j_ref}).
     */
    static
    double
    pairwise_deviation_score(const SeqEntry &a1,
			     const SeqEntry &a2,
			     const SeqEntry &ref1,
			     const SeqEntry &ref2
			     );

public:

    /** 
     * @brief Print contents of object to stream
     * @param out output stream
     */
    void
    write_debug(std::ostream &out=std::cout) const;
};
    
    /**
     * @brief Write multiple alignment to stream
     * @param out output stream
     * @param ma multiple alignment
     * @return output stream
     */
    std::ostream &
    operator << (std::ostream &out, const MultipleAlignment &ma);

} // end namespace

#endif // LOCARNA_MULTIPLE_ALIGNMENT_HH
