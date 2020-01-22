
from distutils.spawn import find_executable
import subprocess
import statistics
import random
import gzip
import uuid
import sys
import re
import os

"""

Run doctests:

python3 -m doctest gplib.py


"""


################################################################################

def graphprot_predictions_get_median(predictions_file):
    """
    Given a GraphProt .predictions file, read in site scores and return 
    the median value.

    >>> test_file = "test-data/test.predictions"
    >>> graphprot_predictions_get_median(test_file)
    0.571673

    """
    # Site scores list.
    sc_list = []
    with open(predictions_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            score = float(cols[2])
            sc_list.append(score)
    f.close()
    # Return the median.
    return statistics.median(sc_list)


################################################################################

def graphprot_profile_get_top_scores_median(profile_file,
                                            profile_type="profile",
                                            avg_profile_extlr=5):

    """
    Given a GraphProt .profile file, extract for each site (identified by 
    column 1 ID) the top (= highest) score. Then return the median of these 
    top scores.
    
    profile_type can be either "profile" or "avg_profile".
    "avg_profile means that the position-wise scores will first get smoothed 
    out by calculating for each position a new score through taking a 
    sequence window -avg_profile_extlr to +avg_profile_extlr of the position 
    and calculate the mean score over this window and assign it to the position.
    After that, the maximum score of each site is chosen, and the median over 
    all maximum scores is returned.
    "profile" leaves the position-wise scores as they are, directly extracting 
    the maximum for each site and then reporting the median.
    
    >>> test_file = "test-data/test.profile"
    >>> graphprot_profile_get_top_scores_median(test_file)
    3.2

    """
    # Dictionary of lists, with list of scores (value) for each site (key).
    lists_dic = {}
    with open(profile_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            score = float(cols[2])
            if seq_id in lists_dic:
                lists_dic[seq_id].append(score)
            else:
                lists_dic[seq_id] = []
                lists_dic[seq_id].append(score)
    f.close()
    # For each site, extract maximum and store in new list.
    max_list = []
    for seq_id in lists_dic:
        if profile_type == "profile":
            max_sc = max(lists_dic[seq_id])
            max_list.append(max_sc)
        elif profile_type == "avg_profile":
            # Convert profile score list to average profile scores list.
            aps_list = list_moving_window_average_values(lists_dic[seq_id],
                                                         win_extlr=avg_profile_extlr)
            max_sc = max(aps_list)
            max_list.append(max_sc)
        else:
            assert 0, "invalid profile_type argument given: \"%s\"" %(profile_type)
    # Return the median.
    return statistics.median(max_list)


################################################################################
#!/usr/bin/env python3
def list_moving_window_average_values(in_list, 
                                      win_extlr=5,
                                      method=1):
    """
    Take a list of numeric values, and calculate for each position a new value, 
    by taking the mean value of the window of positions -win_extlr and 
    +win_extlr. If full extension is not possible (at list ends), it just 
    takes what it gets.
    Two implementations of the task are given, chose by method=1 or method=2.

    >>> test_list = [2, 3, 5, 8, 4, 3, 7, 1]
    >>> list_moving_window_average_values(test_list, win_extlr=2, method=1)
    [3.3333333333333335, 4.5, 4.4, 4.6, 5.4, 4.6, 3.75, 3.6666666666666665]
    >>> list_moving_window_average_values(test_list, win_extlr=2, method=2)
    [3.3333333333333335, 4.5, 4.4, 4.6, 5.4, 4.6, 3.75, 3.6666666666666665]

    """
    l_list = len(in_list)
    assert l_list, "Given list is empty"
    new_list = [0] * l_list
    if win_extlr == 0:
        return l_list
    if method == 1:
        for i in range(l_list):
            s = i - win_extlr
            e = i + win_extlr + 1
            if s < 0:
                s = 0
            if e > l_list:
                e = l_list
            # Extract portion and assign value to new list.
            new_list[i] = statistics.mean(in_list[s:e])
    elif method == 2:
        for i in range(l_list):
            s = i - win_extlr
            e = i + win_extlr + 1
            if s < 0:
                s = 0
            if e > l_list:
                e = l_list
            l = e-s
            sc_sum = 0
            for j in range(l):
                sc_sum += in_list[s+j]
            new_list[i] = sc_sum / l
    else:
        assert 0, "invalid method ID given (%i)" %(method)
    return new_list


################################################################################

def echo_add_to_file(echo_string, out_file):
    """
    Add a string to file, using echo command.

    """
    check_cmd = 'echo "%s" >> %s' % (echo_string, out_file)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "echo is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def count_fasta_headers(fasta_file):
    """
    Count number of FASTA headers in fasta_file using grep.

    >>> test_file = "test-data/test.fa"
    >>> count_fasta_headers(test_file)
    2
    >>> test_file = "test-data/empty_file"
    >>> count_fasta_headers(test_file)
    0

    """
    check_cmd = 'grep -c ">" ' + fasta_file
    output = subprocess.getoutput(check_cmd)
    row_count = int(output.strip())
    return row_count


################################################################################

def make_file_copy(in_file, out_file):
    """
    Make a file copy by copying in_file to out_file.

    """
    check_cmd = "cat " + in_file + " > " + out_file
    assert in_file != out_file, "cat does not like to cat file into same file (%s)" %(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "cat did not like your input (in_file: %s, out_file: %s):\n%s" %(in_file, out_file, output)


################################################################################

def split_fasta_into_test_train_files(in_fasta, test_out_fa, train_out_fa, 
                                      test_size=500):
    """
    Split in_fasta .fa file into two files (e.g. test, train).

    """
    # Read in in_fasta.
    seqs_dic = read_fasta_into_dic(in_fasta)
    # Shuffle IDs.
    rand_ids_list = random_order_dic_keys_into_list(seqs_dic)
    c_out = 0
    TESTOUT = open(test_out_fa, "w")
    TRAINOUT = open(train_out_fa, "w")
    for seq_id in rand_ids_list:
        seq = seqs_dic[seq_id]
        if (c_out >= test_size):
            TRAINOUT.write(">%s\n%s\n" % (seq_id, seq))
        else:
            TESTOUT.write(">%s\n%s\n" % (seq_id, seq))
        c_out += 1
    TESTOUT.close()
    TRAINOUT.close()


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        read_dna=False,
                        reject_lc=False,
                        convert_to_uc=False,
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, convert to RNA, store in dictionary 
    and return dictionary.
    
    >>> test_fasta = "test-data/test.fa"
    >>> read_fasta_into_dic(test_fasta)
    {'seq1': 'acguACGUacgu', 'seq2': 'ugcaUGCAugcaACGUacgu'}
    >>> test_fasta = "test-data/test2.fa"
    >>> read_fasta_into_dic(test_fasta)
    {}
    >>> test_fasta = "test-data/test.ensembl.fa"
    >>> read_fasta_into_dic(test_fasta, read_dna=True)
    {'ENST00000415118': 'GAAATAGT', 'ENST00000448914': 'ACTGGGGGATACGAAAA'}

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""
    seq = ""

    # Go through FASTA file, extract sequences.
    if re.search(".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else: 
        f = open(fasta_file, "r")
    for line in f:
        if re.search(">.+", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            # If there is a ".", take only first part of header.
            # This assumes ENSEMBL header format ">ENST00000631435.1 cdna ..."
            if re.search(".+\..+", seq_id):
                m = re.search("(.+?)\..+", seq_id)
                seq_id = m.group(1)
            assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)
            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
        elif re.search("[ACGTUN]+", line, re.I):
            if seq_id in seqs_dic:
                m = re.search("([ACGTUN]+)", line, re.I)
                seq = m.group(1)
                if reject_lc:
                    assert not re.search("[a-z]", seq), "lowercase characters detected in sequence \"%i\" (reject_lc=True)" %(seq_id)
                if convert_to_uc:
                    seq = seq.upper()
                # If sequences with N nucleotides should be skipped.
                if skip_n_seqs:
                    if "n" in m.group(1) or "N" in m.group(1):
                        print ("WARNING: \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id))
                        del seqs_dic[seq_id]
                        continue
                # Convert to RNA, concatenate sequence.
                if read_dna:
                    seqs_dic[seq_id] += m.group(1).replace("U","T").replace("u","t")
                else:
                    seqs_dic[seq_id] += m.group(1).replace("T","U").replace("t","u")
    f.close()
    return seqs_dic


################################################################################

def random_order_dic_keys_into_list(in_dic):
    """
    Read in dictionary keys, and return random order list of IDs.

    """
    id_list = []
    for key in in_dic:
        id_list.append(key)
    random.shuffle(id_list)
    return id_list


################################################################################

def graphprot_get_param_string(params_file):
    """
    Get parameter string from GraphProt .params file.

    >>> test_params = "test-data/test.params"
    >>> graphprot_get_param_string(test_params)
    '-epochs 20 -lambda 0.01 -R 1 -D 3 -bitsize 14 -onlyseq '

    """
    param_string = ""
    with open(params_file) as f:
        for line in f:
            cols = line.strip().split(" ")
            param = cols[0]
            setting = cols[1]
            if re.search(".+:", param):
                m = re.search("(.+):", line)
                par = m.group(1)
                if re.search("pos_train.+", line):
                    continue
                if par == "model_type":
                    if setting == "sequence":
                        param_string += "-onlyseq "
                else:
                    param_string += "-%s %s " %(par, setting)
            else:
                assert 0, "pattern matching failed for string \"%s\"" %(param)
    return param_string


################################################################################

def seqs_dic_count_uc_nts(seqs_dic):
    """
    Count number of uppercase nucleotides in sequences stored in sequence 
    dictionary.
    
    >>> seqs_dic = {'seq1': "acgtACGTacgt", 'seq2': 'acgtACacgt'}
    >>> seqs_dic_count_uc_nts(seqs_dic)
    6
    >>> seqs_dic = {'seq1': "acgtacgt", 'seq2': 'acgtacgt'}
    >>> seqs_dic_count_uc_nts(seqs_dic)
    0

    """
    assert seqs_dic, "Given sequence dictionary empty"
    c_uc = 0
    for seq_id in seqs_dic:
        c_uc += len(re.findall(r'[A-Z]', seqs_dic[seq_id]))
    return c_uc


################################################################################

def seqs_dic_count_lc_nts(seqs_dic):
    """
    Count number of lowercase nucleotides in sequences stored in sequence 
    dictionary.
    
    >>> seqs_dic = {'seq1': "gtACGTac", 'seq2': 'cgtACacg'}
    >>> seqs_dic_count_lc_nts(seqs_dic)
    10
    >>> seqs_dic = {'seq1': "ACGT", 'seq2': 'ACGTAC'}
    >>> seqs_dic_count_lc_nts(seqs_dic)
    0

    """
    assert seqs_dic, "Given sequence dictionary empty"
    c_uc = 0
    for seq_id in seqs_dic:
        c_uc += len(re.findall(r'[a-z]', seqs_dic[seq_id]))
    return c_uc


################################################################################

def count_file_rows(in_file):
    """
    Count number of file rows for given input file.
    
    >>> test_file = "test-data/test1.bed"
    >>> count_file_rows(test_file)
    7
    >>> test_file = "test-data/empty_file"
    >>> count_file_rows(test_file)
    0

    """
    check_cmd = "cat " + in_file + " | wc -l"
    output = subprocess.getoutput(check_cmd)
    row_count = int(output.strip())
    return row_count


################################################################################

def bed_check_six_col_format(bed_file):
    """
    Check whether given .bed file has 6 columns.

    >>> test_bed = "test-data/test1.bed"
    >>> bed_check_six_col_format(test_bed)
    True
    >>> test_bed = "test-data/empty_file"
    >>> bed_check_six_col_format(test_bed)
    False

    """

    six_col_format = False
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) == 6:
                six_col_format = True
            break
    f.closed
    return six_col_format


################################################################################

def bed_check_unique_ids(bed_file):
    """
    Check whether .bed file (6 column format with IDs in column 4) 
    has unique column 4 IDs.
    
    >>> test_bed = "test-data/test1.bed"
    >>> bed_check_unique_ids(test_bed)
    True
    >>> test_bed = "test-data/test2.bed"
    >>> bed_check_unique_ids(test_bed)
    False

    """

    check_cmd = "cut -f 4 " + bed_file + " | sort | uniq -d"
    output = subprocess.getoutput(check_cmd)
    if output:
        return False
    else:
        return True


################################################################################

def get_seq_lengths_from_seqs_dic(seqs_dic):
    """
    Given a dictionary of sequences, return dictionary of sequence lengths.
    Mapping is sequence ID -> sequence length.
    """
    seq_len_dic = {}
    assert seqs_dic, "sequence dictionary seems to be empty"
    for seq_id in seqs_dic:
        seq_l = len(seqs_dic[seq_id])
        seq_len_dic[seq_id] = seq_l
    return seq_len_dic


################################################################################

def bed_get_region_lengths(bed_file):
    """
    Read in .bed file, store and return region lengths in dictionary.
    key   :  region ID (.bed col4)
    value :  region length (.bed col3-col2)

    >>> test_file = "test-data/test4.bed"
    >>> bed_get_region_lengths(test_file)
    {'CLIP1': 10, 'CLIP2': 10}

    """
    id2len_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_l = site_e - site_s
            assert site_id not in id2len_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(bed_file)
            id2len_dic[site_id] = site_l
    f.closed
    assert id2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_bed)
    return id2len_dic


################################################################################

def graphprot_get_param_dic(params_file):
    """
    Read in GraphProt .params file and store in dictionary.
    key = parameter
    value = parameter value

    >>> params_file = "test-data/test.params"
    >>> graphprot_get_param_dic(params_file)
    {'epochs': '20', 'lambda': '0.01', 'R': '1', 'D': '3', 'bitsize': '14', 'model_type': 'sequence', 'pos_train_ws_pred_median': '0.760321', 'pos_train_profile_median': '5.039610', 'pos_train_avg_profile_median_1': '4.236340', 'pos_train_avg_profile_median_2': '3.868431', 'pos_train_avg_profile_median_3': '3.331277', 'pos_train_avg_profile_median_4': '2.998667', 'pos_train_avg_profile_median_5': '2.829782', 'pos_train_avg_profile_median_6': '2.626623', 'pos_train_avg_profile_median_7': '2.447083', 'pos_train_avg_profile_median_8': '2.349919', 'pos_train_avg_profile_median_9': '2.239829', 'pos_train_avg_profile_median_10': '2.161676'}

    """
    param_dic = {}
    with open(params_file) as f:
        for line in f:
            cols = line.strip().split(" ")
            param = cols[0]
            setting = cols[1]
            if re.search(".+:", param):
                m = re.search("(.+):", line)
                par = m.group(1)
                param_dic[par] = setting
    f.close()
    return param_dic


################################################################################

def graphprot_filter_predictions_file(in_file, out_file,
                                      sc_thr=0):
    """
    Filter GraphProt .predictions file by given score thr_sc.
    """
    OUTPRED = open(out_file, "w")
    with open(in_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            score = float(cols[2])
            if score < sc_thr:
                continue
            OUTPRED.write("%s\n" %(row))
    f.close()
    OUTPRED.close()


################################################################################

def fasta_read_in_ids(fasta_file):
    """
    Given a .fa file, read in header IDs in order appearing in file, 
    and store in list.

    >>> test_file = "test-data/test3.fa"
    >>> fasta_read_in_ids(test_file)
    ['SERBP1_K562_rep01_544', 'SERBP1_K562_rep02_709', 'SERBP1_K562_rep01_316']

    """
    ids_list = []
    with open(fasta_file) as f:
        for line in f:
            if re.search(">.+", line):
                m = re.search(">(.+)", line)
                seq_id = m.group(1)
                ids_list.append(seq_id)
    f.close()
    return ids_list


################################################################################

def graphprot_profile_calculate_avg_profile(in_file, out_file,
                                            ap_extlr=5,
                                            seq_ids_list=False,
                                            method=1):
    """
    Given a GraphProt .profile file, calculate average profiles and output 
    average profile file.
    Average profile means that the position-wise scores will get smoothed 
    out by calculating for each position a new score, taking a sequence 
    window -ap_extlr to +ap_extlr relative to the position 
    and calculate the mean score over this window. The mean score then 
    becomes the new average profile score at this position.
    Two different implementations of the task are given:
    method=1 (new python implementation, slower + more memory but easy to read)
    method=2 (old perl implementation, faster and less memory but more code)

    >>> in_file = "test-data/test2.profile"
    >>> out_file1 = "test-data/test2_1.avg_profile"
    >>> out_file2 = "test-data/test2_2.avg_profile"
    >>> out_file4 = "test-data/test2_3.avg_profile"
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file1, ap_extlr=2, method=1)
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file2, ap_extlr=2, method=2)
    >>> diff_two_files_identical(out_file1, out_file2)
    True
    >>> test_list = ["s1", "s2", "s3", "s4"]
    >>> out_file3_exp = "test-data/test3_added_ids_exp.avg_profile"
    >>> out_file3 = "test-data/test3_added_ids_out.avg_profile"
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file3, ap_extlr=2, method=1, seq_ids_list=test_list)
    >>> diff_two_files_identical(out_file3_exp, out_file3)
    True

    """
    if method == 1:
        # Dictionary of lists, with list of scores (value) for each site (key).
        lists_dic = {}
        site_starts_dic = {}
        with open(in_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                site_id = int(cols[0])
                pos = int(cols[1]) # 0-based.
                score = float(cols[2])
                # Store first position of site.
                if site_id not in site_starts_dic:
                    site_starts_dic[site_id] = pos
                if site_id in lists_dic:
                    lists_dic[site_id].append(score)
                else:
                    lists_dic[site_id] = []
                    lists_dic[site_id].append(score)
        f.close()
        # Check number of IDs (# FASTA sequence IDs has to be same as # site IDs).
        if seq_ids_list:
            c_seq_ids = len(seq_ids_list)
            c_site_ids = len(site_starts_dic)
            assert c_seq_ids == c_site_ids, "# sequence IDs != # site IDs (%i != %i)" %(c_seq_ids, c_site_ids)
        OUTPROF = open(out_file, "w")
        # For each site, calculate average profile scores list.
        max_list = []
        for site_id in lists_dic:
            # Convert profile score list to average profile scores list.
            aps_list = list_moving_window_average_values(lists_dic[site_id],
                                                         win_extlr=ap_extlr)
            start_pos = site_starts_dic[site_id]
            # Get original FASTA sequence ID.
            if seq_ids_list:
                site_id = seq_ids_list[site_id]
            for i, sc in enumerate(aps_list):
                pos = i + start_pos + 1 # make 1-based.
                OUTPROF.write("%s\t%i\t%f\n" %(site_id, pos, sc))
        OUTPROF.close()
    elif method == 2:
        OUTPROF = open(out_file, "w")
        # Old site ID.
        old_id = ""
        # Current site ID.
        cur_id = ""
        # Scores list.
        scores_list = []
        site_starts_dic = {}
        with open(in_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                cur_id = int(cols[0])
                pos = int(cols[1]) # 0-based.
                score = float(cols[2])
                # Store first position of site.
                if cur_id not in site_starts_dic:
                    site_starts_dic[cur_id] = pos
                # Case: new site (new column 1 ID).
                if cur_id != old_id:
                    # Process old id scores.
                    if scores_list:
                        aps_list = list_moving_window_average_values(scores_list,
                                                                     win_extlr=ap_extlr)
                        start_pos = site_starts_dic[old_id]
                        seq_id = old_id
                        # Get original FASTA sequence ID.
                        if seq_ids_list:
                            seq_id = seq_ids_list[old_id]
                        for i, sc in enumerate(aps_list):
                            pos = i + start_pos + 1 # make 1-based.
                            OUTPROF.write("%s\t%i\t%f\n" %(seq_id, pos, sc))
                        # Reset list.
                        scores_list = []
                    old_id = cur_id
                    scores_list.append(score)
                else:
                    # Add to scores_list.
                    scores_list.append(score)
        f.close()
        # Process last block.
        if scores_list:
            aps_list = list_moving_window_average_values(scores_list,
                                                         win_extlr=ap_extlr)
            start_pos = site_starts_dic[old_id]
            seq_id = old_id
            # Get original FASTA sequence ID.
            if seq_ids_list:
                seq_id = seq_ids_list[old_id]
            for i, sc in enumerate(aps_list):
                pos = i + start_pos + 1 # make 1-based.
                OUTPROF.write("%s\t%i\t%f\n" %(seq_id, pos, sc))
        OUTPROF.close()


################################################################################

def graphprot_profile_extract_peak_regions(in_file, out_file,
                                           max_merge_dist=0,
                                           sc_thr=0):
    """
    Extract peak regions from GraphProt .profile file.
    Store the peak regions (defined as regions with scores >= sc_thr) 
    as to out_file in 6-column .bed.

    TODO:
    Add option for genomic coordinates input (+ - polarity support).
    Output genomic regions instead of sequence regions.

    >>> in_file = "test-data/test4.avg_profile"
    >>> out_file = "test-data/test4_out.peaks.bed"
    >>> exp_file = "test-data/test4_out_exp.peaks.bed"
    >>> exp2_file = "test-data/test4_out_exp2.peaks.bed"
    >>> empty_file = "test-data/empty_file"
    >>> graphprot_profile_extract_peak_regions(in_file, out_file)
    >>> diff_two_files_identical(out_file, exp_file)
    True
    >>> graphprot_profile_extract_peak_regions(in_file, out_file, sc_thr=10)
    >>> diff_two_files_identical(out_file, empty_file)
    True
    >>> graphprot_profile_extract_peak_regions(in_file, out_file, max_merge_dist=2)
    >>> diff_two_files_identical(out_file, exp2_file)
    True

    """

    OUTPEAKS = open(out_file, "w")
    # Old site ID.
    old_id = ""
    # Current site ID.
    cur_id = ""
    # Scores list.
    scores_list = []
    site_starts_dic = {}
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            cur_id = cols[0]
            pos = int(cols[1]) # 0-based.
            score = float(cols[2])
            # Store first position of site.
            if cur_id not in site_starts_dic:
                # If first position != zero, we assume positions are 1-based.
                if pos != 0:
                    # Make index 0-based.
                    site_starts_dic[cur_id] = pos - 1
                else:
                    site_starts_dic[cur_id] = pos
            # Case: new site (new column 1 ID).
            if cur_id != old_id:
                # Process old id scores.
                if scores_list:
                    # Extract peaks from region.
                    peak_list = list_extract_peaks(scores_list, 
                                                   max_merge_dist=max_merge_dist,
                                                   coords="bed",
                                                   sc_thr=sc_thr)
                    start_pos = site_starts_dic[old_id]
                    # Print out peaks in .bed format.
                    for l in peak_list:
                        peak_s = start_pos + l[0]
                        peak_e = start_pos + l[1]
                        site_id = "%s,%i" %(old_id, l[2])
                        OUTPEAKS.write("%s\t%i\t%i\t%s\t%f\t+\n" %(old_id, peak_s, peak_e, site_id, l[3]))
                    # Reset list.
                    scores_list = []
                old_id = cur_id
                scores_list.append(score)
            else:
                # Add to scores_list.
                scores_list.append(score)
    f.close()
    # Process last block.
    if scores_list:
        # Extract peaks from region.
        peak_list = list_extract_peaks(scores_list, 
                                       max_merge_dist=max_merge_dist,
                                       coords="bed",
                                       sc_thr=sc_thr)
        start_pos = site_starts_dic[old_id]
        # Print out peaks in .bed format.
        for l in peak_list:
            peak_s = start_pos + l[0]
            peak_e = start_pos + l[1]
            site_id = "%s,%i" %(old_id, l[2]) # best score also 1-based.
            OUTPEAKS.write("%s\t%i\t%i\t%s\t%f\t+\n" %(old_id, peak_s, peak_e, site_id, l[3]))
    OUTPEAKS.close()


################################################################################

def list_extract_peaks(in_list,
                       max_merge_dist=0,
                       coords="list",
                       sc_thr=0):
    """
    Extract peak regions from list. 
    Peak region is defined as region >= score threshold.
    
    coords=bed  :  peak start 0-based, peak end 1-based.
    coords=list :  peak start 0-based, peak end 0-based.
    
    >>> test_list = [-1, 0, 2, 4.5, 1, -1, 5, 6.5]
    >>> list_extract_peaks(test_list)
    [[1, 4, 3, 4.5], [6, 7, 7, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=2)
    [[2, 3, 3, 4.5], [6, 7, 7, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=2, coords="bed")
    [[2, 4, 4, 4.5], [6, 8, 8, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=10)
    []
    >>> test_list = [2, -1, 3, -1, 4, -1, -1, 6, 9]
    >>> list_extract_peaks(test_list, max_merge_dist=2)
    [[0, 4, 4, 4], [7, 8, 8, 9]]
    >>> list_extract_peaks(test_list, max_merge_dist=3)
    [[0, 8, 8, 9]]

    """
    # Check.
    assert len(in_list), "Given list is empty"
    # Peak regions list.
    peak_list = []
    # Help me.
    inside = False
    pr_s = 0
    pr_e = 0
    pr_top_pos = 0
    pr_top_sc = -100000
    for i, sc in enumerate(in_list):
        # Part of peak region?
        if sc >= sc_thr:
            # At peak start.
            if not inside:
                pr_s = i
                pr_e = i
                inside = True
            else:
                # Inside peak region.
                pr_e = i
            # Store top position.
            if sc > pr_top_sc:
                pr_top_sc = sc
                pr_top_pos = i
        else:
            # Before was peak region?
            if inside:
                # Store peak region.
                #peak_infos = "%i,%i,%i,%f" %(pr_s, pr_e, pr_top_pos, pr_top_sc)
                peak_infos = [pr_s, pr_e, pr_top_pos, pr_top_sc]
                peak_list.append(peak_infos)
                inside = False
                pr_top_pos = 0
                pr_top_sc = -100000
    # If peak at the end, also report.
    if inside:
        # Store peak region.
        peak_infos = [pr_s, pr_e, pr_top_pos, pr_top_sc]
        peak_list.append(peak_infos)
    # Merge peaks.
    if max_merge_dist and len(peak_list) > 1:
        iterate = True
        while iterate:
            merged_peak_list = []
            added_peaks_dic = {}
            peaks_merged = False
            for i, l in enumerate(peak_list):
                if i in added_peaks_dic:
                    continue
                j = i + 1
                # Last element.
                if j == len(peak_list):
                    if i not in added_peaks_dic:
                        merged_peak_list.append(peak_list[i])
                    break
                # Compare two elements.
                new_peak = []
                if (peak_list[j][0] - peak_list[i][1]) <= max_merge_dist:
                    peaks_merged = True
                    new_top_pos = peak_list[i][2]
                    new_top_sc = peak_list[i][3]
                    if peak_list[i][3] < peak_list[j][3]:
                        new_top_pos = peak_list[j][2]
                        new_top_sc = peak_list[j][3]
                    new_peak = [peak_list[i][0], peak_list[j][1], new_top_pos, new_top_sc]
                # If two peaks were merged.
                if new_peak:
                    merged_peak_list.append(new_peak)
                    added_peaks_dic[i] = 1
                    added_peaks_dic[j] = 1
                else:
                    merged_peak_list.append(peak_list[i])
                    added_peaks_dic[i] = 1
            if not peaks_merged:
                iterate = False
            peak_list = merged_peak_list
            peaks_merged = False
    # If peak coordinates should be in .bed format, make peak ends 1-based.
    if coords == "bed":
        for i in range(len(peak_list)):
            peak_list[i][1] += 1
            peak_list[i][2] += 1 # 1-base best score position too.
    return peak_list


################################################################################

def bed_peaks_to_genomic_peaks(peak_file, genomic_peak_file, genomic_sites_bed, print_rows=False):
    """
    Given a .bed file of sequence peak regions (possible coordinates from 
    0 to length of s), convert peak coordinates to genomic coordinates.
    Do this by taking genomic regions of sequences as input.

    >>> test_in = "test-data/test.peaks.bed"
    >>> test_exp = "test-data/test_exp.peaks.bed"
    >>> test_out = "test-data/test_out.peaks.bed"
    >>> gen_in = "test-data/test.peaks_genomic.bed"
    >>> bed_peaks_to_genomic_peaks(test_in, test_out, gen_in)
    >>> diff_two_files_identical(test_out, test_exp)
    True

    """
    # Read in genomic region info.
    id2row_dic = {}

    with open(genomic_sites_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_id = cols[3]
            assert site_id not in id2row_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(args.genomic_sites_bed)
            id2row_dic[site_id] = row
    f.close()

    # Read in peaks file and convert coordinates.
    OUTPEAKS = open(genomic_peak_file, "w")
    with open(peak_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id2 = cols[3]
            site_sc = float(cols[4])
            assert re.search(".+,.+", site_id2), "regular expression failed for ID \"%s\"" %(site_id2)
            m = re.search(".+,(\d+)", site_id2)
            sc_pos = int(m.group(1)) # 1-based.
            assert site_id in id2row_dic, "site ID \"%s\" not found in genomic sites dictionary" %(site_id)
            row = id2row_dic[site_id]
            rowl = row.split("\t")
            gen_chr = rowl[0]
            gen_s = int(rowl[1])
            gen_e = int(rowl[2])
            gen_pol = rowl[5]
            new_s = site_s + gen_s
            new_e = site_e + gen_s
            new_sc_pos = sc_pos + gen_s
            if gen_pol == "-":
                new_s = gen_e - site_e
                new_e = gen_e - site_s
                new_sc_pos = gen_e - sc_pos + 1 # keep 1-based.
            new_row = "%s\t%i\t%i\t%s,%i\t%f\t%s" %(gen_chr, new_s, new_e, site_id, new_sc_pos, site_sc, gen_pol)
            OUTPEAKS.write("%s\n" %(new_row))
            if print_rows:
                print(new_row)
    OUTPEAKS.close()


################################################################################

def diff_two_files_identical(file1, file2):
    """
    Check whether two files are identical. Return true if diff reports no 
    differences.
    
    >>> file1 = "test-data/file1"
    >>> file2 = "test-data/file2"
    >>> diff_two_files_identical(file1, file2)
    True
    >>> file1 = "test-data/test1.bed"
    >>> diff_two_files_identical(file1, file2)
    False

    """
    same = True
    check_cmd = "diff " + file1 + " " + file2
    output = subprocess.getoutput(check_cmd)
    if output:
        same = False
    return same


################################################################################


