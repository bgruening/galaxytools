
from galaxy.tools.parameters import DataToolParameter

def get_matrix_header( input_dataset ):
    """
        Not used currently, because the reload of the ckeckboxes did not work.
    """
    input_handle = open( input_dataset.file_name )
    first_header = input_handle.readline()
    second_header = input_handle.readline()
    return [('%s::%s' % (cname2,cname1), str(int(col_num) + 1), False) for col_num, (cname2, cname1) in enumerate(zip(second_header.split()[1:],first_header.split()[1:])) ]


def validate_input( trans, error_map, param_values, page_param_map ):
    """
        Validates the user input, before execution.
    """
    factors = param_values['rep_factorName']
    factor_name_list = []
    factor_duplication = False
    level_duplication = False
    overlapping_selection = False


    for factor in factors:
        # factor names should be unique
        fn = factor['factorName']
        if fn in factor_name_list:
            factor_duplication = True
            break
        factor_name_list.append( fn )

        level_name_list = list()
        factor_index_list = list()

        for level in factor['rep_factorLevel']:
            # level names under one factor should be unique
            fl = level['factorLevel']
            if fl in level_name_list:
                level_duplication = True
            level_name_list.append( fl )

            fi = level['factorIndex']
            if fi:
                # the checkboxes should not have an overlap
                for check in fi:
                    if check in factor_index_list:
                        overlapping_selection = True
                    factor_index_list.append( check )

        if level_duplication:
            error_map['rep_factorName'] = [ dict() for t in factors ]
            for i in range( len( factors ) ):
                error_map['rep_factorName'][i]['rep_factorLevel'] = [ {'factorLevel': 'Factor levels for each factor need to be unique'} for t in factor['rep_factorLevel'] ]
            break
        if overlapping_selection:
            error_map['rep_factorName'] = [ dict() for t in factors ]
            for i in range( len( factors ) ):
                error_map['rep_factorName'][i]['rep_factorLevel'] = [ {'factorIndex': 'The samples from different factors are not allowed to overlap'} for t in factor['rep_factorLevel'] ]
            break

    if factor_duplication:
        error_map['rep_factorName'] = [ dict() for t in factors ]
        for i in range( len( factors ) ):
            error_map['rep_factorName'][i]['factorName'] = 'Factor names need to be unique'

