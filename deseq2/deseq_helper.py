
from galaxy.tools.parameters import DataToolParameter

def get_matrix_header( input_dataset ):
    """
        Not used currently, because the reload of the ckeckboxes did not work.
    """
    input_handle = open( input_dataset.file_name )
    first_header = input_handle.readline()
    second_header = input_handle.readline()
    return [('%s::%s' % (cname2,cname1), str(int(col_num) + 1), False) for col_num, (cname2, cname1) in enumerate(zip(second_header.split()[1:],first_header.split()[1:])) ]



def _construct_error_map( error_map, rep_dict, rep_parent, child, error_value ):
    """
        Its no so easy to create a propper error_map for repetitions in Galaxy.
        This is a helper function.
    """

    error_map[ rep_parent ] = [ dict() for t in rep_dict ]
    for i in range( len( rep_dict ) ):
        error_map[ rep_parent ][i][ child ] = error_value



def validate_input( trans, error_map, param_values, page_param_map ):
    """
        Validates the user input, before execution.
    """
    factors = param_values['rep_factorName']
    factor_name_list = []
    factor_duplication = False
    level_duplication = False
    overlapping_selection = False

    first_condition = True
    factor_indieces = list()

    for factor in factors:
        # factor names should be unique
        fn = factor['factorName']
        if fn in factor_name_list:
            factor_duplication = True
            break
        factor_name_list.append( fn )

        level_name_list = list()
        factor_index_list = list()

        if first_condition and len( factor['rep_factorLevel'] ) < 2:
            # first condition needs to have at least 2 levels
            _construct_error_map( error_map, factors, 'rep_factorName', 'rep_factorLevel', [ {'factorLevel': 'The first condition should have at least 2 factor'} for t in factor['rep_factorLevel'] ] )

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

            print set(factor_index_list)
            print factor_indieces
            if set(factor_index_list) in factor_indieces:
                _construct_error_map( error_map, factors, 'rep_factorName', 'rep_factorLevel', [ {'factorLevel': 'It is not allowed to have two identical factors, that means two factors with the same toggeled checked boxes. '} for t in factor['rep_factorLevel'] ] )
            else:
                factor_indieces.append( set(factor_index_list) )



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

        first_condition = False

    if factor_duplication:
        _construct_error_map( error_map, factors, 'rep_factorName', 'factorName', 'Factor names need to be unique' )
        """
        error_map['rep_factorName'] = [ dict() for t in factors ]
        for i in range( len( factors ) ):
            error_map['rep_factorName'][i]['factorName'] = 'Factor names need to be unique'
        """
