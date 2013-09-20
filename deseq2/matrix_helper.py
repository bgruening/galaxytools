def get_matrix_header( input_dataset ):
     input_handle = open( input_dataset.file_name )
     first_header = input_handle.readline()
     second_header = input_handle.readline()
     return [('%s::%s' % (cname2,cname1), str(int(col_num) + 1), False) for col_num, (cname2, cname1) in enumerate(zip(second_header.split()[1:],first_header.split()[1:])) ]

#print get_matrix_header('/home/videmp/projects/galaxy/galaxy-central/tools/deseq2/countmatrix.txt')
