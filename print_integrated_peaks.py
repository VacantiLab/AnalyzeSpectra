def print_integrated_peaks(file_directory,files,fragment_list,file_data):

    import pdb

    output_text_file = file_directory + 'integrated_peaks.txt'
    file_object_text = open(output_text_file,'w')

    files_string = '\t'.join(files)
    file_object_text.write(' \t' + files_string)
    file_object_text.write('\n')

    file_object_text.write('Retention Times\n')
    for fragment in fragment_list:
        file_object_text.write(fragment)
        file_object_text.write('\t')
        for filename in files:
            rt = str(file_data[filename]['fragments'][fragment]['rt'][0])
            file_object_text.write(rt)
            file_object_text.write('\t')
        file_object_text.write('\n')

    file_object_text.write('\n')
    file_object_text.write('Peak Areas\n')
    for fragment in fragment_list:
        file_object_text.write(fragment)
        file_object_text.write('\t')
        for filename in files:
            rt = str(file_data[filename]['fragments'][fragment]['tot_area'])
            file_object_text.write(rt)
            file_object_text.write('\t')
        file_object_text.write('\n')

    file_object_text.write('\n')
    file_object_text.write('MIDs\n')
    for fragment in fragment_list:
        mid_length = len(file_data[filename]['fragments'][fragment]['mid_c'])
        mid_members = range(0,mid_length)
        for M in mid_members:
            fragment_mi_name = fragment + ' ' + 'M' + str(mid_members[M])
            file_object_text.write(fragment_mi_name)
            file_object_text.write('\t')
            for filename in files:
                current_mi = str(file_data[filename]['fragments'][fragment]['mid_c'][M])
                file_object_text.write(current_mi)
                file_object_text.write('\t')
            file_object_text.write('\n')
        file_object_text.write('\n')


    file_object_text.close()
    return()
