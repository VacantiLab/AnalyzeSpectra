def create_output_directory():
    import os #allows to get working directory
    import re #allows regular expression use

    current_directory = os.getcwd()
    output_directory = re.sub('/(?!.*/).*$','/output',current_directory)
        # '/(?!.*/)' -> '/' not followed by [(?!)] (any number of characters [.] any number of times [*] followed by '/')
        # This followed by any character [.] any number of times [*] followed by the end of the string [$]
        # Replaces the last file name with 'output'
    output_plot_directory = output_directory + '/' + 'output_plot.html'
    return(output_plot_directory,output_directory)
