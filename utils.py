import os
import gzip


def handle_gzip(look_up_directory):
    files = [f for f in os.listdir(look_up_directory) if
             os.path.isfile(os.path.join(look_up_directory, f)) and f.endswith(".gz")]
    files_sorted_by_time = sorted(files, key=lambda x: os.path.getmtime(os.path.join(look_up_directory, x)),
                                  reverse=True)
    file = files_sorted_by_time[0]

    # dispatch file to cwd
    input_file_path = os.path.join(look_up_directory, file)
    os.rename(input_file_path, input_file_path.replace(" ", ""))
    input_file_path = input_file_path.replace(" ", "")
    output_filename = ".".join(file.split(r"\ ")[-1].split(".")[:-1])

    # dump content
    with gzip.open(input_file_path, 'rb') as f_in:
        with open(output_filename, 'wb') as f_out:
            # Read and decompress the content from the input .gz file and write it to the output file
            f_out.write(f_in.read())
    return output_filename
