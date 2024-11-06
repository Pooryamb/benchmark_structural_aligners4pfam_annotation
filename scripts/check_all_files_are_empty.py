import os
import argparse



def check_and_cleanup(directory, output_path):
    non_empty_files = []
    try:
        # Ensure the directory exists
        if not os.path.isdir(directory):
            print("Directory does not exist.")
            return

        # List all files in the directory
        files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

        # Check if all files are empty
        all_empty = True
        for file in files:
            if os.path.getsize(file) > 0:
                non_empty_files.append(file)
                all_empty = False

        # Act based on whether all files are empty
        if all_empty:
            print("All files are empty. Removing all files.")
            for file in files:
                os.remove(file)
        else:
            print("There are non-empty files in the directory.")
            with open(output_path, 'w') as outf:
                outf.write("\n".join(non_empty_files))

    except Exception as e:
        print(f"An error occurred: {e}")




if __name__=="__main__":
    parser = argparse.ArgumentParser(description = """This script will check if all the files within a directory are empty. If so, it will remove them all. 
                                                                               Otherwise, it will write the name of nonempty files into another file""")
    parser.add_argument("--dir", type=str, required=True, help="The path to the directory to check its files")
    parser.add_argument("--output", type=str, required=True, help="The path to the place to store the list of non-empty files")

    args = parser.parse_args()
    directory = args.dir
    output_path = args.output

    check_and_cleanup(directory, output_path)
