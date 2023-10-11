import pandas as pd

# Load the data into a DataFrame
data = pd.read_csv('chrom_names', sep=' ', names=['file_name', 'chromosome_name'])

# Define a function to check if a file contains a string
def file_contains(file_name, string):
    file_name = f'2bit_files/{file_name}/chrom.sizes'
    try:
        with open(file_name, 'r') as f:
            return string in f.read()
    except FileNotFoundError:
        print(f'File not found: {file_name}')
        

# Create a new column to store the result
data['contains'] = data.apply(lambda row: file_contains(row['file_name'], row['chromosome_name']), axis=1)

# Save the result to a text file
data.to_csv('check_chromnames_result.txt', sep=' ', index=False)


'''
import pandas as pd

# Load the data into a DataFrame
data = pd.read_csv('chrom_names1', sep=' ', names=['file_name', 'chromosome_name'])

# Define a function to check if a file contains a string
def file_contains(file_name, string):
    file_name = f'2bit_files/{file_name}/chrom.sizes'
    try:
        with open(file_name, 'r') as f:
            for line in f:
                if line.strip() == string:
                    return True
            return False
    except FileNotFoundError:
        print(f'File not found: {file_name}')
        return False

# Create a new column to store the result
data['contains'] = data.apply(lambda row: file_contains(row['file_name'], row['chromosome_name']), axis=1)

# Save the result to a text file
data.to_csv('result.txt', sep=' ', index=False)

'''