#!/usr/bin/env python

# Written by Karolis Uziela in 2014

import sys

################################ Usage ################################

script_name = sys.argv[0]

usage_string = """
Description:

The script checks if the psiblast matrix file (.psi) does not contain all-zero matrix. If the matrix is all-zero, it is replaced by the matrix
that would correspond to the input sequence itself (100%% probability of obtaining amino acid that is in the actual sequence). If the input psiblast
matrix file does not contain all-zero matrix, the output matrix will be the same as the input matrix.

Usage:

    %s [Parameters]
    
Parameters:

    <input-file> - input psiblast matrix file (.psi)
    <sequence-fasta> - input sequence file in fasta format. Note that the sequence has to be in one line (not text wrapped)
    <output-file> - output psiblast matrix file
""" % script_name

if len(sys.argv) != 4:
    sys.exit(usage_string)

################################ Functions ################################

def read_fasta(filename):
    f = open(filename)
    f.readline()
    n = 0
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        my_seq = line.rstrip('\n')
        n += 1
    f.close()
    if n > 1:
        sys.exit("ERROR: fasta file contains more than 2 lines. Note that input sequence should not be text-wrapped. Fasta file should contain 1 line for fasta ID and one line for sequence. ")
    if n == 0:
        sys.exit("ERROR: fasta file contains only 1 line. Fasta file should contain 1 line for fasta ID and one line for sequence.")
    return my_seq

def check_psiblast_matrix(input_matrix_file, output_matrix_file, fasta_seq):
    
    f = open(input_matrix_file)
    
    line_begin = []
    line_end = []
    matrix = []
    validate_seq = ""
    file_begin = ""
    file_end = ""
    file_middle = ""
    n = 0
    while True:
        line = f.readline()
        if len(line) == 0: 
            break
        n += 1
        #line = line.rstrip('\n')
        bits = line.split()
        #print len(bits)
        if len(bits) == 44:
            line_begin.append(line[:70])
            line_end.append(line[150:])
            matrix.append(map(int,bits[22:42]))
            validate_seq += bits[1]
            file_middle += line
        elif n <= 3:
            file_begin += line
        else:
            file_end += line

    if fasta_seq != validate_seq:
        sys.exit("ERROR: Sequence from fasta file does not match sequence from psiblast matrix file. \n Fasta sequence: \n %s\n Psiblast matrix sequence: \n %s" % (fasta_seq, validate_seq) )
    
    matrix_sum = sum(sum(x) for x in matrix)
    #print "Sum:"
    #print matrix_sum
    out_str = ""
    if matrix_sum == 0:
        new_matrix = create_new_matrix(fasta_seq)
        matrix_string = format_matrix_string(new_matrix, line_begin, line_end)
        out_str = file_begin + matrix_string + file_end
    else:
        out_str = file_begin + file_middle + file_end
    
    write_data(output_matrix_file, out_str)

    f.close()
    
def create_new_matrix(fasta_seq):
    amino_dict = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, "S": 15, "T": 16, "W": 17, "Y": 18, "V":19}
    new_matrix = []
    for letter in fasta_seq:
        freq_list = [0] * 20
        index = amino_dict[letter]
        freq_list[index] = 100
        new_matrix.append(freq_list)
    return new_matrix

def format_matrix_string(new_matrix, line_begin, line_end):
    matrix_string = ""
    for i in range(0,len(new_matrix)):
        freq_list = new_matrix[i]
        matrix_string += line_begin[i]
        for freq in freq_list:
            matrix_string += "%4d" % freq
        matrix_string += line_end[i]
    return matrix_string

def write_data(output_file, out_str):
    f = open(output_file,"w")  
    f.write(out_str)    
    f.close()


################################ Variables ################################

# Input files/directories
input_matrix_file = sys.argv[1]
input_seq_file = sys.argv[2]

# Output files/directories
output_matrix_file = sys.argv[3]

# Constants
# N/A

# Global variables
# N/A

################################ Main script ################################
    
sys.stderr.write("%s is running with arguments: %s\n" % (script_name, str(sys.argv[1:])))

my_seq = read_fasta(input_seq_file)

check_psiblast_matrix(input_matrix_file, output_matrix_file, my_seq)

sys.stderr.write("%s done.\n" % script_name)



