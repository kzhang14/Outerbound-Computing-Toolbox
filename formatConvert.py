#!/usr/bin/python

'''
Input: a .txt file containing the:
1st line: random raviables
    Can be the following symbols: alphanumeric (a-z, A-Z, 0-9) or one of these symbols: ! ' # $ % & / . ; ? @ _ ` ' { } ~.
2nd line to end: constraints

Output:
    CPLEX LP file format file
'''


import time
from imposeSym import *
from itertools import combinations
import fileinput
import numpy as np
import sys
from scipy.sparse import csc_matrix


def split(line):
    '''
    :param line: 'H(X1, X2|X3) - I(X1; W1 | X2) >= 0'
    :return: coef: [1, -1]
             entity: ['H(X1, X2|X3)', 'I(X1; W1 | W2)']
             sign: 'geq'
             constant: 0
    '''
    i = 0
    coef = []
    entity = []
    if line[0] == 'H' or line[0] == 'I' or line[0] == '_':
        line = '1' + line

    # seperate the inequality to single entities
    while (line[i] not in ['<', '>', '=']):
        j = i+1
        while line[j] != 'H' and line[j] != 'I' and line[j] != '_':
            j += 1

        coef_tmp = float(line[i:j].strip().replace(' ', '')) if line[i:j].strip() != '-' and line[i:j].strip() != '+' else float(line[i:j].strip() + '1')

        i = j
        while line[j] != ')' and line[j] != '@':
            j += 1
        entity_tmp = line[i:j+1].strip()

        # change it to canonical form:
        toCanonical(coef_tmp, entity_tmp, coef, entity)

        i = j + 1
        while line[i] == ' ': i += 1

    # reaches to the end
    if line[i+1] == '=':
        sign = line[i: i+2]
        i = i + 2
    else:
        sign = line[i: i+1]
        i = i + 1

    constant = float(line[i:].strip())
    return coef, entity, sign, constant


def splitObj(line):
    '''
    :param line: the objective function 'H(X1, X2|X3) - I(X1; W1 | X2)'
    :return: coef: [1, -1]
             entity: ['H(X1, X2|X3)', I(X1; W1 | W2)]
    '''

    line += '<0'
    coef, entity, _, _ = split(line)
    return coef, entity


def toCanonical(coef, entity, coef_vec, entity_vec):
    '''
    :param line: an information quantity, say: -2, I(X1; X2)
                since its canonical form is: I(X1; X2) = H(X1) + H(X2) - H(X1, X2)
    :return: coef_vec: [-2, -2, 2]
             entity_vec: ['H(X1), H(X2), H(X1, X2)']
    '''
    if entity[0] == 'H':
        is_conditional = entity.find('|')
        if is_conditional == -1:
            coef_vec.append(coef)
            entity_vec.append(entity)

        else:
            coef_vec += [coef, -coef]
            entity_vec.append(entity.replace('|', ','))
            entity_vec.append('H(' + entity[is_conditional + 1:])

    elif entity[0] == 'I':
        is_conditional = entity.find('|')
        if is_conditional == -1:
            coef_vec += [coef, coef, -coef]
            pos_semi_colon = entity.find(';')
            entity_vec.append('H' + entity[1 : pos_semi_colon] + ')')
            entity_vec.append('H('+ entity[pos_semi_colon + 1 :])
            entity_vec.append('H' + entity[1:].replace(';', ','))

        else:
            coef_vec += [coef, coef, -coef, -coef]
            pos_semi_colon = entity.find(';')
            entity_vec.append('H' + entity[1 : pos_semi_colon] + ',' + entity[is_conditional + 1: ])
            entity_vec.append('H('+ entity[pos_semi_colon + 1 :is_conditional].strip() + ',' + entity[is_conditional + 1: ])
            entity_vec.append('H('+ entity[is_conditional + 1:])
            entity_vec.append('H' + entity[1:].replace(';', ',').replace('|', ','))

def jointEtrptoPos(jointEtrp, varList):
    '''
    :param jointEtrp: a string like 'H(X,Y)'
    :param varList: a list of variables [X, Y, Z]
    :return: 3, since 011->3
    '''
    pos = 0
    for i, v in enumerate(varList):
        if v in jointEtrp:
            pos += 2 ** i
    return pos


def lexsort_based(data):
    sorted_data = data[np.lexsort(data.T), :]
    row_mask = np.append([True], np.any(np.diff(sorted_data, axis=0), 1))
    return sorted_data[row_mask]

def powerset(items):
    combo = [0] * (2 ** len(items) - 1)
    i = 0
    for r in range(1, len(items) + 1):
        #use a list to coerce a actual list from the combinations generator
        pset = list(combinations(items, r))
        combo[i: i + len(pset)] = [sum(x) for x in pset]
        i += len(pset)
    return combo

"""
parser = argparse.ArgumentParser()
parser.add_argument('problemFile',  help='Provide the problem file, e.g. input.txt')
parser.add_argument('symmetryFile', help='Provide the permutation group file, e.g. perm.txt')
parser.add_argument('noConstraints',  help='The number of problem specific constraints, e.g. 16')
parser.add_argument('outputFormat',  help='The output file format, e.g. MPS')
args = parser.parse_args()


i_file = sys.argv[1]
s_file = sys.argv[2]
num_of_prob_spec_constraints = 1 + int(sys.argv[3]) # To initialize the coefficient matrix size
o_format = sys.argv[4]
"""

i_file = 'input_caching_2_3.txt'
s_file = 'perm_caching_2_3.txt'
num_of_prob_spec_constraints = 1 + 37 # To initialize the coefficient matrix size
o_format = 'MPS'



''' ====================== Process 'input.txt' ====================================== '''
for line in fileinput.input(i_file):
    # escape the blank line
    if line == '\n':
        continue

    line = line.strip()

    if line.startswith('Variables'):
        lineState = 'V'
        continue
    elif line.startswith('Maximize'):
        lineState = 'O-max'
        continue
    elif line.startswith('Minimize'):
        lineState = 'O-min'
        continue
    elif line.startswith('Subject To'):
        lineState = 'S'
        continue


    if lineState == 'V':

        varList = [item.strip() for item in line.rstrip('\n').split(',')]
        num_of_var = len(varList)

        # Initialize the coefficient matrix
        MPScoefMatrix1_row = np.empty([1])
        MPScoefMatrix1_col = np.empty([1])
        MPScoefMatrix1_data = np.zeros([1])                            # Do not set to empty, it'll produce random number initially, we need to use +=
        MPSsignVector1 = ['None'] * num_of_prob_spec_constraints       # vector of string
        MPSbVector1 = [0] * num_of_prob_spec_constraints               # vector of float

        # currently we are at the beginning
        non_sparse_pos = -1
        current_constraint = -1

    elif lineState == 'O-max' or lineState == 'O-min':

        # print('Objective function is :', line)
        coef, entity = splitObj(line)
        current_constraint += 1

        c = [0] * (2 ** num_of_var - 1)

        for i in range(len(coef)):
            lp_var = entity[i]
            pos = jointEtrptoPos(lp_var, varList)
            c[pos - 1] += coef[i]

            non_sparse_pos += 1
            # Expand size to twice large
            if non_sparse_pos + 1 > len(MPScoefMatrix1_row):
                MPScoefMatrix1_row = np.append(MPScoefMatrix1_row, np.empty([len(MPScoefMatrix1_row)]))
                MPScoefMatrix1_col = np.append(MPScoefMatrix1_col, np.empty([len(MPScoefMatrix1_col)]))
                MPScoefMatrix1_data = np.append(MPScoefMatrix1_data, np.zeros([len(MPScoefMatrix1_data)]))

            MPScoefMatrix1_row[non_sparse_pos] = current_constraint
            MPScoefMatrix1_col[non_sparse_pos] = pos - 1
            MPScoefMatrix1_data[non_sparse_pos] = coef[i]

        MPSsignVector1[current_constraint] = 'NaN'
        MPSbVector1[current_constraint] = np.NaN

    elif lineState == 'S':

        # print('Information inequality: ', line)
        line = line.strip() # This is a constraint imposed by symmetry
        coef, entity, sign, constant = split(line)

        if sign == '=':

            current_constraint += 1
            for i in range(len(coef)):
                lp_var = entity[i]
                pos = jointEtrptoPos(lp_var, varList)

                non_sparse_pos += 1
                # Expand size to twice large
                if non_sparse_pos + 1 > len(MPScoefMatrix1_row):
                    MPScoefMatrix1_row = np.append(MPScoefMatrix1_row, np.empty([len(MPScoefMatrix1_row)]))
                    MPScoefMatrix1_col = np.append(MPScoefMatrix1_col, np.empty([len(MPScoefMatrix1_col)]))
                    MPScoefMatrix1_data = np.append(MPScoefMatrix1_data, np.zeros([len(MPScoefMatrix1_data)]))

                MPScoefMatrix1_row[non_sparse_pos] = current_constraint
                MPScoefMatrix1_col[non_sparse_pos] = pos - 1
                MPScoefMatrix1_data[non_sparse_pos] = coef[i]

            b_eq = constant
            MPSsignVector1[current_constraint] = sign
            MPSbVector1[current_constraint] = b_eq

        else:

            current_constraint += 1
            for i in range(len(coef)):
                lp_var = entity[i]
                pos = jointEtrptoPos(lp_var, varList)

                non_sparse_pos += 1
                # Expand size to twice large
                if non_sparse_pos + 1 > len(MPScoefMatrix1_row):
                    MPScoefMatrix1_row = np.append(MPScoefMatrix1_row, np.empty([len(MPScoefMatrix1_row)]))
                    MPScoefMatrix1_col = np.append(MPScoefMatrix1_col, np.empty([len(MPScoefMatrix1_col)]))
                    MPScoefMatrix1_data = np.append(MPScoefMatrix1_data, np.zeros([len(MPScoefMatrix1_data)]))

                MPScoefMatrix1_row[non_sparse_pos] = current_constraint
                MPScoefMatrix1_col[non_sparse_pos] = pos - 1
                MPScoefMatrix1_data[non_sparse_pos] = coef[i]


            b = constant

            MPSsignVector1[current_constraint] = sign
            MPSbVector1[current_constraint] = b

start_lp = time.time()

MPScoefMatrix1_row = MPScoefMatrix1_row[0 : non_sparse_pos + 1]
MPScoefMatrix1_col = MPScoefMatrix1_col[0 : non_sparse_pos + 1]
MPScoefMatrix1_data = MPScoefMatrix1_data[0 : non_sparse_pos + 1]

MPScoefMatrix1_sparse = csc_matrix((MPScoefMatrix1_data, (MPScoefMatrix1_row, MPScoefMatrix1_col)), shape=(num_of_prob_spec_constraints, 2**num_of_var - 1))

if o_format == 'LP':
    ''' ======================================== Write LP File ============================================== '''
    print('Write LP File')

    try:
       # open file stream
       fileoutputLP = open('outputLP.lp', 'w')
    except IOError:
       print('There was an error writing to outputLP.lp')
       sys.exit()


    fileoutputLP.write('Maximize\n' if lineState == 'O-max' else 'Minimize\n')
    fileoutputLP.write(' obj: ')

    # 1st line: objective function
    c = MPScoefMatrix1_sparse.getrow(0)
    fileoutputLP.write('\t' + str(c.data[0]) + ' x' + str(c.indices[0] + 1))
    for i in c.indices[1:]:
        if c[0, i] > 0:
            fileoutputLP.write(' + ' + str(c[0, i]) + ' x' + str(i + 1))
        else:
            fileoutputLP.write(' - ' + str(abs(c[0, i])) + ' x' + str(i + 1))
    fileoutputLP.write('\n')
    fileoutputLP.write('Subject To\n')

    # Starting from the 2nd line, it is the constraints
    for i in range(1, num_of_prob_spec_constraints):
        A = MPScoefMatrix1_sparse.getrow(i)
        sign = MPSsignVector1[i]
        b = MPSbVector1[i]

        fileoutputLP.write('\t' + str(A.data[0]) + ' x' + str(A.indices[0] + 1))
        for j in A.indices[1:]:
            if A[0, j] > 0:
                fileoutputLP.write(' + ' + str(A[0, j]) + ' x' + str(j + 1))
            else:
                fileoutputLP.write(' - ' + str(abs(A[0, j])) + ' x' + str(j + 1))
        fileoutputLP.write(' ' + sign + ' ' + str(b) + '\n')


    fileoutputLP.write('End')
    fileoutputLP.close()

    end_lp = time.time()
    print('Finish. Used {0:.4} seconds\n'.format(end_lp - start_lp))

elif o_format == 'MPS':

    ''' ====================== Write MPS file ======================'''
    start_mps = time.time()
    print('Write MPS File')

    try:
       # open file stream
       fileoutputMPS = open('outputMPS.mps', 'w')
    except IOError:
       print('There was an error writing to outputMPS.mps')
       sys.exit()


    fileoutputMPS.write('{:<14}{}\n'.format('NAME', 'ALPHA'))
    fileoutputMPS.write('ROWS\n')

    for i, item in enumerate(MPSsignVector1):
        if item == 'NaN':
            fileoutputMPS.write('{}{:<3}{:<10}\n'.format(' ', 'N', 'COST'))
        elif item == '<=':
            fileoutputMPS.write('{}{:<3}{:<10}\n'.format(' ', 'L', 'R' + str(i)))
        elif item == '>=':
            fileoutputMPS.write('{}{:<3}{:<10}\n'.format(' ', 'G', 'R' + str(i)))
        elif item == '=':
            fileoutputMPS.write('{}{:<3}{:<10}\n'.format(' ', 'E', 'R' + str(i)))

    fileoutputMPS.write('COLUMNS')

    for i, j in zip(MPScoefMatrix1_sparse.indptr, np.append(MPScoefMatrix1_sparse.indptr[1:], None)):
        nextRow = True
        for data in MPScoefMatrix1_sparse.data[i:j]:
            if nextRow == True:
                fileoutputMPS.write('\n{}{:<10}{:<15}{:<+10.4f}'.format(4*' ', 'C'+str(j), 'R' + str(i) if i != 0 else 'COST', data))
                nextRow = False
            elif nextRow == False:
                fileoutputMPS.write('{:<15}{:<+10.4f}'.format('R' + str(i) if i != 0 else 'COST', data))
                nextRow = True


    fileoutputMPS.write('\nRHS')
    nextRow = True
    for i in range(1, len(MPSbVector1)):
        if nextRow == True:
            fileoutputMPS.write('\n{}{:<10}{:<15}{:<+10.4f}'.format(4*' ', 'RHS1', 'R' + str(i), MPSbVector1[i]))
            nextRow = False
        elif nextRow == False:
            fileoutputMPS.write('{:<15}{:<+10.4f}'.format('R' + str(i), MPSbVector1[i]))
            nextRow = True

    fileoutputMPS.write('\nENDATA')

    end_mps = time.time()
    print('Finish. Use {0:.4} seconds\n'.format(end_mps - start_mps))




