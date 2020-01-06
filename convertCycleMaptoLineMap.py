
import fileinput

def permToMap(perm, varList):
    '''
    Read a permutation from the perm.txt and convert it to a map
    :param perm: (X)(Y Z)(W)
    :param varList: ['X','Y','Z','W']
    :return:
        perm_map: {1: 1, 2: 2, 3: 4, 4: 3, 5: 5, 8: 8, 6: 7, 7: 6}
        length - 1: the length of the longest cycle
    '''
    perm_map = {}

    perm = perm.strip().split(')(') # result:  ['(X', 'Y, Z', 'W)']
    perm[0] = perm[0].lstrip('(')
    perm[-1] = perm[-1].rstrip(')') # result:  ['(X', 'Y, Z', 'W']
    # length = 1
    for subperm in perm: # example: 'Y, Z'
        elements = subperm.split(' ') # example: ['Y', 'Z']
        # length = lcm(length, len(elements))
        for index, element in enumerate(elements):
            perm_map[varList.index(element) + 1] = varList.index(elements[index + 1 if index <= len(elements) - 2 else 0]) + 1
    return perm_map  # , length - 1  # since the last map will map back to the 1st one


for line in fileinput.input('input caching 3,2.txt'):
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
        break

foutput = open('perm caching 3,2.txt', 'w')


with open('perm.txt') as filePerm:
    for perm in filePerm:  # for each permutation
        if perm == '\n':
            continue

        # for each permutation, convert it to map
        perm_map = permToMap(perm, varList)
        print(perm_map)
        for i in range(len(varList)):
            foutput.write(varList[perm_map[i+1]-1])
            foutput.write(' ')
        foutput.write('\n')
        print(perm_map)



