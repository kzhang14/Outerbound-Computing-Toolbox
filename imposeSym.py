def binaryToDecimal(num_of_var):
    binaryToDecimal_map = {}
    decimalToBinary_map = {}

    for j in range(1, 2 ** num_of_var):
        binary = bin(j)[2:].zfill(num_of_var)[::-1]
        binaryToDecimal_map[int(binary)] = j
        decimalToBinary_map[j] = binary
    return binaryToDecimal_map, decimalToBinary_map

    # binaryToDecimal_map is a binary to integer mapping
    # {10000000: 1, 1000000: 2, 11000000: 3, 100000: 4, 10100000: 5, 1100000: 6, 11100000: 7, 10000: 8, etc.


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
    for index, var in enumerate(perm.strip().split()):
        perm_map[index + 1] = varList.index(var) + 1
    return perm_map

def mapping(map_from, map_to_binary, perm_map, num_of_var, binaryToDecimal_map, decimalToBinary_map):
    '''
        :param map_from: 1~15
        :return:  1~15
    '''

    map_from_binary = decimalToBinary_map[map_from]
    for j in range(1, num_of_var + 1):
        map_to_binary[perm_map[j] - 1] = map_from_binary[j - 1]
    tmp = ''.join(map_to_binary)
    map_to = binaryToDecimal_map[int(tmp)]
    return map_to


def genEquivalentClass(varList, num_of_var, binaryToDecimal_map, decimalToBinary_map, symmetry_file):
    perm_group = []
    with open(symmetry_file) as filePerm:
        for perm in filePerm:  # for each permutation
            if perm == '\n':
                continue

            # for each permutation, convert it to map
            perm_map = permToMap(perm, varList)
            perm_group.append(perm_map)

    all_entropy_label = [0] * (2 ** num_of_var - 1)
    jointentropy_map = dict()
    numofVar_reduce = 0
    reduced_vars = dict()

    i = 0
    benchmark = int((2 ** num_of_var - 1) / 20)
    while i != -1:
        if i % benchmark < 20:
             print("{0:.2%}".format(i / (2**num_of_var - 1)))

        leader = i + 1

        all_entropy_label[leader - 1] = 1
        jointentropy_map[leader] = leader
        numofVar_reduce += 1
        reduced_vars[leader] = numofVar_reduce - 1
        test = 0
        for perm_map in perm_group:
            test += 1

            cycle = [leader]

            while 1:
                map_from = cycle[-1]
                map_to_binary = ['0'] * num_of_var
                map_to = mapping(map_from, map_to_binary, perm_map, num_of_var, binaryToDecimal_map, decimalToBinary_map)
                if map_to in cycle:
                    break
                cycle.append(map_to)
                jointentropy_map[map_to] = leader
                all_entropy_label[map_to - 1] = 1
        try:
            i = all_entropy_label.index(0)
        except ValueError:
            i = -1
    return jointentropy_map, numofVar_reduce, reduced_vars
