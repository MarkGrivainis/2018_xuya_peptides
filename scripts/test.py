def validity_check(sequence):
    bases = 'AaCcGgTt'
    return all(char in bases for char in sequence)

print(validity_check('ACgTu'))
