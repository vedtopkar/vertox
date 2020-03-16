def linearfold(test_sequence):
    output = subprocess.check_output('echo "{}" | /usr/local/src/LinearFold/linearfold'.format(test_sequence), shell=True)
    output = output.split('\n')[1]
    structure, energy = output.split(' ')
    energy = float(energy.replace('(', '').replace(')', ''))
    return structure, energy
