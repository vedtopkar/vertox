import subprocess

def linearfold(test_sequence, linearfold_executable_path, vienna_mode=False):
    """
    Given a test_sequence:
    Run linearfold and return the predicted secondary structure and folding energy

    Can optionally run linearfold in Vienna mode (default is Contrafold mode)
    """
    if vienna_mode:
        vienna_flag = ' -V'
    else:
        vienna_flag = ''

    output = subprocess.check_output('echo "{}" | {}{}'.format(test_sequence,
                                                               linearfold_executable_path,
                                                               vienna_flag),
                                     shell=True, encoding='utf8')

    print(output)
    output = output.split('\n')[1]
    structure, energy = output.split(' ')
    energy = float(energy.replace('(', '').replace(')', ''))
    return structure, energy
