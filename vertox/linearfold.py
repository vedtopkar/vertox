import subprocess
from .path import PATH_LINEARFOLD_EXECUTABLE

def linearfold(test_sequence, vienna_mode=False):
    """
    Given a test_sequence:
    Run linearfold and return the predicted secondary structure and folding energy

    Can optionally run linearfold in Vienna mode (default is Contrafold mode)
    """
    if vienna_mode:
        vienna_flag = ' -V'
    else:
        vienna_flag = ''

    command = 'echo "{}" | {}{}'.format(test_sequence, PATH_LINEARFOLD_EXECUTABLE, vienna_flag)
    output = subprocess.check_output(command, shell=True, encoding='utf8')

    output = output.split('\n')[1]
    structure, energy = output.split(' ')
    energy = float(energy.replace('(', '').replace(')', ''))
    return structure, energy
