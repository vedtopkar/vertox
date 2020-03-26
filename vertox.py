import configparser
import sys
import os

import pandas as pd

from vertox.mutate_rescue import generate_mutate_rescue_library

def parse_config(config_filepath):
    config = configparser.ConfigParser()
    config.read(config_filepath)
    master_sequence = config['Data Input']['master_sequence']
    master_structure = config['Data Input']['master_structure']
    regions = config['Data Input']['regions'].split('\n')

    parsed_regions = []
    for region in regions:
        region = region.split(':')
        assert len(region) == 2

        parsed_regions.append([int(region[0]), int(region[1])])

    data_input = {'master_sequence': master_sequence,
                 'master_structure': master_structure,
                 'regions': parsed_regions}

    library_parameters = {}
    library_parameters['total_iterations'] = int(config['Library Parameters']['total_iterations'])
    library_parameters['min_helix_size'] = int(config['Library Parameters']['min_helix_size'])
    library_parameters['n_stochastic_results'] = int(config['Library Parameters']['n_stochastic_results'])
    library_parameters['five_prime_adapter'] = config['Library Parameters']['five_prime_adapter']
    library_parameters['three_prime_adapter'] = config['Library Parameters']['three_prime_adapter']

    io_parameters = config['IO Parameters']

    return data_input, library_parameters, io_parameters


def main(config_filepath):
    data_input, library_parameters, io_parameters = parse_config(config_filepath)

    master_sequence = data_input['master_sequence']
    master_structure = data_input['master_structure']
    regions = data_input['regions']

    output_folder = io_parameters['output_folder']
    experiment_name = io_parameters['experiment_name']

    sequences = []
    structures = []

    for region in regions:
        sequence = master_sequence[region[0]:region[1] + 1]
        structure = master_structure[region[0]:region[1] + 1]

        sequences.append(sequence)
        structures.append(structure)

    print('Generating libraries for:')
    print('Regions: {}'.format(regions))
    print('Sequences: {}'.format(sequences))
    print('Structures: {}'.format(structures))

    all_results = {}
    for region, sequence, structure in zip(regions, sequences, structures):
        all_results[str(region)] = generate_mutate_rescue_library(sequence,
                                                                   structure,
                                                                   total_iterations=library_parameters['total_iterations'],
                                                                   min_helix_size=library_parameters['min_helix_size'],
                                                                   stochastic_results=library_parameters['n_stochastic_results'])

    print(all_results)

    # Open an xls writer and write each region to a different datasheet
    writer = pd.ExcelWriter(os.path.join(output_folder, experiment_name + '_mutate_rescue_library.xlsx'), engine='xlsxwriter')
    for key in all_results.keys():
        all_results[key].to_excel(writer, sheet_name=key.replace('[', '').replace(']', '').replace(', ', '..'), index=False)

    master_results_sheet = pd.concat(list(all_results.values()))
    master_results_sheet.to_excel(writer, sheet_name='All Results', index=False)

    writer.save()

    # Create a separate file that can be easily uploaded to IDT for opool downloading
    ordering_sheet = master_results_sheet[['Sequence']]

    # Add adapters
    ordering_sheet['Sequence'] = library_parameters['five_prime_adapter'] + \
                                 ordering_sheet['Sequence'].astype(str) + \
                                 library_parameters['three_prime_adapter']

    # Make sure everything is uppercase and that there are Ts not Us
    ordering_sheet['Sequence'] = ordering_sheet['Sequence'].str.upper().replace('U', 'T')

    # Add pool name as experiment name and export
    ordering_sheet['Pool name'] = experiment_name
    ordering_sheet.to_excel(os.path.join(output_folder, experiment_name + '_mutate_rescue_OPOOL_ORDERING.xlsx'), index=False)


if __name__ == '__main__':
    main(sys.argv[1])

