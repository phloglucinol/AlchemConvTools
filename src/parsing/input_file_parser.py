import types

'''
###Input file format:
[Basic_settings]
simulation_software = openmm
file_directory = .
annihilation_group_num = None
file_prefix = 'state_'
file_suffix = '.csv'
fraction = 1.0
output_csv_filename = 'free_ene.csv'
energy_unit = 'kcal/mol'
calculation_windows = 'all'
std_mode = bar_std
plot_du = False

[Convergence_analysis]
Time_series = True
Boltzmann_weight_PdU = True
Curve_Fitting_Method = True
Bennett_Overlap_Value = True
Pai_Value = True
MeasurementOfConvergence_Value = True
Reweighting_SelfConsistent = True

[Reweighting_analysis]
reweight_use_wins_step = 0.1
reweight_use_wins_num = 1
ifforward = True
ifdiagnum = True
ifplotheatmap = True
ifreweight_time_serial = True
'''

class InputParser:
    def __init__(self, filename):
        # Define the expected keys for each section and their default values
        self.def_sections()

        # Initialize variables to hold the parsed data
        self.data = {}
        for section, options in self.sections.items():
            self.data[section] = {}
            for option, default_value in options.items():
                self.data[section][option] = default_value

        self.current_section = None

        # Parse the input file
        with open(filename, 'r') as f:
            for line in f:
                # Remove inline comments
                line = line.split('#')[0].strip()
                if not line:  # Skip empty lines and comments
                    continue
                if line.startswith('[') and line.endswith(']'):
                    self.current_section = line[1:-1]
                elif self.current_section is not None:
                    for item in line.split(','):
                        key, value = item.split('=')
                        key = key.strip()
                        value = value.strip()
                        if key in self.sections[self.current_section]:
                            if value.isdigit():
                                self.data[self.current_section][key] = int(value)
                            elif value.replace('.', '').isdigit():
                                self.data[self.current_section][key] = float(value)
                            elif value.lower() == 'true':
                                self.data[self.current_section][key] = True
                            elif value.lower() == 'false':
                                self.data[self.current_section][key] = False
                            elif value.lower() == 'none':
                                self.data[self.current_section][key] = None
                            else:
                                self.data[self.current_section][key] = value.strip("'").strip('"')


        # Define the section-related methods dynamically
        for section in self.sections:
            def get_section_data(self, section=section):
                return self.data.get(section, {})
            setattr(self, f'get_{section}', types.MethodType(get_section_data, self))

    def def_sections(self):
        # Define the expected keys for each section and their default values
        self.sections = {
            'Basic_settings': {
                'simulation_software' : 'openmm', # openmm, gromacs, amber
                'file_directory': '.', 
                'file_prefix': 'state_',
                'file_suffix': '.csv',
                'fraction': 1.0,
                'output_csv_filename': 'free_ene.csv',
                'energy_unit': 'kcal/mol',
                'calculation_windows': 'all',
                'subsample': False,
                'std_mode': 'bar_std',
                'plot_du': True,
                },
            'Convergence_analysis': {
                'Time_series': True,
                'Boltzmann_weight_PdU': True,
                'Curve_Fitting_Method': True,
                'Bennett_Overlap_Value': True,
                'Pai_Value': True,
                'MeasurementOfConvergence_Value': True,
                'Reweighting_SelfConsistent': True,
                },
            'Reweighting_analysis': {
                'reweight_use_wins_step': 0.1,
                'reweight_use_wins_num': 1,
                'ifforward': True,
                'ifdiagnum': True,
                'ifplotheatmap': True,
                'ifreweight_time_serial': True,
                'reweight_err_bar_max': 2.0,
                },
            'Time_series_analysis':{
                'useFEP_check_everywins_Time_series': False,
                'plot_plan': 'forward|reverse|moving',
                'std_step_ratio': 0.1,
                'divided_ratio': 0.01,
                'dGdlambda_divided_ratio': 0.1,
                'width_ratio': 0.2,
                'ifplot_bar_time_series': True,
                'ifplot_bar_dG_dlambda':False,
                'ifplot_fep_dG_dlambda':False,
                'iflog_fep_csv': True,
                'iflog_bar_csv': True,
                'ifplot_fep_time_series': True,
                'use_forward': True,
            }
        }



if __name__ == '__main__':
    def main():
        filename = 'input.txt'
        parser = InputParser(filename)
        basic_settings = parser.get_Basic_settings()
        print(basic_settings)
        convergence_analysis_settings = parser.get_Convergence_analysis()
        print(convergence_analysis_settings)
        if convergence_analysis_settings['Reweighting_SelfConsistent'] == True:
            reweighting_analysis_settings = parser.get_Reweighting_analysis()
            print(reweighting_analysis_settings)

    main()

