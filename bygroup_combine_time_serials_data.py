import pandas as pd
from src.plotting.plotting_tools import PLOTTING
import numpy as np
import os
import math
from glob import glob
from optparse import OptionParser
import copy
import re

def extract_tuple_and_substract(string):
    # pattern = r'\((\d+\.?\d*),\s*(\d+\.?\d*),\s*(\d+\.?\d*)\)'
    pattern = r'\(([\d.]+(?:,\s*[\d.]+)*)\)\s*to\s*\(([\d.]+(?:,\s*[\d.]+)*)\)'
    matches = re.findall(pattern, string)
    lambda_float_diff_eq1_num = 0
    # print(matches)
    if len(matches) == 1:
        # [('0.0, 1.0, 0.2', '0.0, 1.0, 0.25')]
        tuple1 = tuple(float(x) for x in matches[0][0].split(','))
        tuple2 = tuple(float(x) for x in matches[0][1].split(','))
        subtracted_tuple = tuple(abs(a - b) for a, b in zip(tuple1, tuple2))
        for lambda_float_diff in subtracted_tuple:
            if lambda_float_diff == 1:
                lambda_float_diff_eq1_num+=1
    else:
        raise ValueError(f"Cannot find the two lambda tuple in the string: {string}")
    return lambda_float_diff_eq1_num

def get_newest_file(pattern):
    # Use glob to find all files matching the pattern
    files = glob(pattern)

    # Check if at least one file was found
    if files:
        # Find the file with the latest modification time
        # latest_file = max(files, key=os.path.getmtime)
        # Find the file with the largest lambda_float_diff_eq1_num
        latest_file = max(files, key=extract_tuple_and_substract)

        # Get the absolute path of the latest file
        latest_file_path = os.path.abspath(latest_file)

        # Print or return the path
        # print(latest_file_path)
    else:
        print("No files found matching the pattern.")
        latest_file_path = None
    return latest_file_path

class CombineMeta(type):
    def __new__(cls, name, bases, attrs):
        if 'combine_way' not in attrs:
            raise AttributeError(f"{name} class must define a combine_way method.")
        return super().__new__(cls, name, bases, attrs)


class DataFrameCombiner(metaclass=CombineMeta):
    def __init__(self, df_name_dir, delimiter=','):
        self.df_name_dir = df_name_dir
        self.df_dir = {}
        for key, value in df_name_dir.items():
            # print(value)
            self.df_dir[key] = pd.read_csv(value, delimiter=delimiter)
        self.total_df = None

    def combine_way(self):
        pass

         
class SumTimeFEDataFrameCombiner(DataFrameCombiner):
    def combine_way(self):
        total_df = None
        one_df = list(self.df_dir.values())[0]
        for key, df in self.df_dir.items():
            if total_df is None:
                total_df = pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
            else:
                total_df += pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
        # print(total_df)
        total_df.index = one_df.loc[:, 'time_ratio']
        self.total_df = total_df
    
    def cal_rolling_std(self, std_step_ratio):
        data_df = self.total_df
        std_step = int(np.floor(len(data_df)*std_step_ratio))
        std_columns = ['Rolling_STD(kcal/mol)',]
        std_df = data_df.rolling(std_step).std()
        for i in range(std_step):
            std_df.iloc[i,0]=std_df.iloc[std_step, 0]
        std_df.columns = std_columns
        total_df_with_std = pd.concat([data_df, std_df], axis=1)
        total_df_with_std.index = data_df.index
        return total_df_with_std

class SumFEDataFrameCombiner(DataFrameCombiner):
    def combine_way(self, std_column):
        fe_df = None
        std_df = None
        index = None
        for key, df in self.df_dir.items():
            if index is None:
                index = df.loc[:, 'delta_A_what_to_what']
                # print(index)
            else:
                assert df['delta_A_what_to_what'].equals(index), "DataFrames column 'delta_A_what_to_what' does not equal to each other."
            if fe_df is None:
                fe_df = pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
            else:
                fe_df += pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
            if std_df is None:
                std_df = pd.DataFrame(df.loc[:, std_column])
            else:
                std_df = (std_df**2+pd.DataFrame(df.loc[:, std_column])**2)**(0.5)
        total_df_with_std = pd.concat([fe_df, std_df], axis=1)
        total_df_with_std.index = index
        return total_df_with_std

class FreeEnergyDataFrameCombiner():
    def __init__(self, df_basic_name, sub_string_dir):
        self.df_basic_name = df_basic_name # free_ene.csv
        # self.base_path = base_path # /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out
        self.sub_string_dir = sub_string_dir # {'group1': 'complex/group1', 'group2': 'complex/group2', 'group3': 'complex/group3', 'group4': 'complex/group4'}
        self.df_name_dir = self.gen_dir_name_dir(sub_string_dir)
        self.fe = 0
        self.fe_std = 0
        for df_name in self.df_name_dir.values():
            if not os.path.exists(df_name):
                raise FileNotFoundError(f'{df_name} does not exist!')
            else:
                fe, std = self.get_one_df_fe_std(df_name)
                self.fe += fe
                self.fe_std = (self.fe_std**2 + std**2)**(0.5)

    def get_one_df_fe_std(self, df_name):
        df = pd.read_csv(df_name, delimiter='|')
        fe = float(df.iloc[-1, 1])
        std = float(df.iloc[-1, 2])
        return fe, std

    def gen_dir_name_dir(self, sub_string_dir, ):
        df_name_dir = {}
        for key, value in sub_string_dir.items():
            prefix_path = value
            df_name_dir[key] = f'{prefix_path}/sample_csv_data/fe_cal_out/{self.df_basic_name}'
        return df_name_dir

class TimeSeriesDataFrameCombiner():  
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_reverse_esti.csv
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_moving_esti.csv
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_forward_esti.csv
    # base_path: /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out/time_serial_check
    
    def __init__(self, df_basic_pattern, base_path, plotting_plan=['moving','forward', 'reverse']):
        self.df_basic_pattern = df_basic_pattern # _bar_
        self.base_path = base_path # /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/complex
        self.sub_string_dir = self.get_side_thermo_path(self.base_path, ['restraints*', 'electrostatics*', 'sterics*', 'release*']) 
        self.plot_fe_he_std_dir = {}
        self.FE_combine_obj = FreeEnergyDataFrameCombiner('free_ene.csv', self.sub_string_dir)
        
        for aly in plotting_plan:
            df_name_dir = self.gen_dir_name_dir(self.sub_string_dir, aly)
            combine_obj = SumTimeFEDataFrameCombiner(df_name_dir) # different combiner
            combine_obj.combine_way()
            self.plot_fe_he_std_dir[aly] = combine_obj.cal_rolling_std(0.1)
        self.plot_fe_he_std_dir['fe'] = self.FE_combine_obj.fe
        self.plot_fe_he_std_dir['fe_std'] = self.FE_combine_obj.fe_std
    
    def gen_dir_name_dir(self, sub_string_dir, type_):
        # setattr(self, f'{type_}_df_name', self.df_basic_name.replace('SUBSTI', type_))
        df_name_dir = {}
        for key, value in sub_string_dir.items():
            prefix_path = value # /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/complex/electrostatics
            df_name = get_newest_file(os.path.join(prefix_path, 'sample_csv_data', 'fe_cal_out', 'time_serial_check', f'{self.df_basic_pattern}*_{type_}_esti.csv'))
            df_name_dir[key] = df_name
        return df_name_dir
    
    def plotting_combined_time_serial(self, plot_fe_he_std_dir, png_name):
        plot_obj = PLOTTING()
        plot_obj.plot_fe_time_serial(png_name, **plot_fe_he_std_dir)

    def get_side_thermo_path(self, base_path, glob_patterns = ['restraints*', 'electrostatics*', 'sterics*']):
        dirs_lst = []
        for glob_pattern in glob_patterns:
            dirs_lst.extend(glob(os.path.join(base_path, glob_pattern)))
        dir_dict = {}
        for dir_ in dirs_lst:
            dir_dict[os.path.basename(dir_)] = dir_
        return dir_dict

def get_aly_restr_lig(res_info_csv):
    K = 8.314472*0.001  # Gas constant in kJ/mol/K
    V = 1.66            # standard volume in nm^3
    T = 300.0           # Temperature in Kelvin
    res_df = pd.read_csv(res_info_csv, )
    r0 = float(res_df.iloc[0, 4])/10 # distance in A
    K_r = 4184.0 # force constant for distance (kJ/mol/nm^2)

    thA = float(res_df.iloc[0, 5]) # Angle in rad
    K_thA = 41.84 # force constant for angle (kJ/mol/rad^2)

    thB = float(res_df.iloc[0, 6]) # Angle in rad
    K_thB = 41.84 # force constant for angle (kJ/mol/rad^2)

    K_phiA = 41.84 # force constant for angle (kJ/mol/rad^2)
    K_phiB = 41.84 # force constant for angle (kJ/mol/rad^2)
    K_phiC = 41.84 # force constant for angle (kJ/mol/rad^2)

    arg =(
        (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(thA) * math.sin(thB))
        *
        (
            ( (K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC)**0.5 ) / ( (2.0 * math.pi * K * T)**(3.0) )
        )
    )
    dG = - K * T * math.log(arg)
    return abs(dG)/4.184

class ABFE_whole_mol_time_series_combiner():
    
    def __init__(self, res_aly, com_plot_fe_he_std_dir, lig_plot_fe_he_std_dir, png_name, plotting_plan=['moving','forward', 'reverse']):
        self.plot_fe_he_std_dir = {}
        com_fe = com_plot_fe_he_std_dir['fe']
        com_fe_std = com_plot_fe_he_std_dir['fe_std']
        lig_fe = lig_plot_fe_he_std_dir['fe']
        lig_fe_std = lig_plot_fe_he_std_dir['fe_std']
        whole_mol_fe = res_aly + lig_fe - com_fe
        whole_mol_fe_std = (com_fe_std**2 + lig_fe_std**2)**(0.5)
        for aly in plotting_plan:
            com_df_aly = com_plot_fe_he_std_dir[aly]
            lig_df_aly = lig_plot_fe_he_std_dir[aly]
            com_fe_aly = com_df_aly['free_energy(kcal/mol)']
            lig_fe_aly = lig_df_aly['free_energy(kcal/mol)']
            # print(com_fe_aly)
            # print(lig_fe_aly)
            whole_mol_df_aly = lig_fe_aly + res_aly - com_fe_aly
            whole_mol_df_aly.index = com_df_aly.index # assume the index of com_df_aly is the same as lig_df_aly
            # print()
            # print(whole_mol_df_aly)
            std_step_ratio = 0.1
            std_step = int(np.floor(len(whole_mol_df_aly)*std_step_ratio))
            std_columns = ['Rolling_STD(kcal/mol)',]
            std_df = whole_mol_df_aly.rolling(std_step).std()
            # print(std_df)
            # print(std_step)
            ori_std_df = copy.deepcopy(std_df)
            for i in range(std_step):
                # print(i, std_step, ori_std_df.iloc[std_step])
                std_df.iloc[i]=ori_std_df.iloc[std_step]
            std_df.columns = std_columns
            total_df_with_std = pd.concat([whole_mol_df_aly, std_df], axis=1)
            total_df_with_std.index = whole_mol_df_aly.index
            self.plot_fe_he_std_dir[aly] = total_df_with_std
        self.plot_fe_he_std_dir['fe'] = whole_mol_fe
        self.plot_fe_he_std_dir['fe_std'] = whole_mol_fe_std
        self.gen_whole_fe_time_serials_csv(self.plot_fe_he_std_dir)
        plot_obj = PLOTTING()
        plot_obj.plot_fe_time_serial(png_name, **self.plot_fe_he_std_dir)

    def gen_whole_fe_time_serials_csv(self, plot_fe_he_std_dir):
        for key, value in plot_fe_he_std_dir.items():
            time_serials_type = key
            if time_serials_type in ['forward', 'reverse', 'moving']:
                value.to_csv(f'{time_serials_type}_whole_mol_fe_time_serial.csv', sep=',', index=False)

class optParser():  
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-s', '--system_name', dest='system_name', help="The system name of the fe calculation.", default='1ke7')
        parser.add_option('-p', '--root_path', dest='root_path', help="The root path.",)
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

if __name__ == "__main__":
    opts = optParser('')
    system_name = str(opts.option.system_name)
    root_path = str(opts.option.root_path)
    df_name_pattern = '_bar_'
    lig_lst = [system_name]
    plotting_plan=['moving', 'forward', 'reverse']
    for lig_name in lig_lst:
        for side in ['complex', 'ligand']:
            base_path = os.path.join(root_path, lig_name, 'openmm_run', side)
            png_name = f'{lig_name}_{side}_fe_time_serial.png'
            # try:
            ts = TimeSeriesDataFrameCombiner(df_name_pattern, base_path, plotting_plan)
            plot_fe_he_std_dir = ts.plot_fe_he_std_dir
            ts.plotting_combined_time_serial(plot_fe_he_std_dir, png_name)
            if side == 'complex':
                com_plot_fe_he_std_dir = ts.plot_fe_he_std_dir
            else:
                lig_plot_fe_he_std_dir = ts.plot_fe_he_std_dir
    # res_aly = get_aly_restr_lig(os.path.join(root_path, lig_name, 'openmm_run', 'complex', 'res_databystd.csv'))
    with open('aly_res_dg.txt', 'r') as f:
        res_aly = float(f.read())
    # print(com_plot_fe_he_std_dir)
    # print(lig_plot_fe_he_std_dir)
    whole_mol_png_name = f'{lig_name}_whole_mol_fe_time_serial.png'
    abfe_combiner = ABFE_whole_mol_time_series_combiner(res_aly, com_plot_fe_he_std_dir, lig_plot_fe_he_std_dir, whole_mol_png_name)

            # plot whole mol time serial

            # except:
            #     print(f'{lig_name}_{side} plotting error!')

    
    
    


