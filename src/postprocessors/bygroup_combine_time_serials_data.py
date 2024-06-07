import pandas as pd
from plotting_tools import PLOTTING
import numpy as np
import os
from glob import glob

def get_newest_file(pattern):
    # Use glob to find all files matching the pattern
    files = glob(pattern)

    # Check if at least one file was found
    if files:
        # Find the file with the latest modification time
        latest_file = max(files, key=os.path.getmtime)

        # Get the absolute path of the latest file
        latest_file_path = os.path.abspath(latest_file)

        # Print or return the path
        print(latest_file_path)
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
            self.df_dir[key] = pd.read_csv(value, delimiter)
        self.total_df = None

    def combine_way(self):
        pass

         
class SumTimeFEDataFrameCombiner(DataFrameCombiner):
    def combine_way(self):
        total_df = None
        for key, df in self.df_dir.items():
            if total_df is None:
                total_df = pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
            else:
                total_df += pd.DataFrame(df.loc[:, 'free_energy(kcal/mol)'])
        total_df.index = df.loc[:, 'time_ratio']
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
        self.combine_obj = SumFEDataFrameCombiner(self.df_name_dir, "|")
        self.total_df_with_std = self.combine_obj.combine_way('bar_std')
        self.fe = self.total_df_with_std.iloc[-1, 0]
        self.fe_std = self.total_df_with_std.iloc[-1, 1]       

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
        self.sub_string_dir = self.get_side_thermo_path(self.base_path, ['restraints*', 'electrostatics*', 'sterics*']) 
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
        return dirs_lst

if __name__ == "__main__":
    df_name_pattern = '_bar_'
    root_path = "/nfs/export3_25T/Bygroup_FEP_data/CDK2_bygroup/2nd_batch"
    lig_lst = ['1ke7']
    plotting_plan=['moving', 'forward', 'reverse']
    for lig_name in lig_lst:
        for side in ['complex', 'ligand']:
            base_path = os.path.join(root_path, lig_name, 'openmm_run', side)
            png_name = f'{lig_name}_{side}_fe_time_serial.png'
            try:
                ts = TimeSeriesDataFrameCombiner(df_name_pattern, base_path, plotting_plan)
                plot_fe_he_std_dir = ts.plot_fe_he_std_dir
                ts.plotting_combined_time_serial(plot_fe_he_std_dir, png_name)
            except:
                print(f'{lig_name}_{side} plotting error!')

    
    
    


