import pandas as pd
from plotting_tools import PLOTTING
import numpy as np
import os
from glob import glob

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
    def __init__(self, df_basic_name, base_path, sub_string_dir):
        self.df_basic_name = df_basic_name # free_ene_SUBSTI.csv
        self.base_path = base_path # /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out
        self.sub_string_dir = sub_string_dir # {'group1': 'complex/group1', 'group2': 'complex/group2', 'group3': 'complex/group3', 'group4': 'complex/group4'}
        self.df_name_dir = self.gen_dir_name_dir(sub_string_dir)
        self.combine_obj = SumFEDataFrameCombiner(self.df_name_dir, "|")
        self.total_df_with_std = self.combine_obj.combine_way('bar_std')
        self.fe = self.total_df_with_std.iloc[-1, 0]
        self.fe_std = self.total_df_with_std.iloc[-1, 1]       

    def gen_dir_name_dir(self, sub_string_dir, ):
        df_name_dir = {}
        for key, value in sub_string_dir.items():
            df_name = self.df_basic_name.replace('SUBSTI', key)
            prefix_path = self.base_path.replace('SUBSTI', value)
            df_name_dir[key] = f'{prefix_path}/{df_name}'
        return df_name_dir

class TimeSeriesDataFrameCombiner():  
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_reverse_esti.csv
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_moving_esti.csv
    # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_forward_esti.csv
    # base_path: /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out/time_serial_check
    
    def __init__(self, df_basic_name, base_path, sub_string_dir, plotting_plan=['moving','forward', 'reverse']):
        self.df_basic_name = df_basic_name # _bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_SUBSTI_esti.csv
        self.base_path = base_path # /HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out
        self.sub_string_dir = sub_string_dir # {'group1': 'complex/group1', 'group2': 'complex/group2', 'group3': 'complex/group3', 'group4': 'complex/group4'}
        self.plot_fe_he_std_dir = {}
        self.FE_combine_obj = FreeEnergyDataFrameCombiner('free_ene_SUBSTI.csv', self.base_path, self.sub_string_dir)
        
        for aly in plotting_plan:
            df_name_dir = self.gen_dir_name_dir(self.sub_string_dir, aly)
            combine_obj = SumTimeFEDataFrameCombiner(df_name_dir) # different combiner
            combine_obj.combine_way()
            self.plot_fe_he_std_dir[aly] = combine_obj.cal_rolling_std(0.1)
        self.plot_fe_he_std_dir['fe'] = self.FE_combine_obj.fe
        self.plot_fe_he_std_dir['fe_std'] = self.FE_combine_obj.fe_std
    
    def gen_dir_name_dir(self, sub_string_dir, type_):
        # setattr(self, f'{type_}_df_name', self.df_basic_name.replace('SUBSTI', type_))
        df_name = self.df_basic_name.replace('SUBSTI', type_)
        df_name_dir = {}
        for key, value in sub_string_dir.items():
            prefix_path = self.base_path.replace('SUBSTI', value)
            df_name_dir[key] = f'{prefix_path}/time_serial_check/{df_name}'
        return df_name_dir
    
    def plotting_combined_time_serial(self, plot_fe_he_std_dir, png_name):
        plot_obj = PLOTTING()
        plot_obj.plot_fe_time_serial(png_name, **plot_fe_he_std_dir)

def generate_sub_string_dir(side_name, group_num, if_start_zero):
    sub_string_dir = {}
    for i in range(1, group_num + 1):
        group_name = f"group{i}"
        sub_string = f"{side_name}/{group_name}"
        if if_start_zero:
            group_name = f"group{i-1}"
        sub_string_dir[group_name] = sub_string
    return sub_string_dir

def generate_input_parameter(root_path, lig_name, side, if_start_zero):
    # root_path: /HOME/scz1641/run/BACE1_bygroup/3rd_batch/
    group_path = os.path.join(root_path,lig_name, 'openmm_run', side)
    pattern = 'group[0-9]*'
    group_directories = glob(os.path.join(group_path, pattern))
    group_num = len(group_directories)
    base_path = f"{root_path}/{lig_name}/openmm_run/SUBSTI/fe_cal_out"
    sub_string_dir = generate_sub_string_dir(side, group_num, if_start_zero)
    png_name = f"{lig_name}_{side}_{group_num}groups_vdw_time_serial.png"
    return base_path, sub_string_dir, png_name


if __name__ == "__main__":
    com_df_basic_name = "_bar__(1.0, 1.0, 0.0) to (1.0, 1.0, 1.0)_SUBSTI_esti.csv"
    lig_df_basic_name = "_bar__(0.0, 1.0, 0.0) to (0.0, 1.0, 1.0)_SUBSTI_esti.csv"
    root_path = "/HOME/scz1641/run/BACE1_bygroup/3rd_batch"
    lig_lst = ['lig_02', 'lig_03','lig_04','lig_05','lig_06', 'lig_07', 'lig_1a', 'lig_11', 'lig_36', 'lig_41', 'lig_45', 'lig_67', 'lig_13', 'lig_16', 'lig_32', 'lig_74', 'lig_81', 'lig_69', 'lig_03_5group_plip_seq', 'lig_03_5group_res_dist_based_seq', 'lig_03_9group_plip_seq', 'lig_03_9group_res_dist_based_seq', 'lig_05_8groups_manu', 'lig_07_8groups_manu', 'lig_13_11groups_manu', 'lig_1a_9groups_manu', 'lig_32_10groups_manu', 'lig_36_13groups_manu', 'lig_69_12groups_manu', 'lig_81_12groups_manu']
    lig_lst = ['lig_05_8groups_manu', 'lig_07_8groups_manu', 'lig_13_11groups_manu', 'lig_1a_9groups_manu', 'lig_32_10groups_manu', 'lig_36_13groups_manu', 'lig_69_12groups_manu', 'lig_81_12groups_manu']
    # base_path = "/HOME/scz1641/run/BACE1_bygroup/3rd_batch/lig_81/openmm_run/SUBSTI/fe_cal_out"
    # sub_string_dir = generate_sub_string_dir('complex', 8, False)
    plotting_plan=['moving','forward', 'reverse']
    # png_name = 'lig_81_com_8groups_vdw_time_serial.png'
    for lig_name in lig_lst:
        for side in ['complex', 'ligand']:
            base_path, sub_string_dir, png_name = generate_input_parameter(root_path, lig_name, side, False)
            if side == 'complex':
                try:
                    ts = TimeSeriesDataFrameCombiner(com_df_basic_name, base_path, sub_string_dir, plotting_plan)
                    plot_fe_he_std_dir = ts.plot_fe_he_std_dir
                    ts.plotting_combined_time_serial(plot_fe_he_std_dir, png_name)
                except:
                    print(f'{lig_name}_{side} plotting error!')
            else:
                try:
                    ts = TimeSeriesDataFrameCombiner(lig_df_basic_name, base_path, sub_string_dir, plotting_plan)
                    plot_fe_he_std_dir = ts.plot_fe_he_std_dir
                    ts.plotting_combined_time_serial(plot_fe_he_std_dir, png_name)
                except:
                    print(f'{lig_name}_{side} plotting error!')
    
    
    


