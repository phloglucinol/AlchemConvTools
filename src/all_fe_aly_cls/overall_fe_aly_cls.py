import pandas as pd
from ..plotting.plotting_tools import PLOTTING
import multiprocessing
import os 
from ..common_tools.get_aly_restr_lig import get_aly_restr_lig
from ..parsing import read_openmm_out, read_amber_out, read_gmx_out
from ..convergence.FE_cal_convergence import CAL_FREE_ENERGY
from ..convergence.Time_series_convergence import ANA_FEP_TOOLS
import re
from PIL import Image
from contextlib import contextmanager
from glob import glob

@contextmanager
def bf_af_plot(mkdir_, lambda_info):
    pwd = os.getcwd()
    if os.path.exists(mkdir_) and os.path.isdir(mkdir_):
        pass
    else:
        os.makedirs(mkdir_)
    os.chdir(mkdir_)
    print(f"Be aware that the sequence numbers in the output images(stored in {mkdir_}) indicate the corresponding window index information, which remains consistent with the input's wanted_win_lst.")
    with open('lambda_info', 'w') as f:
        print(lambda_info, file=f)
    yield
    os.chdir(pwd) 


class ANA_OVERALL():
    def __init__(self,):
        pass
    
    @staticmethod
    def get_data_after_time_serial_anal(ana_fep_tools, lambda_, get_which_two_columns, to_do_strategy=['forward', 'reverse', 'moving']):
        fe_columns_name = get_which_two_columns[0]
        std_columns_name = get_which_two_columns[1]
        fe_he_std_dir = {}
        for need_do_stra in to_do_strategy:
            df = getattr(ana_fep_tools, f'{need_do_stra}_esti_obj').df_dict_with_std[lambda_]
            two_col_df = pd.concat([df[fe_columns_name], df[std_columns_name]],axis=1)
            fe_he_std_dir[need_do_stra] = two_col_df
        return fe_he_std_dir

    def single_side_time_serial_anal(self, side_name, iflog_csv, to_do_strategy=['forward', 'reverse', 'moving'], ifFEP=False, fep_forward=True, std_step_ratio=0.1, divided_ratio=0.01, block_width_ratio_formoving=0.2, lambda_seq_index = 'last'):
        if ifFEP:
            if fep_forward:
                get_which_two_columns=['FEP_forward_bythislambda(kcal/mol)', 'FEP_forward_Rolling_STD(kcal/mol)']
            else:
                get_which_two_columns=['FEP_reverse_bythislambda(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
        else:
            get_which_two_columns=['free_energy(kcal/mol)', 'Rolling_STD(kcal/mol)']
        ana_fep_tools = getattr(self, f'{side_name}_ANA_FEP_TOOLS')
        for spec_stra in to_do_strategy:
            lambda_dict = ana_fep_tools.generate_data_dict(False, 1) #fraction=1, no need to assign the fraction because the self.f'{side_name}_ANA_FEP_TOOLS' has been assigned the fraction 
            # print(lambda_dict)
            stra_estimate_method = getattr(ana_fep_tools, f'{spec_stra}_estimate')
            if spec_stra == 'moving':
                stra_estimate_method(lambda_dict, std_step_ratio, divided_ratio, block_width_ratio_formoving, ifFEP)
            else:
                # print(lambda_dict, 'forwarddddddddddddddddddddddddddd')
                stra_estimate_method(lambda_dict, std_step_ratio, divided_ratio, ifFEP)
            if iflog_csv:
                if ifFEP:
                    getattr(ana_fep_tools, f'{spec_stra}_esti_obj').output_csv('_fep_')
                else:
                    getattr(ana_fep_tools, f'{spec_stra}_esti_obj').output_csv('_bar_')
            if lambda_seq_index == 'last':
                lambda_seq_ele = getattr(ana_fep_tools, f'{spec_stra}_esti_obj').lambda_seq[-1] 
            else:
                lambda_seq_ele = getattr(ana_fep_tools, f'{spec_stra}_esti_obj').lambda_seq[lambda_seq_index]
        fe_and_std_dir = type(self).get_data_after_time_serial_anal(ana_fep_tools, lambda_seq_ele, get_which_two_columns, to_do_strategy)

        return fe_and_std_dir
    
    @staticmethod
    def get_festd_df_cal_all_fe(side_df_dict, code_fe_df, code_std_df):
        side_df_dict = side_df_dict
        fe_df = eval(code_fe_df)
        std_df = eval(code_std_df)
        all_estimated_FE_df = pd.concat([fe_df, std_df], axis = 1)
        return all_estimated_FE_df
    
    
    def plot_time_serial_convergence(self, plot_plan, fe=None, fe_std=None, png_name='bar_sum_all_plot_fbm.png'):
        plot_obj = PLOTTING()
        if len(plot_plan) == 3:
            plot_fe_he_std_dir = {
                'forward': self.forward_df_allfestd,
                'reverse': self.reverse_df_allfestd,
                'moving':  self.moving_df_allfestd,
                'fe': fe,
                'fe_std': fe_std,
            }
        elif len(plot_plan) == 2: # ['forward', 'reverse'] or ['forward', 'moving'] or ['reverse', 'moving'] will be leading to ['forward', 'reverse']
            plot_fe_he_std_dir = {
                'forward': self.forward_df_allfestd,
                'reverse': self.reverse_df_allfestd,
                'fe': fe,
                'fe_std': fe_std,
            }
        elif len(plot_plan) == 1:
            plot_fe_he_std_dir = {
                f'{plot_plan[0]}':  getattr(self, f'{plot_plan[0]}_df_allfestd'),
                'fe': fe,
                'fe_std': fe_std,
            }
        plot_obj.plot_fe_time_serial(png_name, **plot_fe_he_std_dir)

# class ANA_ABFE_OVERALL(ANA_OVERALL): 维护中，暂时弃用
#     def __init__(self, sys_path,  index_col, delimiter):
#         self.sys_path = sys_path
#         self.com_path_pattern = os.path.join(sys_path, 'complex', 'state_s*.csv')
#         self.lig_path_pattern = os.path.join(sys_path, 'ligand', 'state_s*.csv')
#         self.res_info_path = os.path.join(sys_path, 'complex', 'res_databystd.csv')
#         self.aly_restr_lig_ene = get_aly_restr_lig(self.res_info_path)# in kcal/mol
#         processes_num = multiprocessing.cpu_count()

#         self.com_u_nk_pd = READ_PROD_OUT(self.com_path_pattern).extract_data(processes_num, index_col, delimiter)
#         self.lig_u_nk_pd = READ_PROD_OUT(self.lig_path_pattern).extract_data(processes_num, index_col, delimiter)
        
#         self.com_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.com_u_nk_pd)
#         self.lig_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.lig_u_nk_pd)
        
#     def normal_abfe_time_serial_analysis(self, iflog_csv, output_strategy=['forward', 'reverse', 'moving'], ifPlot=False, fraction=1):
#         for side in ['com', 'lig']:
#             setattr(self, f'{side}_fe_and_std_dir', self.single_side_time_serial_anal(side, iflog_csv, output_strategy, fraction))
#             for i in output_strategy:
#                 fe_and_std_dir = getattr(self, f'{side}_fe_and_std_dir')
#                 setattr(self, f'{i}_df_{side}', fe_and_std_dir[i])
                
#                 if iflog_csv:
#                     getattr(self, f'{i}_df_{side}').to_csv(f'{i}_df_{side}.csv', sep='|')
#         for i in output_strategy:
#             spec_stra_com_festd_df = getattr(self, f'{i}_df_com')
#             spec_stra_lig_festd_df = getattr(self, f'{i}_df_lig')
#             spec_stra_side_df_dict = {'com': spec_stra_com_festd_df, 
#                                       'lig': spec_stra_lig_festd_df,
#                                       'aly_restr':self.aly_restr_lig_ene,
#                                      }
#             code_fe_df = "side_df_dict['lig'].iloc[:,0]+side_df_dict['aly_restr']-side_df_dict['com'].iloc[:,0]"
#             code_std_df = "(side_df_dict['lig'].iloc[:,1]**2+side_df_dict['com'].iloc[:,1]**2)**0.5"
#             spec_stra_all_festd_df = type(self).get_festd_df_cal_all_fe(spec_stra_side_df_dict, code_fe_df, code_std_df)
#             setattr(self, f'{i}_df_allfestd', spec_stra_all_festd_df)
#             if iflog_csv:
#                 getattr(self, f'{i}_df_allfestd').to_csv(f'{i}_df_allfestd.csv', sep='|')
#         if ifPlot:
#             # assert output_strategy == ['forward', 'reverse', 'moving'], "Plot option only support when output_strategy equal to ['forward', 'reverse', 'moving']"
#             if output_strategy == ['forward', 'reverse']:
#                 plot_plan = 2
#             elif output_strategy == ['forward', 'reverse', 'moving']:
#                 plot_plan = 3
#             elif output_strategy == ['moving']:
#                 plot_plan = 1
#             self.plot_time_serial_convergence(plot_plan)



class singleside_ANA_ABFE_OVERALL(ANA_OVERALL):
    def __init__(self, sys_path, simulation_pack, file_prefix, file_suffix, index_col=[0,1,2,3], delimiter='|', ifsubsample=False):
        self.sys_path = sys_path
        self.simulation_pack=simulation_pack
        self.file_prefix = file_prefix
        self.file_suffix = file_suffix
        if self.simulation_pack == 'openmm':
            self.openmm_init(self.sys_path, file_prefix, file_suffix, index_col, delimiter)
        elif self.simulation_pack == 'amber':
            self.amber_init(self.sys_path, file_prefix, file_suffix,)
        elif self.simulation_pack == 'gmx':
            self.gmx_init(self.sys_path, file_prefix, file_suffix,)
        self.com_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.com_u_nk_pd)
        if ifsubsample:
            self.subsample_correlated_data()

        self.final_fe = 0.0 # will be updated by self.cal_fe
        self.final_fe_std = 0.0 # will be updated by self.cal_fe
    
    def plot_deltaU(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.plot_simulation_deltaU(50)
        

    def subsample_correlated_data(self,):
        lambda_info_list = list(self.com_u_nk_pd.index.names)
        lambda_info_list.pop(0)
        self.simulation_win_lst = [i for i,j in self.com_u_nk_pd.groupby(lambda_info_list, sort=False)]
        subsample_lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(self.simulation_win_lst, 1, True)
        self.subsample_lambda_unk_df = pd.concat(list(subsample_lambda_dict.values()))
        self.com_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.subsample_lambda_unk_df)

    def openmm_init(self, sys_path, file_prefix, file_suffix, index_col, delimiter, ):
        processes_num = multiprocessing.cpu_count()
        read_openmm_obj = read_openmm_out.READ_PROD_OUT(sys_path, file_prefix, file_suffix)
        self.com_u_nk_pd = read_openmm_obj.extract_data(processes_num, read_openmm_obj.index_col, delimiter)

    def amber_init(self, sys_path, file_prefix, file_suffix, ):
        self.com_u_nk_pd = read_amber_out.ReadProdOut(sys_path, file_prefix, file_suffix).extract_data()
    
    def gmx_init(self, sys_path, file_prefix, file_suffix,):
        com_u_nk_pd = read_gmx_out.ReadProdOut(sys_path, file_prefix, file_suffix).extract_data()
        com_u_nk_pd.columns = [(i[2],i[0],i[1]) for i in com_u_nk_pd.columns]
        tuples = [(i[0], i[3], i[1], i[2]) for i in com_u_nk_pd.index]
        simulation_lambda_names = ['time', 'lambda_restraints', 'lambda_electrostatics', 'lambda_sterics']
        com_u_nk_pd.index = pd.MultiIndex.from_tuples(tuples, name=simulation_lambda_names)
        self.com_u_nk_pd = com_u_nk_pd

    # def reweight_analysis(self, use_wins_num, unit='kcal/mol', ifforward=True, fraction=1):
    #     reweight_obj = REWEIGHT(self.com_u_nk_pd, True, fraction)
    #     reweight_obj.cal_FE_once_reweight(use_wins_num, unit, forward=ifforward)

    def check_MeasurementOfConvergence(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.check_MeasurementOfConvergence_overall("MeasurementOfConvergence.csv")

    def pai_checking(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.check_pai_overall("pai_check.csv")

    def reweight_entropy_checking(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.check_reweight_entropy_overall("reweight_entropy_check.csv")

    def cfm_checking(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.check_cfm_overall(True)

    def dc_overlap_checking(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.check_DC_MBAR_overlap_overall()
            
    def mbar_overlap_matrix_checking(self, wanted_win_lst, output_file_dir, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.plot_overlap_matrix('mbar_matrix.png')   
 
    def Boltzmann_weight_PdU_checking(self, wanted_win_lst, output_file_dir, fraction=1):       
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_cal_FEP_TOOLS = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        with bf_af_plot(output_file_dir, wanted_win_lst):
            com_cal_FEP_TOOLS.plot_Boltzmann_weight_PdU()
 
    def cal_fe(self, output_result_dir, output_filename, wanted_win_lst, unit, fe_mode, std_mode, ifplot_du, fraction=1):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        # print(lambda_dict)
        with bf_af_plot(output_result_dir, wanted_win_lst):
            fe_df = self.com_ANA_FEP_TOOLS.output_fe_and_std(lambda_dict, unit, fe_mode, std_mode, ifdetail=True, ifplot=ifplot_du)
            fe_df.to_csv(output_filename, sep="|")
            self.final_fe_std = fe_df.iloc[-1,-1]
            self.final_fe = fe_df.iloc[-1,-2]
            # print(self.final_fe, self.final_fe_std)

    def cal_fe_reweighting(self, output_result_dir, wanted_win_lst, use_wins_step, use_wins_num, unit, forward, ifdiagnum, postfix, fraction=1, ifplot=False, reweight_err_bar_max=2.0):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        with bf_af_plot(output_result_dir, wanted_win_lst):
            fe_df = self.com_ANA_FEP_TOOLS.output_fe_reweighting(lambda_dict, use_wins_step, use_wins_num, unit, forward, ifdiagnum, postfix, ifplot, reweight_err_bar_max)
    
    def time_serial_reweight(self, output_result_dir, wanted_win_lst, fraction, divided_ratio=0.01, block_width_ratio_formoving=0.2, plot_plan=['forward', 'reverse',], png_prefix='FEP', wins_step=0.1, wins_num=3, forward=True, ifdiagnum=True, reweight_err_bar_max=2.0):
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        with bf_af_plot(output_result_dir, wanted_win_lst):
            self.com_ANA_FEP_TOOLS.reweight_check_time_serial(lambda_dict, plot_plan, png_prefix, divided_ratio, block_width_ratio_formoving, wins_step, wins_num, forward, ifdiagnum, reweight_err_bar_max)

    def time_serial_check(self, output_result_dir, wanted_win_lst, fraction, std_step_ratio=0.1, divided_ratio=0.01, block_width_ratio_formoving=0.2, plot_strategy=['forward', 'reverse',], ifplot_bar_time_series=True, ifuse_fep_check_everywins=True, iflog_fep_csv=True, iflog_bar_csv=True, ifplot_fep_time_series=True, use_forward=True):
        print('Warning: The fe and fe_std should be calculated by self.cal_fe, before you check time_serial_convergence!')
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        # print('---------------------------------------------',lambda_unk_df,'----------')
        # print(lambda_unk_df.columns, type(lambda_unk_df.columns[0]))
        ori_self_com_ANA_FEP_TOOLS = self.com_ANA_FEP_TOOLS
        self.com_ANA_FEP_TOOLS = ANA_FEP_TOOLS(lambda_unk_df)# will be used by the self.single_side_time_serial_anal
        plot_strategy_str = '_'.join(plot_strategy)
        with bf_af_plot(output_result_dir, wanted_win_lst):
            if ifuse_fep_check_everywins:
                for i in range(len(wanted_win_lst)):
                    png_name = f'fep_time_series_{plot_strategy_str}_{i}.png'
                    self.single_side_time_serial_analysis(iflog_csv=iflog_fep_csv, side_name="com", output_strategy=plot_strategy, ifPlot=ifplot_fep_time_series, png_name=png_name,ifFEP=True, iffep_forward=use_forward, std_step_ratio=std_step_ratio, divided_ratio=divided_ratio, block_width_ratio_formoving=block_width_ratio_formoving, lambda_seq_index = i)
            png_name = f'bar_time_series_{plot_strategy_str}_sum_all.png'
            self.single_side_time_serial_analysis(iflog_csv=iflog_bar_csv, side_name="com", output_strategy=plot_strategy, ifPlot=ifplot_bar_time_series, png_name=png_name,ifFEP=False, iffep_forward=use_forward, std_step_ratio=std_step_ratio, divided_ratio=divided_ratio, block_width_ratio_formoving=block_width_ratio_formoving, lambda_seq_index = 'last')
        self.com_ANA_FEP_TOOLS = ori_self_com_ANA_FEP_TOOLS

    def single_side_time_serial_analysis(self, iflog_csv, side_name, output_strategy=['forward', 'reverse', 'moving'], ifPlot=False, png_name='time_series.png',ifFEP=False, iffep_forward=True, std_step_ratio=0.1, divided_ratio=0.01, block_width_ratio_formoving=0.2, lambda_seq_index = 'last'):
        setattr(self, f'{side_name}_fe_and_std_dir', self.single_side_time_serial_anal(side_name, iflog_csv, output_strategy, ifFEP, iffep_forward, std_step_ratio, divided_ratio, block_width_ratio_formoving, lambda_seq_index))
        for i in output_strategy:
            fe_and_std_dir = getattr(self, f'{side_name}_fe_and_std_dir')
            setattr(self, f'{i}_df_{side_name}', fe_and_std_dir[i])
            
            if iflog_csv:
                getattr(self, f'{i}_df_{side_name}').to_csv(f'{i}_df_{side_name}.csv', sep='|')
        for i in output_strategy:
            spec_stra_com_festd_df = getattr(self, f'{i}_df_{side_name}')
            spec_stra_side_df_dict = {f'{side_name}': spec_stra_com_festd_df, }
            code_fe_df = f"side_df_dict['{side_name}'].iloc[:,0]"
            # print("--------------",f'spec:{spec_stra_side_df_dict}','-----------------\n',f'code_fe:{code_fe_df}')
            code_std_df = f"side_df_dict['{side_name}'].iloc[:,1]"
            spec_stra_all_festd_df = type(self).get_festd_df_cal_all_fe(spec_stra_side_df_dict, code_fe_df, code_std_df)
            setattr(self, f'{i}_df_allfestd', spec_stra_all_festd_df)
        if ifPlot:
            if lambda_seq_index == 'last':
                self.plot_time_serial_convergence(output_strategy, self.final_fe, self.final_fe_std, png_name)
            else:
                self.plot_time_serial_convergence(output_strategy, None, None, png_name)

    def dG_dlambda_time_serial_analysis(self, output_result_dir, wanted_win_lst, fraction, std_step_ratio=0.1, divided_ratio=0.1, ifplot_bar_dG_dlambda=True, ifplot_fep_dG_dlambda=True, use_forward=True):
        def save_gif(fe_mode):
            images = []
            cur_path = os.getcwd()
            pngfile_list = sorted(glob(os.path.join(cur_path, f'*{fe_mode}_dG_dlambda*.png')),key=os.path.getmtime)
            print(pngfile_list)
            for pngfile in pngfile_list:
                image = Image.open(pngfile)
                images.append(image)
            output_file = f'{fe_mode}_dG_dlambda.gif'
            images[0].save(output_file, save_all=True, append_images=images[1:], duration=500, loop=1)
        lambda_dict = self.com_ANA_FEP_TOOLS.generate_data_dict(wanted_win_lst, fraction)
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        com_ANA_FEP_TOOLS = ANA_FEP_TOOLS(lambda_unk_df)
        with bf_af_plot(output_result_dir, wanted_win_lst):
            if ifplot_fep_dG_dlambda:
                com_ANA_FEP_TOOLS.forward_estimate(lambda_dict, std_step_ratio, divided_ratio, iffep=True, ifplottime_serials=False, ifplotdG_dlambda=True, ifuseforward=use_forward)
                save_gif('FEP')
            if ifplot_bar_dG_dlambda:
                com_ANA_FEP_TOOLS.forward_estimate(lambda_dict, std_step_ratio, divided_ratio, iffep=False, ifplottime_serials=False, ifplotdG_dlambda=True, ifuseforward=use_forward)
                save_gif('BAR')
        