# from email.policy import default
from optparse import OptionParser
from src.parsing.input_file_parser import InputParser as InputParser
from src.all_fe_aly_cls.overall_fe_aly_cls import singleside_ANA_ABFE_OVERALL as singleside_ANA_ABFE_OVERALL
import os




class optParser():  
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-i', '--input', dest='input', help="The file name of input, which recording the analysis settings. Default: 'input.txt'", default='input.txt')
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

opts = optParser('')
#####################################################################################################
# keep for jupyter notebook test
# fakeArgs = '-i input.txt'
# opts = optParser(fakeArgs.strip().split())
#####################################################################################################
input_parser = InputParser(opts.option.input)
basic_settings = input_parser.get_Basic_settings()
# print(basic_settings)
convergence_analysis_settings = input_parser.get_Convergence_analysis()
# print(convergence_analysis_settings)


# Basic settings
simulation_pack = str(basic_settings['simulation_software'])
file_str = str(basic_settings['file_directory']) # For openmm, it is the directory that stores the state_s*.csv or state_g*_s*.csv files.
file_prefix = str(basic_settings['file_prefix'])
file_suffix = str(basic_settings['file_suffix'])
fraction = float(basic_settings['fraction'])
output_csv_filename = str(basic_settings['output_csv_filename'])
energy_unit = str(basic_settings['energy_unit'])
std_mode = str(basic_settings['std_mode'])
cal_wins = str(basic_settings['calculation_windows'])
ifplot_du = basic_settings['plot_du']
ifsubsample = basic_settings['subsample']

# Convergence analysis settings
ifTime_serials_aly = convergence_analysis_settings['Time_series']
ifBoltzmann_weight_PdU = convergence_analysis_settings['Boltzmann_weight_PdU']
ifCurve_Fitting_Method = convergence_analysis_settings['Curve_Fitting_Method']
ifBennett_Overlap_Value = convergence_analysis_settings['Bennett_Overlap_Value']
ifPai_Value = convergence_analysis_settings['Pai_Value']
ifMeasurementOfConvergence_Value = convergence_analysis_settings['MeasurementOfConvergence_Value']
ifReweighting_SelfConsistent = convergence_analysis_settings['Reweighting_SelfConsistent']

# Reweighting analysis settings
if ifReweighting_SelfConsistent:
    reweighting_analysis_settings = input_parser.get_Reweighting_analysis()
    reweight_use_wins_step = float(reweighting_analysis_settings['reweight_use_wins_step'])
    reweight_use_wins_num = int(reweighting_analysis_settings['reweight_use_wins_num'])
    ifforward = reweighting_analysis_settings['ifforward']
    ifdiagnum = reweighting_analysis_settings['ifdiagnum']
    ifreweight_time_serial = reweighting_analysis_settings['ifreweight_time_serial']
    ifplotheatmap = reweighting_analysis_settings['ifplotheatmap']
    reweight_err_bar_max = reweighting_analysis_settings['reweight_err_bar_max']

else:
    reweight_use_wins_step = 0.0
    reweight_use_wins_num = 0
    ifforward = False
    ifdiagnum = False
    ifreweight_time_serial = False
    ifplotheatmap = False
    reweight_err_bar_max = 2


time_serials_analysis_settings = input_parser.get_Time_series_analysis()
ifuseFEP_check_everywins_Time_series = time_serials_analysis_settings['useFEP_check_everywins_Time_series']
plot_plan = time_serials_analysis_settings['plot_plan'].strip().split('|')
std_step_ratio = time_serials_analysis_settings['std_step_ratio']
divided_ratio = time_serials_analysis_settings['divided_ratio']
dGdlambda_divided_ratio = time_serials_analysis_settings['dGdlambda_divided_ratio']
block_width_ratio_formoving = time_serials_analysis_settings['width_ratio']
ifplot_bar_time_series = time_serials_analysis_settings['ifplot_bar_time_series']
ifplot_bar_dG_dlambda = time_serials_analysis_settings['ifplot_bar_dG_dlambda']
ifplot_fep_dG_dlambda = time_serials_analysis_settings['ifplot_fep_dG_dlambda']
print(ifplot_bar_time_series)
print(type(ifplot_bar_time_series))
iflog_fep_csv = time_serials_analysis_settings['iflog_fep_csv']
iflog_bar_csv = time_serials_analysis_settings['iflog_bar_csv']
ifplot_fep_time_series = time_serials_analysis_settings['ifplot_fep_time_series']
use_forward = time_serials_analysis_settings['use_forward']




fe_cls = singleside_ANA_ABFE_OVERALL(file_str, simulation_pack, file_prefix, file_suffix, [0,1,2,3], '|', ifsubsample)

if cal_wins == 'all':
    lambda_info_list = list(fe_cls.com_u_nk_pd.index.names)
    lambda_info_list.pop(0)
    simulation_win_lst = [i for i,j in fe_cls.com_u_nk_pd.groupby(lambda_info_list, sort=False)]
    force_field_cal_win_lst = list(fe_cls.com_u_nk_pd.columns)
    cal_win_lst = simulation_win_lst
    cal_win_lst = [cal_win_lst, ]
else:
    with open(cal_wins, 'r') as f:
        f_content = f.readlines()
        cal_win_lst = [eval(line.strip()) for line in f_content]
        import sys 
        if not isinstance(cal_win_lst, list):
            print('Error: The content of the win_lst file is not a list!')
            sys.exit()
# print(cal_win_lst[0])

def jobs(fe_cal_cls, store_path, cal_win, energy_unit, std_mode, fraction, ifplot_du, iftime_serials_aly, ifuseFEP_check_everywins_Time_series, plot_plan, std_step_ratio, divided_ratio, block_width_ratio_formoving, ifplot_bar_time_series, iflog_fep_csv, iflog_bar_csv, ifplot_fep_time_series, use_forward, ifBoltzmann_weight_PdU, ifCurve_Fitting_Method, ifBennett_Overlap_Value, ifPai_Value, ifMeasurementOfConvergence_Value, ifReweighting_SelfConsistent, reweight_use_wins_step, reweight_use_wins_num, ifforward, ifdiagnum, ifplotheatmap, ifreweight_time_serial, reweight_err_bar_max, dGdlambda_divided_ratio, ifplot_bar_dG_dlambda, ifplot_fep_dG_dlambda):
    fe_cal_cls.cal_fe(store_path, output_csv_filename, cal_win, energy_unit, 'BAR', std_mode, False, fraction)
    if ifplot_du:
        fe_cal_cls.plot_deltaU(cal_win, os.path.join(store_path,'du_distribution'), fraction)
    if ifBoltzmann_weight_PdU:
        fe_cal_cls.Boltzmann_weight_PdU_checking(cal_win, os.path.join(store_path,'Boltzmann_weight_PdU'), fraction)
    if ifCurve_Fitting_Method:
        fe_cal_cls.cfm_checking(cal_win, os.path.join(store_path,'cfm_checking'), fraction)
    if ifBennett_Overlap_Value:
        fe_cal_cls.dc_overlap_checking(cal_win, os.path.join(store_path,'dc_overlap_checking'), fraction)
    if ifPai_Value:
        fe_cal_cls.pai_checking(cal_win, os.path.join(store_path,'pai_checking'), fraction)
    if ifMeasurementOfConvergence_Value:
        fe_cal_cls.check_MeasurementOfConvergence(cal_win, os.path.join(store_path,'MeasurementOfConvergence_checking'), fraction)
    if iftime_serials_aly:
        fe_cal_cls.time_serial_check(os.path.join(store_path,'time_serial_check'), cal_win, fraction, std_step_ratio, divided_ratio, block_width_ratio_formoving, plot_plan, ifplot_bar_time_series, ifuseFEP_check_everywins_Time_series, iflog_fep_csv, iflog_bar_csv, ifplot_fep_time_series, use_forward)
    if ifplot_bar_dG_dlambda or ifplot_fep_dG_dlambda:
        fe_cal_cls.dG_dlambda_time_serial_analysis(os.path.join(store_path,'dG_dlambda_time_serials'), cal_win, fraction, std_step_ratio, dGdlambda_divided_ratio, ifplot_bar_dG_dlambda, ifplot_fep_dG_dlambda, use_forward)
    if ifReweighting_SelfConsistent:
        fe_cal_cls.cal_fe_reweighting(store_path, cal_win, reweight_use_wins_step, reweight_use_wins_num, energy_unit, ifforward, ifdiagnum, 'timeall', fraction, ifplotheatmap, reweight_err_bar_max)
        if ifreweight_time_serial:
            fe_cal_cls.time_serial_reweight(store_path, cal_win, fraction, divided_ratio, block_width_ratio_formoving, plot_plan, 'reweighting_selfconsistent', reweight_use_wins_step, reweight_use_wins_num, ifforward, ifdiagnum, reweight_err_bar_max)
    

if len(cal_win_lst) == 1:
    print(ifplot_fep_time_series)
    print(ifplot_fep_time_series)
    jobs(fe_cls, 'fe_cal_out', cal_win_lst[0], energy_unit, std_mode, fraction, ifplot_du, ifTime_serials_aly, ifuseFEP_check_everywins_Time_series, plot_plan, std_step_ratio, divided_ratio, block_width_ratio_formoving, ifplot_bar_time_series, iflog_fep_csv, iflog_bar_csv, ifplot_fep_time_series, use_forward, ifBoltzmann_weight_PdU, ifCurve_Fitting_Method, ifBennett_Overlap_Value, ifPai_Value, ifMeasurementOfConvergence_Value, ifReweighting_SelfConsistent, reweight_use_wins_step, reweight_use_wins_num, ifforward, ifdiagnum, ifplotheatmap, ifreweight_time_serial, reweight_err_bar_max, dGdlambda_divided_ratio, ifplot_bar_dG_dlambda, ifplot_fep_dG_dlambda)
else:
    for idx_ in range(0, len(cal_win_lst)):
        win_lst = cal_win_lst[idx_]
        print(win_lst)
        jobs(fe_cls, f'fe_cal_out_win_lst_{idx_}', win_lst, energy_unit, std_mode, fraction, ifplot_du, ifTime_serials_aly, ifuseFEP_check_everywins_Time_series, plot_plan, std_step_ratio, divided_ratio, block_width_ratio_formoving, ifplot_bar_time_series, iflog_fep_csv, iflog_bar_csv, ifplot_fep_time_series, use_forward, ifBoltzmann_weight_PdU, ifCurve_Fitting_Method, ifBennett_Overlap_Value, ifPai_Value, ifMeasurementOfConvergence_Value, ifReweighting_SelfConsistent, reweight_use_wins_step, reweight_use_wins_num, ifforward, ifdiagnum, ifplotheatmap, ifreweight_time_serial, reweight_err_bar_max, dGdlambda_divided_ratio, ifplot_bar_dG_dlambda, ifplot_fep_dG_dlambda)
