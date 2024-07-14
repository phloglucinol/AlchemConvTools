import math
import pandas as pd
import numpy as np
try:
    from pymbar import BAR as BAR_
except:
    from pymbar import bar as BAR_
from pymbar import timeseries as timeseries
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('agg')
from ..plotting.plotting_tools import PLOTTING
from PIL import Image
from .FE_cal_convergence import CAL_FREE_ENERGY
from ..plotting.dG_dlambda import dG_dlambda


class TIME_SERIAL_DATA():
    def __init__(self, time_serial_type, lambda_seq=[], columns=[]):
        self.time_serial_type = time_serial_type
        self.data_dict ={}
        self.lambda_seq = lambda_seq
        if lambda_seq==[]:
            pass
        else:
            print(f"The time_serial_type is {self.time_serial_type}")
            print('This is my lambda_seq:{}'.format(str(self.lambda_seq)))
        self.columns_list=columns
        self.columns_list.insert(0, 'time_ratio')
        for i in lambda_seq:
            self.data_dict[i]=[]
        self.df_dict={}
        self.df_dict_with_std = {}
        # print(self.data_dict)
        
    
    def update_data(self, time_ratio, data_):
        '''Update the free energy data, note that please send me the free energy data with unit of kcal*mol^(-1)
        
        '''
        lambda_seq = data_.index
        lambda_seq = lambda_seq.unique()
        # print(lambda_seq)
        # print(len(lambda_seq))
        for i in range(0, len(lambda_seq)):
            data_lst = list(data_.iloc[i,:])
            # print(data_lst)
            data_lst.insert(0, time_ratio)
            lambda_ = lambda_seq[i]
            # print(f'lambda_:{lambda_}; type:{type(lambda_)}')
            self.data_dict[lambda_].append(data_lst)
    
    def generate_df(self,):
        for i in self.data_dict.keys():
            onetime_data_df = pd.DataFrame(self.data_dict[i])
            onetime_data_df.columns=self.columns_list
            onetime_data_df.index = onetime_data_df.iloc[:,0]
#             print(onetime_data_df)
            self.df_dict[i] = onetime_data_df

    def update_data_reweight(self, time_ratio, data_):
        self.df_dict[time_ratio] = data_

    # def generate_df_reweight(self,):
    #     all_data_df = pd.DataFrame()
    #     for i in self.data_dict.keys():
    #         onetime_data_df = pd.DataFrame(self.data_dict[i])
    #         onetime_data_df.columns=self.columns_list
    #         onetime_data_df.index = [[i]*len(onetime_data_df.index)]
    #         all_data_df = pd.concat([all_data_df, onetime_data_df])
    
    def cal_rolling_std(self, std_step_ratio, iffep=True):        
        for key in self.df_dict.keys():
            onetime_data_df = self.df_dict[key]
            # print(onetime_data_df)
            std_step = int(np.floor(len(onetime_data_df)*std_step_ratio))
            if iffep:
                to_cal_df = onetime_data_df.loc[:, 'FEP_forward_bythislambda(kcal/mol)':'FEP_reverse_bythislambda(kcal/mol)']
                std_columns = ['FEP_forward_Rolling_STD(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            else:
                to_cal_df = pd.DataFrame(onetime_data_df.loc[:, 'free_energy(kcal/mol)'])
                std_columns = ['Rolling_STD(kcal/mol)',]
            # print(to_cal_df)
            std_df = to_cal_df.rolling(std_step).std()
#             print(std_df)
            for i in range(std_step):
                if iffep:
                    std_df.iloc[i,0]=std_df.iloc[std_step, 0]
                    std_df.iloc[i,1]=std_df.iloc[std_step, 1]
                else:
                    std_df.iloc[i,0]=std_df.iloc[std_step, 0]
                
            std_df.columns = std_columns
            onetime_data_df = pd.concat([onetime_data_df,std_df], axis = 1)
            self.df_dict_with_std[key] = onetime_data_df
            
    def output_csv(self, csv_prefix, ):
        for key in self.df_dict_with_std.keys():
            onetime_data_df = self.df_dict_with_std[key]
            csv_name = csv_prefix+'_'+str(key)+'_'+self.time_serial_type+'_esti'+'.csv'
            onetime_data_df.to_csv(csv_name)

    def output_reweight_csv(self, csv_prefix,):
        csv_name = csv_prefix+'_reweight'+self.time_serial_type+'_esti'+'.csv'
        self.all_data_df.to_csv(csv_name)


class ANA_FEP_TOOLS():
    
    def __init__(self, all_df, ):
        '''The Class used to analyze the one tail process's free energy convergence.
        Parameters
        ----------
        all_df: pd.DataFrame
            The dataframe obtained by alchemlyb extraction from amber prod.out
               
        Key properties
        ----------
        self.forward_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the forward estimate
        self.reverse_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the reverse estimate
        self.moving_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the moving estimate
        '''
        self.all_df = all_df
        self.all_df_fe_obj = CAL_FREE_ENERGY(self.all_df, False, 1,)
        ori_index_list = list(self.all_df.index.names)
        ori_index_list.pop(0)
        self.lambda_list = ori_index_list
        self.forward_esti_obj = TIME_SERIAL_DATA('forward')
        self.reverse_esti_obj = TIME_SERIAL_DATA('reverse')
        self.moving_esti_obj = TIME_SERIAL_DATA('moving')
  

    def generate_data_dict(self, wanted_win_lst, scale_f, ifsubsample=False):
        '''Apply the specific window list and scaling factor to get the data dict that used for time_serial analysis or free energy calculation
        Parameters
        ----------
        wanted_win_lst: list
            A list assigns simulation windows' data used in the following calculation
        scale_f: float
            To tail the percentage of the every dataframe data used for the calculation
        ifsubsample: bool
            To determine if subsampleCorraltedData.
        
        Return
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        '''
        ###apply scaling###
        tail_df_list = [j.tail(n=math.floor(j.shape[0]*scale_f)) for i,j in self.all_df.groupby(self.lambda_list, sort=False)]
        tail_df = pd.concat(tail_df_list,axis=0)
        # print(tail_df)
        groupobj = tail_df.groupby(self.lambda_list, sort=False)
        # print(groupobj)
        ###
        ###apply partial lambda###
        lambda_dict = {}
        K=0
        if wanted_win_lst is False:
            for i,j in groupobj:
                every_dataframe = pd.DataFrame(j)
                if ifsubsample:
                    every_dataframe = self.subsample_onestate_data(every_dataframe, i)
                    print(len(every_dataframe))
                lambda_dict[K] = every_dataframe
                K+=1
        else:
            for i,j in groupobj:
                # print(type(i),type(j))
                if wanted_win_lst[K]==i:
                    every_dataframe = pd.DataFrame(j)
                    if ifsubsample:
                        every_dataframe = self.subsample_onestate_data(every_dataframe, i)
                        print(len(every_dataframe))
                    lambda_dict[K] = every_dataframe
                    K+=1
                    if K == len(wanted_win_lst):
                        break
                else:
                    pass
        ###
        lambda_dict = lambda_dict
        return lambda_dict

    def subsample_onestate_data(self, onestate_df, column_value,):
        '''Use pymbar.timeseries.subsampleCorrelatedData to generate the uncorrelated delta_U data
        Parameters
        ----------
        onestate_df: pd.DataFrame
            The dataframe containing the U values of many states that sampled at a specific state.
        column_value:
            The column value of the onestate_df, which identifies the specific state.
        '''
        state_index = list(onestate_df.columns).index(column_value)
        A_n = np.array(onestate_df.iloc[:, state_index])
        a_indices = timeseries.subsampleCorrelatedData(A_n)
        if len(a_indices) < 100:
            new_onestate_df = onestate_df
            print('Warning: the subsample data len number is less than 100, will not update the subsample data.')
        else:
            new_onestate_df = onestate_df.iloc[a_indices, :]
        return new_onestate_df
    
    def plot_time_serial(self, png_prefix, lambda_, get_which_two_columns, plot_plan, ):#弃用
        '''Using the PLOTTING object to plot the time serial analysis for specifc lambda_ (lambda_ may be 'sum_of_all' or '0.0 to 1.0')
        Parameters
        ----------
        png_prefix: str
            The string that given for the output png's name, usually contain the information about the target, run-1 or run-2, which pair(38-60), which tail(cM2A)
        lambda_: str or float
            To specify which window for analysis
            For EXP(FEP) estimator, lambda_ can be one of [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 'sum_of_all']
            For BAR estimator, lambda_ can be one of ['0.0 to 0.05', '0.05 to 0.1', '0.1 to 0.2', '0.2 to 0.3', '0.3 to 0.4', '0.4 to 0.5', '0.5 to 0.6', 
            '0.6 to 0.7', '0.7 to 0.8', '0.8 to 0.9', '0.9 to 0.95', '0.95 to 1.0', '0.0 to 1.0']
        get_which_two_columns: [str, str]
            To specify which free energy and std columns for analysis
            For EXP(FEP) estimator, get_which_two_columns can be ['FEP_forward_bythislambda(kcal/mol)', 'FEP_forward_Rolling_STD(kcal/mol)'] 
            or ['FEP_reverse_bythislambda(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            For BAR estimator, get_which_two_columns will be ['free_energy(kcal/mol)', 'Rolling_STD(kcal/mol)']
        '''
        plot_obj = PLOTTING()
        fe_columns_name = get_which_two_columns[0]
        std_columns_name = get_which_two_columns[1]
#         print(self.forward_esti_obj.df_dict_with_std)
        if plot_plan == 2:
            for_df = self.forward_esti_obj.df_dict_with_std[lambda_]
            rev_df = self.reverse_esti_obj.df_dict_with_std[lambda_]
            forward_two_col_df = pd.concat([for_df[fe_columns_name], for_df[std_columns_name]],axis=1)
            reverse_two_col_df = pd.concat([rev_df[fe_columns_name], rev_df[std_columns_name]],axis=1)
            fe_he_std_dir = {'forward':forward_two_col_df, 
                             'reverse':reverse_two_col_df,
                             'fe': None,
                             'fe_std': None}
            plot_obj.plot_fe_time_serial('{}_{}_2.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        elif plot_plan == 3:
            mov_df = self.moving_esti_obj.df_dict_with_std[lambda_]
            for_df = self.forward_esti_obj.df_dict_with_std[lambda_]
            rev_df = self.reverse_esti_obj.df_dict_with_std[lambda_]
            forward_two_col_df = pd.concat([for_df[fe_columns_name], for_df[std_columns_name]],axis=1)
            reverse_two_col_df = pd.concat([rev_df[fe_columns_name], rev_df[std_columns_name]],axis=1)
            moving_two_col_df = pd.concat([mov_df[fe_columns_name], mov_df[std_columns_name]], axis=1)
            fe_he_std_dir = {'forward':forward_two_col_df, 
                             'reverse':reverse_two_col_df,
                             'moving':moving_two_col_df,
                             'fe': None,
                             'fe_std': None}
            plot_obj.plot_fe_time_serial('{}_{}_3.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        elif plot_plan == 1:
            mov_df = self.moving_esti_obj.df_dict_with_std[lambda_]
            moving_two_col_df = pd.concat([mov_df[fe_columns_name], mov_df[std_columns_name]], axis=1)
            fe_he_std_dir = {'moving':moving_two_col_df, 
                             'fe': None,
                             'fe_std': None}
            plot_obj.plot_fe_time_serial('{}_{}_1.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        else:
            print('Error! No such plot plan!')
        
    
    def lambda_dict_to_df(self, lambda_dict):
        '''Convert the lambda_dict to lambda_df, which is the key input for the CAL_FREE_ENERGY object.
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        
        Return 
        ----------
        processed_df: pd.DataFrame
            The dataframe generated by concating every windows data 
        '''
        
        df_list = []
        for single_df in lambda_dict.values():
            df_list.append(single_df)
        processed_df = pd.concat(df_list, )
        return processed_df
    
    def single_interval_estimate(self, lambda_dict, start_frame, end_frame, iffep=True, ):
        '''Calculate free energy according to the given lambda_dict, time_interval of [start_frame, end_frame]
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        start_frame: int
            The integer (actually the index) to represent the start frame of the trajactory
        end_frame: int
            The integer (actually the index) to represent the end frame of the trajactory
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: True
            
        Return
        ----------
        result_fe: pd.DataFrame
            For iffep==True: result_fe.index may be like [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 'sum_of_all'] with index name of 'lambda_value' 
                             result_fe.columns be like ['FEP_forward_bythislambda(kcal/mol)', 'FEP_reverse_bythislambda(kcal/mol)']
            For iffep==False:result_fe.index may be like ['0.0 to 0.05', '0.05 to 0.1', '0.1 to 0.2', '0.2 to 0.3', '0.3 to 0.4', '0.4 to 0.5', '0.5 to 0.6', 
            '0.6 to 0.7', '0.7 to 0.8', '0.8 to 0.9', '0.9 to 0.95', '0.95 to 1.0', '0.0 to 1.0'] with index name of 'delta_A_what_to_what'
                             result_fe.columns be like ['free_energy(kcal/mol)', 'estimated std']
        '''
        df_list = []
        for single_df in lambda_dict.values():
            single_df = single_df
            new_single_df = single_df.iloc[start_frame: end_frame,:]
            df_list.append(new_single_df)
        # processed_df = pd.concat(df_list, sort=True) # fixed no need sort
        processed_df = pd.concat(df_list,)
#         return processed_df
#         processed_df = self.lambda_dict_to_df(lambda_dict)

        fep_obj = CAL_FREE_ENERGY(processed_df, False, 1,)
        if iffep:
            result_fe = fep_obj.cal_FE_FEP(None, 'kcal/mol')
        else:
            result_fe = fep_obj.cal_FE(None, 'kcal/mol')
        return result_fe
        
   
    def moving_estimate(self, lambda_dict, std_step_ratio, divided_ratio, width_ratio = 0.2, iffep=False):
        '''Check the time serial convergence by moving estimate using the time interval with fixed width to scan the time serial data forward. 
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        width_ratio: float
            To determine the width of the fixed width of the moving time interval
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
        
        Update
        ----------
        self.moving_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the moving estimate
        '''
        static_ratio = divided_ratio
        count_=0
        while(divided_ratio <= 1):
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                end_frame = max(int(np.floor(frame_num*divided_ratio)),1)
                width = int(np.floor(frame_num*width_ratio))
                start_frame = end_frame - width
                if start_frame < 0:
                    start_frame = 0
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)
            # print(processed_df)
            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1, True)

            if iffep:
                data_ = fep_obj.cal_FE_FEP(None, 'kcal/mol')
            else:
                data_ = fep_obj.cal_FE(None, 'kcal/mol')

            if count_ == 0:
                lambda_seq = list(data_.index)
                columns = list(data_.columns)
                time_serial_data_obj = TIME_SERIAL_DATA('moving', lambda_seq, columns)
            else:
                pass
            time_serial_data_obj.update_data(divided_ratio,data_)
            # print(f'this time: {divided_ratio}')
            # print(time_serial_data_obj.data_dict)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        time_serial_data_obj.generate_df()

        time_serial_data_obj.cal_rolling_std(std_step_ratio,iffep)

        self.moving_esti_obj = time_serial_data_obj
        

    def forward_estimate(self, lambda_dict, std_step_ratio, divided_ratio=0.1, iffep=False, ifplottime_serials=True, ifplotdG_dlambda=False, ifuseforward=True):
        '''Check the time serial convergence by accumulating data points from the first frame forward with a fixed length
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
        
        Update
        ----------
        self.forward_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the forward estimate
        '''
        static_ratio = divided_ratio
        # print(divided_ratio)
        count_=0
        while(divided_ratio <= 1):
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                # print(f'{frame_num} frame_num----------------------------')
                start_frame = 0
                # print(frame_num*divided_ratio)
                end_frame = max(int(np.floor(frame_num*divided_ratio)),1)
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                # print(start_frame)
                # print(end_frame)
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)
            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1, True)

            if iffep:
                data_ = fep_obj.cal_FE_FEP(None, 'kcal/mol')
                if ifplotdG_dlambda:
                    get_csv_obj = dG_dlambda(data_, 'FEP')
                    if ifuseforward:
                        lambda_col_name = 'delta_what_to_what_forward'
                        fe_col_name = 'FEP_forward_bythislambda(kcal/mol)'
                    else:
                        lambda_col_name = 'delta_what_to_what_reverse'
                        fe_col_name = 'FEP_reverse_bythislambda(kcal/mol)'
                    get_csv_obj.get_dG_dlambda(lambda_col_name, fe_col_name)
                    get_csv_obj.plot_dG_dlambda(get_csv_obj.df_csv,f'FEP_dG_dlambda_{divided_ratio}')
            else:
                # print(f'ana_----------------{divided_ratio}')
                data_ = fep_obj.cal_FE(None, 'kcal/mol')
                if ifplotdG_dlambda:
                    data_for_dG_dlambda=pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kcal/mol)'], index=data_.index)
                    data_for_dG_dlambda['delta_A_what_to_what'] = data_.index
                    data_for_dG_dlambda['free_energy(kcal/mol)'] = data_['free_energy(kcal/mol)']
                    get_csv_obj = dG_dlambda(data_for_dG_dlambda, 'BAR')
                    get_csv_obj.get_dG_dlambda('delta_A_what_to_what', 'free_energy(kcal/mol)')
                    get_csv_obj.plot_dG_dlambda(get_csv_obj.df_csv,f'BAR_dG_dlambda_{divided_ratio}')
                # print(f'----------------{divided_ratio}-------------')
            if ifplottime_serials:
                if count_ == 0:
                    lambda_seq = list(data_.index)
                    columns = list(data_.columns)
                    time_serial_data_obj = TIME_SERIAL_DATA('forward', lambda_seq, columns)
                else:
                    pass
                time_serial_data_obj.update_data(divided_ratio, data_)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        if ifplottime_serials:
            time_serial_data_obj.generate_df()
            time_serial_data_obj.cal_rolling_std(std_step_ratio, iffep)

            self.forward_esti_obj = time_serial_data_obj
    
    def reverse_estimate(self, lambda_dict, std_step_ratio, divided_ratio=0.1, iffep=False):
        '''Check the time serial convergence by accumulating data points from the last frame backward with a fixed length
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
            
        Update
        ----------
        self.reverse_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the reverse estimate
        '''
        static_ratio = divided_ratio
        count_=0
        while(divided_ratio <= 1):
            start_ratio = round(1 - divided_ratio, 5)
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                end_frame = frame_num
                start_frame = int(np.floor(frame_num*start_ratio))
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)

            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1, True)

            if iffep:
                data_ = fep_obj.cal_FE_FEP(None, 'kcal/mol')
            else:
                data_ = fep_obj.cal_FE(None, 'kcal/mol')

            if count_ == 0:
                lambda_seq = list(data_.index)
                columns = list(data_.columns)
                time_serial_data_obj = TIME_SERIAL_DATA('reverse', lambda_seq, columns)
            else:
                pass
            time_serial_data_obj.update_data(divided_ratio, data_)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        time_serial_data_obj.generate_df()
        time_serial_data_obj.cal_rolling_std(std_step_ratio, iffep)

        self.reverse_esti_obj = time_serial_data_obj

        
    def use_fep_check_time_serial(self, lambda_dict, plot_plan, png_prefix, use_forward=True, std_step_ratio=0.1, divided_ratio=0.01, block_width_ratio_formoving=0.2):
        '''Use fep calculated free energy value to check time serial convergence, which will output the figure and update self.forward_esti_obj, self.reverse_esti_obj, self.moving_esti_obj
        '''
        self.moving_estimate(lambda_dict, std_step_ratio, divided_ratio, block_width_ratio_formoving, True)
        self.forward_estimate(lambda_dict, std_step_ratio, divided_ratio, True)
        self.reverse_estimate(lambda_dict, std_step_ratio, divided_ratio, True)
        lambda_lst = self.forward_esti_obj.df_dict_with_std.keys()
        if use_forward:
            get_which_two_columns=['FEP_forward_bythislambda(kcal/mol)', 'FEP_forward_Rolling_STD(kcal/mol)']
            png_prefix = png_prefix+'_useForward'
            for lambda_ in lambda_lst:
                if lambda_ == 1.0:
                    pass
                else:
                    self.plot_time_serial(png_prefix, lambda_, get_which_two_columns, plot_plan, )
        else:
            get_which_two_columns=['FEP_reverse_bythislambda(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            png_prefix = png_prefix+'_useReverse'
            for lambda_ in lambda_lst:
                if lambda_ == 0.0:
                    pass
                else:
                    self.plot_time_serial(png_prefix, lambda_, get_which_two_columns, plot_plan, )
    
    def use_bar_check_time_serial(self, ):
        '''Use fep calculated free energy value to check time serial convergence, which will update self.forward_esti_obj, self.reverse_esti_obj, self.moving_esti_obj
        '''
        lambda_dict = self.generate_data_dict(False, 1)
        self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, False)
        self.forward_estimate(lambda_dict, 0.1, 0.01, False)
        self.reverse_estimate(lambda_dict, 0.1, 0.01, False)
        
                            
    def output_fe_and_std(self, lambda_dict, unit, fe_mode, std_mode, ifdetail=False, ifplot=False):
        '''Output the calculated free energy and standard deviation according to the given data dict
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        fe_mode: str
            'BAR' or 'FEP';
            'BAR': Use Bennett Acceptance Ratio estimator to calculate free energy
            'FEP': Use EXPonential averaging (EXP) estimator to calculate free energy
        std_mode: str
            'time_serial' or 'bar_std' or 'bootstrap_std'
            'time_serial': Output the standard deviation calculated by the last 100*scale_f% moving estimated free energies' std
            'bar_std': Output the standard deviation calculated by BAR estimator. Note: the bar_std may underesitmate the std
            'bootstrap_std': Output the standard deviation calculated by bootstrapping. Resample 200 times and resample size will be the same as the input data.
         '''
#         print(f'lambda_dict: {lambda_dict}')
        
        # lambda_unk_df = pd.concat(list(lambda_dict.values()),sort=True) # fixed no need sort
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        cal_fe_obj = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        if ifplot:
            cal_fe_obj.plot_simulation_deltaU()
        if fe_mode == 'BAR':
            if std_mode == 'time_serial':
                data_df = cal_fe_obj.cal_FE(None, unit,)
                self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, False)
                if ifdetail:
                    fe = data_df.iloc[:,0]
                    lambda_timestd_list = []
                    for lambda_key in self.moving_esti_obj.df_dict_with_std.keys():
                        lambda_timestd_list.append(self.moving_esti_obj.df_dict_with_std[lambda_key].loc[0:1, f'free_energy({unit})'].std()) 
                    std = lambda_timestd_list    
                else:
                    fe = data_df.iloc[:,0][-1]
                    key_for_overall_fe = [i for i in self.moving_esti_obj.df_dict_with_std.keys()][-1] 
                    std = self.moving_esti_obj.df_dict_with_std[key_for_overall_fe].loc[0:1, f'free_energy({unit})'].std()
            else:
                if std_mode == 'bootstrap_std':
                    ifbootstrap = True
                elif std_mode == 'bar_std':
                    ifbootstrap = False
                data_df = cal_fe_obj.cal_FE(None, unit, ifbootstrap_std=ifbootstrap)
                if ifdetail:
                    fe = data_df.iloc[:,0]
                    std = data_df.iloc[:,1]
                else:
                    fe = data_df.iloc[:,0][-1]
                    std = data_df.iloc[:,1][-1]   
        elif fe_mode == 'FEP': 
            print('Note that the FEP so far calculation mode may not calculate the whole process free energy difference, only for convergence test! The ifdetail will not affect this fe_mode calcualtion!')
            fe = cal_fe_obj.cal_FE_FEP(None, unit).iloc[:,0][-1]
            self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, True)
#             print(self.moving_esti_obj.df_dict_with_std['sum_of_all'].loc[:scale_f, 'FEP_forward_bythislambda(kcal/mol)'])
            std = self.moving_esti_obj.df_dict_with_std['sum_of_all'].loc[0:1, f'FEP_forward_bythislambda({unit})'].std()
        if ifdetail:
            resu_df = pd.DataFrame(fe)
            resu_df[std_mode] = std
            return resu_df
        else:
            return fe, std            

    def moving_reweight_estimate(self, lambda_dict, divided_ratio, width_ratio = 0.2, wins_step=0.1, wins_num=3, forward=True, ifdiagnum=True):
        static_ratio = divided_ratio
        count_=0
        while(divided_ratio <= 1):
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                end_frame = int(np.floor(frame_num*divided_ratio))
                width = int(np.floor(frame_num*width_ratio))
                start_frame = end_frame - width
                if start_frame < 0:
                    start_frame = 0
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)

            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1,)
            rewei_data, diff_rewei, diff_percent_rewei = fep_obj.cal_reweight_all(use_wins_step=wins_step, use_wins_num=wins_num, unit='kcal/mol', forward=forward, ifdiagnum=ifdiagnum)
            if count_ == 0:
                diff_time_serial_data_obj = TIME_SERIAL_DATA('moving', list(diff_rewei.index), list(diff_rewei.columns))
                diff_percent_time_serial_data_obj = TIME_SERIAL_DATA('moving', list(diff_rewei.index), list(diff_rewei.columns))
            else:
                pass
            diff_time_serial_data_obj.update_data_reweight(divided_ratio,diff_rewei)
            diff_percent_time_serial_data_obj.update_data_reweight(divided_ratio,diff_percent_rewei)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        # diff_time_serial_data_obj.generate_df_reweight()
        # diff_percent_time_serial_data_obj.generate_df_reweight()

        self.diff_rewei_moving_esti_obj = diff_time_serial_data_obj
        self.diff_percent_moving_esti_obj = diff_percent_time_serial_data_obj

    def forward_reweight_estimate(self, lambda_dict, divided_ratio=0.1, wins_step=0.1, wins_num=3, forward=True, ifdiagnum=True):
        static_ratio = divided_ratio
        count_=0
        while(divided_ratio <= 1):
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                start_frame = 0
                end_frame = int(np.floor(frame_num*divided_ratio))
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)

            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1,)
            rewei_data, diff_rewei, diff_percent_rewei = fep_obj.cal_reweight_all(use_wins_step=wins_step, use_wins_num=wins_num, unit='kcal/mol', forward=forward, ifdiagnum=ifdiagnum)
            if count_ == 0:
                diff_time_serial_data_obj = TIME_SERIAL_DATA('forward', list(diff_rewei.index), list(diff_rewei.columns))
                diff_percent_time_serial_data_obj = TIME_SERIAL_DATA('forward', list(diff_rewei.index), list(diff_rewei.columns))
            else:
                pass
            diff_time_serial_data_obj.update_data_reweight(divided_ratio,diff_rewei)
            diff_percent_time_serial_data_obj.update_data_reweight(divided_ratio,diff_percent_rewei)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        # diff_time_serial_data_obj.generate_df_reweight()
        # diff_percent_time_serial_data_obj.generate_df_reweight()

        self.diff_rewei_forward_esti_obj = diff_time_serial_data_obj
        self.diff_percent_forward_esti_obj = diff_percent_time_serial_data_obj

    def reverse_reweight_estimate(self, lambda_dict, divided_ratio=0.1, wins_step=0.1, wins_num=3, forward=True, ifdiagnum=True):
        static_ratio = divided_ratio
        count_=0
        while(divided_ratio <= 1):
            start_ratio = round(1 - divided_ratio, 5)
            df_list = []
            for single_df in lambda_dict.values():
                frame_num = len(single_df)
                end_frame = frame_num
                start_frame = int(np.floor(frame_num*start_ratio))
                new_single_df = single_df.iloc[start_frame: end_frame,:]
                df_list.append(new_single_df)
            processed_df = pd.concat(df_list, sort=True)

            fep_obj = CAL_FREE_ENERGY(processed_df, False, 1,)
            rewei_data, diff_rewei, diff_percent_rewei = fep_obj.cal_reweight_all(use_wins_step=wins_step, use_wins_num=wins_num, unit='kcal/mol', forward=forward, ifdiagnum=ifdiagnum)
            if count_ == 0:
                diff_time_serial_data_obj = TIME_SERIAL_DATA('reverse',list(diff_rewei.index), list(diff_rewei.columns))
                diff_percent_time_serial_data_obj = TIME_SERIAL_DATA('reverse', list(diff_rewei.index), list(diff_rewei.columns))
            else:
                pass
            diff_time_serial_data_obj.update_data_reweight(divided_ratio,diff_rewei)
            diff_percent_time_serial_data_obj.update_data_reweight(divided_ratio,diff_percent_rewei)
            divided_ratio = round(divided_ratio+static_ratio, 5)
            count_+=1
        # diff_time_serial_data_obj.generate_df_reweight()
        # diff_percent_time_serial_data_obj.generate_df_reweight()

        self.diff_rewei_reverse_esti_obj = diff_time_serial_data_obj
        self.diff_percent_reverse_esti_obj = diff_percent_time_serial_data_obj

    def reweight_check_time_serial(self, lambda_dict, plot_plan, png_prefix, divided_ratio=0.01, width_ratio=0.2, wins_step=0.1, wins_num=3, forward=True, ifdiagnum=True, error_max=2):
        def plot_(time_estimate_obj,direction, df_type):
            df_dict = time_estimate_obj.df_dict
            self.plot_reweight_time_serial(df_dict, f'{png_prefix}_dG_{df_type}_rewei_use{wins_step}wins{wins_num}_{direction}', error_max)
        if 'moving' in plot_plan:
            self.moving_reweight_estimate(lambda_dict, divided_ratio, width_ratio, wins_step, wins_num, forward, ifdiagnum)
            plot_(self.diff_rewei_moving_esti_obj, 'Moving', 'diff')
            plot_(self.diff_percent_moving_esti_obj, 'Moving', 'diff_precent')
        if 'forward' in plot_plan:
            self.forward_reweight_estimate(lambda_dict, divided_ratio, wins_step, wins_num, forward, ifdiagnum)
            plot_(self.diff_rewei_forward_esti_obj, 'Forward', 'diff')
            plot_(self.diff_percent_forward_esti_obj, 'Forward', 'diff_precent')
        if 'reverse' in plot_plan:
            self.reverse_reweight_estimate(lambda_dict, divided_ratio, wins_step, wins_num, forward, ifdiagnum)
            plot_(self.diff_rewei_reverse_esti_obj, 'Reverse', 'diff')
            plot_(self.diff_percent_reverse_esti_obj, 'Reverse', 'diff_precent')


    def plot_reweight_time_serial(self, dict_time_reweight, png_prefix, error_max):
        images = []
        plot_obj = PLOTTING()
        for key, value in dict_time_reweight.items():
            pngfile = f'{png_prefix}{key}.png'
            plot_obj.plot_heatmap_cmap(value, error_max, pngfile)#df, error_max=2, png_file=None
            image = Image.open(pngfile)
            images.append(image)
        output_file = f'{png_prefix}.gif'
        images[0].save(output_file, save_all=True, append_images=images[1:], duration=500, loop=1)


    def output_fe_reweighting(self, lambda_dict, use_wins_step, use_wins_num, unit, forward, ifdiagnum, postfix, ifplot, error_max):
        lambda_unk_df = pd.concat(list(lambda_dict.values()))
        cal_rewei_obj = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
        result_df = cal_rewei_obj.cal_reweight_all(use_wins_step, use_wins_num, unit, forward, True, True, ifdiagnum, postfix, ifplot, error_max)
        return result_df

    def check_part_lambda(self, part_lambda_list, scale_f, fe_mode, std_mode, filename):
        ###generate lambda_dict
        lambda_dict_list = []
        for part_lambda in part_lambda_list:
            lambda_dict_list.append(self.generate_data_dict(part_lambda, scale_f))
#         print(lambda_dict_list)
        fe_data_df_list = []
        for lambda_dict_ in lambda_dict_list:
            fe, std = self.output_fe_and_std(lambda_dict_, scale_f, fe_mode, std_mode,)
            fe_data_df_list.append([fe, std])
        fe_data_df = pd.DataFrame(fe_data_df_list, )
        fe_data_df.columns = ['fe', 'std']
        fe_data_df.index = [str(i) for i in part_lambda_list]
        fe_data_df.index.names = ['lambda_scheme']
        if filename:
            fe_data_df.to_csv(filename)
            
        return fe_data_df
            
    
    def overlap_matrix_each_pair(self,  scale_f):
        
        for i in range(self.all_df_fe_obj.lambda_range-1):
            win_lst = [self.all_df_fe_obj.simulation_lambda[i], self.all_df_fe_obj.simulation_lambda[i+1]]
            self.overlap_matrix_for_specific_pair(scale_f, win_lst)
    
    def overlap_matrix_for_specific_pair(self, scale_f, win_lst):    
        print(win_lst)
        adjacent_fe_obj = CAL_FREE_ENERGY(self.all_df, win_lst, scale_f)
        overlap_matrix = adjacent_fe_obj.plot_overlap_matrix()
        adj_overlap_matrix = np.array([[overlap_matrix.loc[win_lst[0],win_lst[0]], overlap_matrix.loc[win_lst[1],win_lst[0]]], 
                                       [overlap_matrix.loc[win_lst[0],win_lst[1]], overlap_matrix.loc[win_lst[1],win_lst[1]]]])
        print(adj_overlap_matrix)

    def cal_forward_var(self, lambda_0, lambda_1, divided_ratio):
        dU_1to0_at_0 = self.all_df_fe_obj.get_deltaU_in_lambda(lambda_0, (lambda_1, lambda_0), filter_=False)
        dU_1to0_at_1 = self.all_df_fe_obj.get_deltaU_in_lambda(lambda_1, (lambda_1, lambda_0), filter_=False)
        frame_num = len(dU_1to0_at_0)
        step_width = int(np.floor(frame_num*divided_ratio))
        end_frame_list = [i for i in range(0, frame_num, step_width)]
        end_frame_list.pop(0)
        end_frame_list.append(frame_num)
        start_frame=0
        stat_df = pd.DataFrame(columns=['ui_tc', 'ui_bc', 'uii', 'ma','UI_deduce'])
        # uii_df = pd.DataFrame(columns=["uii"])
        # ma_df = pd.DataFrame(columns=["ma"])
        for end_frame in end_frame_list:
            # print(f'end_frame: {end_frame}')
            dU_at_0 = dU_1to0_at_0[start_frame:end_frame]
            dU_at_1 = dU_1to0_at_1[start_frame:end_frame]
            df, dff = BAR_(dU_at_0, -dU_at_1, method='self-consistent-iteration', 
                          maximum_iterations=1000, verbose=False)
            a, b = 1/2, 1/2
            tc = 1/(a + b*np.exp(-dU_at_1 + df))
            bc = 1/(a*np.exp(dU_at_0 - df) + b)
            UI_tc = tc.mean()
            UI_deduce = 1- UI_tc
            UI_bc = bc.mean()
            # print(UI_tc,UI_bc)
            UII = a*((tc**2).mean()) + b*((bc**2).mean())
            ma = (UI_tc - UII)/UI_tc
            stat_df = pd.concat([stat_df, pd.DataFrame(data=[[UI_tc, UII, ma, UI_deduce]], index=[end_frame, ], columns=['ui', 'uii', 'ma', 'UI_deduce'])])

            # print('__________',UII,'_____________')
            # print(ma)
            # print(ma_df)      

        self.stat_df = stat_df
        # return ui_df, uii_df

    def plot(self, lambda_0, lambda_1, divided_ratio=0.05):
        self.cal_forward_var(lambda_0, lambda_1, divided_ratio)
        plt.figure()
        # print(self.ui_df.index())
        x1 = list(self.stat_df.index)
        # print(x1)
        y1 = self.stat_df['ui']
        # x2 = list(self.uii_df.index)
        y2 = self.stat_df["uii"]
        y3 = self.stat_df['ma']
        y4 = self.stat_df['UI_deduce']
        plt.plot(x1, y1, "r")
        plt.plot(x1, y2, "b")
        plt.plot(x1, y3, 'g')
        plt.plot(x1, y4, 'black')
        plt.savefig(f'ma_{lambda_0}_to_{lambda_1}.png')
        plt.clf()
        plt.close()
