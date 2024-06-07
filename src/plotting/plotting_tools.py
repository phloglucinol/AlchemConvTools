import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
matplotlib.use('agg')
from contextlib import contextmanager





@contextmanager
def plot_settings(figsize, font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=None, y_label=None, y_lim_tuple=None, x_lim_tuple=None, png_file_name=None, iflegend=False):
    plt.clf()
    figure, ax = plt.subplots(figsize=figsize)
    plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
    plt.rcParams['font.sans-serif']=[font_family]#设置全局字体

    # 设置坐标轴刻度粗细及字体大小
    plt.minorticks_on()#开启小坐标
    plt.tick_params(which='major',width=3, length=6)#设置大刻度的大小
    plt.tick_params(which='minor',width=2, length=4)#设置小刻度的大小
    plt.yticks(fontproperties=font_family, size=font_size, weight=font_weight)#设置y轴刻度字体大小及加粗
    plt.xticks(fontproperties=font_family, size=font_size, weight=font_weight)#设置x轴刻度字体大小及加粗

    # 设置x轴和y轴的标签
    if x_label is not None:
        plt.xlabel(x_label, fontsize=font_size, fontweight=font_weight)
    if y_label is not None:
        plt.ylabel(y_label, fontsize=font_size, fontweight=font_weight)

    # 设置坐标轴粗细
    ax=plt.gca();#获得坐标轴的句柄
    bottom_lineweight = lineweight_dict['bottom']
    left_lineweight = lineweight_dict['left']
    right_lineweight = lineweight_dict['right']
    top_lineweight = lineweight_dict['top']
    ax.spines['bottom'].set_linewidth(bottom_lineweight);#设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(left_lineweight);#设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(right_lineweight);#设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(top_lineweight);#设置上部坐标轴的粗细

    yield # 在这里执行plotting操作

    # 设置坐标轴范围
    if y_lim_tuple is not None:
        plt.ylim(y_lim_tuple[0], y_lim_tuple[1])#限制y轴的大小
    if x_lim_tuple is not None:
        plt.xlim(x_lim_tuple[0], x_lim_tuple[1])#限制x轴的大小

    # 设置图例
    if iflegend:
        plt.legend(loc="best",fontsize=font_size,scatterpoints=1,shadow=True,frameon=False)#添加图例，loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，shadow为图例边框阴影，frameon控制是否有边框

    plt.tight_layout()#防止因为ticklabel或者title等过长、或者过大导致的显示不全

    # 图片存储相关
    if png_file_name is not None:
        plt.savefig(png_file_name, format='png', bbox_inches='tight', transparent=True, dpi=600)
        # plt.show()
        plt.clf()
    else:
        plt.clf()


class PLOTTING():
    
    def __init__(self,):
        pass
    
    def unpack_esti_dir(self, esti_df):
        frame_ratio = np.around(np.array(esti_df.index), decimals=2)
        fe = np.array(esti_df.iloc[:,0])
        std = np.array(esti_df.iloc[:,1])
        fe_up = fe+std
        fe_down = fe-std
        #print(frame_ratio, fe, fe_up, fe_down)
        return frame_ratio, fe, fe_up, fe_down

    def plot_dU_distribution(self, d_u,png_file=None, ifplt=True, bins_=50):
        if ifplt:
            png_file_name = png_file
        else:
            png_file_name = None
        with plot_settings(figsize=(8, 6), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file_name, iflegend=False):
            count,bin_edges = np.histogram(d_u,bins=bins_,density=True)
            xu = []
            for i in range(len(bin_edges)-1):
                xu.append((bin_edges[i]+bin_edges[i+1])/2)
            plt.plot(xu,count,"o", color = "green", markersize=4)
        return count, xu

    def plot_weighted_dU_distribution(self, d_u, weights_array, png_file=None, ifplt=False, bins_=50):
        if ifplt:
            png_file_name = png_file
        else:
            png_file_name = None
        with plot_settings(figsize=(8, 6), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file_name, iflegend=False):
            count,bin_edges = np.histogram(d_u,bins=bins_,density=True, weights = weights_array)
            xu = []
            for i in range(len(bin_edges)-1):
                xu.append((bin_edges[i]+bin_edges[i+1])/2)
            plt.plot(xu,count,"o", color = "green", markersize=4)
        return count, xu

    def plot_resample_dU_distribution(self, ori_d_u, resample_dU, png_file=None, ifplt=True,bins=100):
        plt.clf()
        with plot_settings(figsize=(8, 6), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            pre_y,pre_bin = self.plot_dU_distribution(ori_d_u,png_file=None,ifplt=False,bins_=bins)
            resample_count, xu = self.plot_dU_distribution(resample_dU,png_file=None, ifplt=False,bins_=bins)
            plt.plot(pre_bin, pre_y,"o",label='original target dU distrib',color="green",markersize=5)
            plt.plot(xu,resample_count,"x",label='resampled target dU distrib',color="blue",markersize=5)
        # ori_d_u_mean, ori_d_u_std = np.mean(ori_d_u), np.std(ori_d_u)
        # resample_d_u_mean, resample_d_u_std = np.mean(resample_dU), np.std(resample_dU)
        # plt.text(union_x[-1],0.5*max(prob_y),'mean_pre_d_u : {:.3f}\nstd_pre_d_u : {:.5f}\nmean_resample_d_u : {:.3f}\nstd_resample_d_u : {:.5f}'.format\(ori_d_u_mean,ori_d_u_std,resample_d_u_mean,resample_d_u_std),\
        #          horizontalalignment='left',verticalalignment='center')#文本注释mean值和std值

    def plot_Boltzmann_weight_PdU_result(self, input_interpolate_bin, plot_input_interpolate_bin, B_wei_times_PdU, interpolate_PdU, plot_Boltzmann_weight, y_lim_tuple, png_file):
        with plot_settings(figsize=(8, 6), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU_bin(kBT)', y_label='P(dU) or Boltzman_weight or Boltzman_weight*P(dU)', y_lim_tuple=y_lim_tuple, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            plt.plot(input_interpolate_bin,B_wei_times_PdU,"o-", color = "tab:red", markersize=5,linewidth=3, label='Boltzman_weight*P(dU)')
            plt.plot(input_interpolate_bin,interpolate_PdU,"o-", color = "tab:blue", markersize=5,linewidth=3, label='P(dU)')
            plt.plot(plot_input_interpolate_bin,plot_Boltzmann_weight,"-", color = "black", linewidth=3, label='Boltzman_weight')

    def plot_cfm_checking_result(self, df, df_std, d_u_bin_filter_nan, diff_filter_nan, input_interpolate_x, png_file):
        xlabel = r'$\Delta U_{i,i+1}\/(k_{B}T)$'
        ylabel = r'$\Delta g_{i+1,i}\/(k_{B}T)$'
        with plot_settings(figsize=(8, 6), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=xlabel, y_label=ylabel, y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            plt.plot(d_u_bin_filter_nan,diff_filter_nan,"o-", color = "tab:red", markersize=5,linewidth=3)
            plt.fill_between(input_interpolate_x, df-df_std, df+df_std, color='tab:blue', zorder=-1)

    def plot_fe_time_serial(self, png_file_name, **fe_he_std_dir,):
        figsize=8,6
        figure, ax = plt.subplots(figsize=figsize)
        plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
        plt.rcParams['font.sans-serif']=['DejaVu Sans']#设置全局字体，可选择需要的字体替换掉‘Times New Roman’
        #plt.rcParams['font.sans-serif']=['SimHei']#使用黑体'SimHei'作为全局字体，可以显示中文
        font1={'family': 'DejaVu Sans', 'weight': 'bold', 'size': 14}#设置字体模板，
        font2={'family': 'DejaVu Sans', 'weight': 'bold', 'size': 20}#wight为字体的粗细，可选 ‘normal\bold\light’等
        font3={'family': 'DejaVu Sans', 'weight': 'light', 'size': 12}#size为字体大小
        plt.minorticks_on()#开启小坐标
        plt.tick_params(which='major',width=3, length=6)#设置大刻度的大小
        plt.tick_params(which='minor',width=2, length=4)#设置小刻度的大小
        ax=plt.gca();#获得坐标轴的句柄
        ax.spines['bottom'].set_linewidth(4);#设置底部坐标轴的粗细
        ax.spines['left'].set_linewidth(4);#设置左边坐标轴的粗细
        ax.spines['right'].set_linewidth(4);#设置右边坐标轴的粗细
        ax.spines['top'].set_linewidth(4);#设置上部坐标轴的粗细
        # print(fe_he_std_dir)
        ######
        if len(fe_he_std_dir) == 5:
            dir_move = fe_he_std_dir['moving']
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_move, fe_move, fe_up_move, fe_down_move = self.unpack_esti_dir(dir_move)
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            ###moving estimate plot
            plt.plot(frame_ratio_move, fe_move, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_move, fe_up_move, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_move, fe_down_move, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_move, fe_down_move, fe_up_move, where=fe_down_move <= fe_up_move,
                     facecolor='#a9f971', interpolate=True,alpha=0.5)
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_forw, fe_up_forw, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_forw, fe_down_forw, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_reve, fe_up_reve, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_reve, fe_down_reve, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)

        elif len(fe_he_std_dir) == 4:
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)

        elif len(fe_he_std_dir) == 3:
            aly_stra = list(fe_he_std_dir.keys())[0]
            dir_ = fe_he_std_dir[aly_stra]
            frame_ratio_, fe_, fe_up_, fe_down_ = self.unpack_esti_dir(dir_)
            y_min = fe_up_[-1]-2
            y_max = fe_up_[-1]+2
            ### estimate plot
            if aly_stra == 'forward':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            elif aly_stra == 'reverse':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a2cffe', interpolate=True, alpha=0.5)
            elif aly_stra == 'moving':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a9f971', interpolate=True,alpha=0.5)
        if fe_he_std_dir['fe'] is None:
            estimate_fe = None
        else:
            fe_x_plot = np.linspace(0,1,10000)
            estimate_fe = np.array([fe_he_std_dir['fe'] for i in range(0, len(fe_x_plot))])
            estimate_std = np.array([fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))])
            estimate_std_range_up = estimate_fe+estimate_std
            estimate_std_range_down = estimate_fe-estimate_std
            # estimate_std_range_up = [fe_he_std_dir['fe']+fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))]
            # estimate_std_range_down = [fe_he_std_dir['fe']-fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))]
            # print(fe_x_plot)
            # print(len(estimate_fe))
            # print(len(estimate_std_range_down))
            # print(estimate_std_range_up)
        ###estimate_fe_horizontal_line
        if estimate_fe is not None:
            plt.plot(fe_x_plot, estimate_fe, '-', lw=2, color='#FFE11A', label='BAR estimate result', alpha=1)
            plt.fill_between(fe_x_plot, estimate_std_range_down, estimate_std_range_up, where=estimate_std_range_down <= estimate_std_range_up,
                    facecolor='#FFE11A', interpolate=True, alpha=0.5)
        ######
        #plt.xlim(-15,-5)#限制x轴的大小
        # print('___________ymin:',y_min,"\n____________",y_max)
        plt.ylim(y_min,y_max)#限制y轴的大小
        #plt.title("Statistical analysis on the restraint strategies used in the test system ",fontdict=font2)#标题
        #plt.xlabel(r'$\Delta G_{exp} $ (kcal/mol)', fontdict=font2)#x轴标签
        #plt.ylabel(r'$\Delta G_{cal}^{MM-GBSA} $ (kcal/mol)',fontdict=font2)#y轴标签
        plt.legend(loc="best",scatterpoints=1,prop=font2,shadow=True,frameon=False)#添加图例，\
        # # loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，\
        # # shadow为图例边框阴影，frameon控制是否有边框
        plt.tick_params(\
            axis='x',#设置x轴
            direction='in',# 小坐标方向，in、out
            which='both',      # 主标尺和小标尺一起显示，major、minor、both
            bottom=True,      #底部标尺打开
            top=False,         #上部标尺关闭
            labelbottom=True, #x轴标签打开
            labelsize=20) #x轴标签大小
        plt.tick_params(\
            axis='y',
            direction='in',
            which='both',
            left=True,
            right=False,
            labelbottom=True,
            labelsize=20)
        plt.yticks(fontproperties='DejaVu Sans', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.xticks(fontproperties='DejaVu Sans', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.ticklabel_format(axis='both',style='sci')#sci文章的风格
        plt.tight_layout(rect=(0,0,1,1))#rect=[left,bottom,right,top]
        plt.tight_layout()#防止因为ticklabel或者title等过长、或者过大导致的显示不全
        plt.savefig(png_file_name,format="png",dpi=600,transparent=True)#存储png
        #plt.show()

    def plot_heatmap(self, df, vmin=-1, vmax=1, center=0, fontsize=10, png_file=None):
        
        # xLabel = df.columns
        # yLabel = df.index
        with plot_settings(figsize=(15,11), font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=None, y_label=None, y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=False):
            cmap = 'coolwarm'
        #     im = ax.imshow(df,cmap = 'gray' , interpolation = 'nearest',origin='lower')
            df = df.astype(float)
            # print(df)
            ax = sns.heatmap(data=df, linewidth=0.8, cmap=cmap, vmin=vmin, vmax=vmax, center=center)
            # ax.set_yticks(range(len(yLabel)))
            # ax.set_yticklabels(yLabel,fontsize=fontsize, rotation=0)
            # # ax.set_xticks(range(len(xLabel)))
            # ax.set_xticklabels(xLabel,fontsize=fontsize, rotation=45) 

    def plot_heatmap_cmap(self, df, error_max=2, png_file=None):
        with plot_settings(figsize=(15,11), font_family='DejaVu Sans', font_weight='bold', font_size=10, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=None, y_label=None, y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=False):
            df = df.dropna(how='all')
            df = df.astype(float)
            # 定义分段的边界值 'tab20c'
            diff_acept = np.around(error_max/df.shape[0], decimals=3)
            thresholds = [-2*diff_acept, -diff_acept, -1/2*diff_acept, 1/2*diff_acept, diff_acept, 2*diff_acept]
            norm = plt.Normalize(vmin=thresholds[0], vmax=thresholds[-1])
            ax = sns.heatmap(data=df, linewidth=0.8, cmap='coolwarm', cbar=True, norm=norm)
            cbar = ax.collections[0].colorbar
            cbar.set_ticks(thresholds)
            cbar.set_ticklabels(thresholds)

