import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
try:
    from glob import glob
except:
    pass



def col_val_count(df, column_name, min_value, max_value):
    count = ((df.loc[:,column_name] >= min_value) & (df.loc[:, column_name] <= max_value)).sum()
    return count

def plot_histogram(counts_dict, png_name):
    # plt.figure(num=1,figsize=(2.7,1.8))#创建画图，序号为1，图片大小为2.7*1.8
    # plt.figure(num=1,figsize=(6,6.5))#创建画图
    figsize=6,6
    figure, ax = plt.subplots(figsize=figsize)
    plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
    plt.rcParams['font.sans-serif']=['Times New Roman']
    #设置全局字体，可选择需要的字体替换掉‘Times New Roman’
    #使用黑体'SimHei'作为全局字体，可以显示中文
    #plt.rcParams['font.sans-serif']=['SimHei']
    font1={'family': 'Times New Roman', 'weight': 'bold', 'size': 14}#设置字体模板，
    font2={'family': 'Times New Roman', 'weight': 'light', 'size': 16}
    #wight为字体的粗细，可选 ‘normal\bold\light’等
    #size为字体大小
    plt.minorticks_on()#开启小坐标
    # plt.title("Statistical analysis on the restraint strategies used in the test system ",fontdict=font2)#标题
    
    label_list = np.array([ i for i in counts_dict.keys()])
    x_=np.arange(len(label_list))+1
    height_list = np.array([ i for i in counts_dict.values()])
    height_sum = np.sum(height_list)
    print(label_list)
    print(x_)
    print(height_list)
    rects1 = plt.bar(x=x_, height=height_list, width=0.45, color='tab:red',)
    plt.ylim(0, 100)

    # plt.ylabel("Count", fontdict=font1)
    # plt.xlabel('Overlap value range',fontdict=font1)
    # plt.legend(loc="best",scatterpoints=1,prop=font1,shadow=True,frameon=False)#添加图例，\
    # # loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，\
    # # shadow为图例边框阴影，frameon控制是否有边框
    # plt.xticks(x_,label_list,size='small',rotation=30)
    plt.xticks(x_,label_list,)
    # for a,b in zip(x_,height_list):
    #     plt.text(a, b+0.5, '%.0f' % b, ha='center', va= 'bottom',fontsize=14)
    # plt.tight_layout(rect=(0,0,1,1))#rect=[left,bottom,right,top]
    ###设置坐标轴的粗细
    ax=plt.gca();#获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(4);###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(4);####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(4);###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(4);####设置上部坐标轴的粗细
    #设置坐标刻度值的大小以及刻度值的字体
    plt.tick_params(labelsize=14)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    [label.set_fontweight('bold') for label in labels]

    zipped = zip(x_, height_list)
    for i, value in zipped:
        text_ = f'{value} {np.around(100*value/height_sum, 1)}%'
        plt.text(i, value, text_, ha='center', va='bottom')

    plt.tick_params(which='major',width=4, length=5)
    plt.tick_params(which='minor',width=3, length=3)
    plt.savefig(png_name,format="png",dpi=600, transparent=True)
    # plt.show()

if __name__ == '__main__':
    import sys
    png_name = sys.argv[1]
    file_pattern = './group*/fe_cal_out/pai_checking/pai_check.csv'
    files = glob(r'{}'.format(file_pattern))
    # print(files)
    df_list = []
    for i in files:
        df_list.append(pd.read_csv(i, ))

    df_all = pd.concat(df_list, axis=0)
    # print(df_all.columns)
    bounds_dict = {'less than 0':[-10000000, -0.000000001], '0 to 0.1':[0, 0.099999], '0.1 to 0.5':[0.1, 0.499999], 'greater than 0.5':[0.5, 1000000]}
    counts_dict = {'less than 0':0, '0 to 0.1': 0, '0.1 to 0.5': 0, 'greater than 0.5': 0}

    for bound_key in bounds_dict.keys():
        for col in ['relative_entrophy_0', 'relative_entrophy_1']:
            bound = bounds_dict[bound_key]
            counts_dict[bound_key]+=col_val_count(df_all, col, bound[0], bound[1])

    plot_histogram(counts_dict, png_name)

