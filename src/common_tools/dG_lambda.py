import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
import matplotlib
matplotlib.use('agg')

class OptParser():
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-f', '--free_ene_file', dest='free_ene_file', help='The file name of the free energy calculation result.', default='free_ene.csv')
        # parser.add_option('-n', '--filename', dest='filename', help='The prefix of the result csv filename. Default:ts.csv', default='ts')
        parser.add_option('-o', '--output_png_file', dest='output_png_file', help='The file name of the output png file.', default='dG_lambda.png')
        parser.add_option('-l', '--delimiter', dest='delimiter', help='The delimiter of free_ene_file.', default='|')
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

class dG_lambda():
    def __init__(self, csv_file, delimiter='|'):
        self.csv_file = csv_file
        self.df_csv = pd.read_csv(self.csv_file, delimiter)
        self.df_csv = self.df_csv.head(self.df_csv.shape[0]-1)
        new_data = []
        cumulative_sum = 0
        for index, row in self.df_csv.iterrows():
            cumulative_sum += row[1]
            new_data.append([str(row[0].split(" to ")[1]), cumulative_sum])
        self.new_df = pd.DataFrame(new_data, columns=["lambda", "Cumulative_Free_Energy"])
    def plot_dG_lambda(self, png_file=None):
        figsize=8,6
        figure, ax = plt.subplots(figsize=figsize)
        plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
        plt.rcParams['font.sans-serif']=['Times New Roman']#设置全局字体，可选择需要的字体替换掉‘Times New Roman’
        #plt.rcParams['font.sans-serif']=['SimHei']#使用黑体'SimHei'作为全局字体，可以显示中文
        font1={'family': 'Times New Roman', 'weight': 'bold', 'size': 14}#设置字体模板，
        font2={'family': 'Times New Roman', 'weight': 'bold', 'size': 20}#wight为字体的粗细，可选 ‘normal\bold\light’等
        font3={'family': 'Times New Roman', 'weight': 'light', 'size': 12}#size为字体大小
        plt.minorticks_on()#开启小坐标
        plt.tick_params(which='major',width=3, length=6)#设置大刻度的大小
        plt.tick_params(which='minor',width=2, length=4)#设置小刻度的大小
        ax=plt.gca();#获得坐标轴的句柄
        ax.spines['bottom'].set_linewidth(4);#设置底部坐标轴的粗细
        ax.spines['left'].set_linewidth(4);#设置左边坐标轴的粗细
        ax.spines['right'].set_linewidth(4);#设置右边坐标轴的粗细
        ax.spines['top'].set_linewidth(4);#设置上部坐标轴的粗细
        
        x = []
        count_ = 0
        for i in self.new_df.iloc[:,0]:
            try:
                i_ = float(i)
                
            except:
                i_ = count_
            x.append(i_)
            count_+=1
        
        y = self.new_df.iloc[:,1]
        
        plt.plot(x, y, '^-', lw=2, color='#736AFF',)
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
        plt.yticks(fontproperties='Times New Roman', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.xticks(fontproperties='Times New Roman', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.ticklabel_format(axis='both',style='sci')#sci文章的风格
        plt.tight_layout(rect=(0,0,1,1))#rect=[left,bottom,right,top]
        plt.tight_layout()#防止因为ticklabel或者title等过长、或者过大导致的显示不全
        if png_file:
            plt.savefig(png_file,format="png",dpi=600,transparent=True)#存储png
        else:
            plt.show()

if __name__ == "__main__": 
    opts = OptParser('')
    free_ene_file = str(opts.option.free_ene_file)
    output_png_file = str(opts.option.output_png_file)
    delimiter = str(opts.option.delimiter)
    ts = dG_lambda(free_ene_file, delimiter)
    ts.plot_dG_lambda(output_png_file)

