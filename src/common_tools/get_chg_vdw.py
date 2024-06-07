import pandas as pd 
from optparse import OptionParser
import sys

class optParser():  
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-d', '--dir', dest='file_directory', help='Directory in which free_ene.csv are stored. Default: current directory', default='./free_ene.csv')
        parser.add_option('-c', '--charge_wins', dest='charge_wins', help='To specify the windows number of the decharge process. Default: 11', default=11)
        parser.add_option('-r', '--restra_wins', dest='restra_wins', help='To specify the windows number of the restraining process (in complex). Default: 2', default=2)
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

if __name__ == '__main__':
    opts = optParser('')
    dir = str(opts.option.file_directory)
    charge_wins = int(opts.option.charge_wins)
    restra_wins = int(opts.option.restra_wins)
    df = pd.read_csv(dir, delimiter="|")

    if restra_wins == 0: 
        res_dA = 0
        chg_start_line = 0
    elif restra_wins >0:
        res_df = df.iloc[0: restra_wins-1, :]
        # print(res_df)
        res_dA = res_df['free_energy(kcal/mol)'].sum()
        chg_start_line = restra_wins-1
    else:
        print("restra_wins should not be less than zero.", file=sys.stderr)
        sys.exit()

    if charge_wins > 0:
        vdw_start_line = chg_start_line+charge_wins-1
        chg_df = df.iloc[chg_start_line: vdw_start_line, :]
        vdw_df = df.iloc[vdw_start_line:-1, :]
    else:
        print("charge_wins must be greater than zero.", file=sys.stderr)
        sys.exit()
    chg_dA = chg_df['free_energy(kcal/mol)'].sum()
    vdw_dA = vdw_df['free_energy(kcal/mol)'].sum()

    print(f"{chg_dA} {vdw_dA} {res_dA}")
