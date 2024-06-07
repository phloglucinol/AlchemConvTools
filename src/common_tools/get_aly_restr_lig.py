import pandas as pd
import math 

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

if __name__ == '__main__':
    print(get_aly_restr_lig('res_databystd.csv'))