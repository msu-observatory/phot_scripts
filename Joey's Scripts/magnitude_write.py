import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def mag(Tcounts, Terror, Ccounts, Cerror, Cmag):

    Tmag = -2.5*np.log10(Tcounts/Ccounts) + Cmag
    Terrmag = np.sqrt((2.5*Terror/(np.log(10.0)*Tcounts))**2 +
    ((2.5*Cerror/(np.log(10.0)*Ccounts))**2))

    return Tmag, Terrmag

def V_var(vvar, vcomp, bvar, bcomp, Vcomp):
    '''
    Calculates the color transformed values for V magnitude

    Parameters:
    vvar: variable star measured V magnitude
    vcomp: comparison star measured V magnitude
    TvBv, Tbv: color transforms
    bvar: variable star measured B magnitude
    bcomp: comparison star measured B magnitude
    Vcomp: published V magnitude of comparison star
    '''

    Tbv = 1.458
    Tvbv = -0.065

    deltav = vvar - vcomp

    deltabv = (bvar-vvar)-(bcomp-vcomp)

    deltaBV = Tbv * deltabv

    return deltav + Tvbv * deltaBV + Vcomp

def B_var(vvar, vcomp, bvar, bcomp, Bcomp):
    '''
    Calculates the color transformed values for B magnitude

    Parameters:
    vvar: variable star measured V magnitude
    vcomp: comparison star measured V magnitude
    TbBv, Tbv: color transforms
    bvar: variable star measured B magnitude
    bcomp: comparison star measured B magnitude
    Bcomp: published B magnitude of comparison star
    '''

    Tbv = 1.458

    Tbbv = 0.250

    deltab = bvar - bcomp

    deltabv = (bvar-vvar)-(bcomp-vcomp)

    deltaBV = Tbv * deltabv

    return deltab + Tbbv * deltaBV + Bcomp


class MagnitudeCalculation():
    '''
    A class for caculating magnitudes of stars. Takes data from an excel file
    (saved as a csv) output by AstroImageJ and creates lightcurves for each of
    the comparison stars, then outputs the calculated magnitudes in AAVSO and
    CBA formats (plots and data output still pending). Will also (eventually)
    have the ability to color correct the calculated magnitudes. Currently
    support for correcting V and B is being added.
    '''
    def __init__(self):
        self.numC = float(input('Number of comparison stars (up to 6):'))
        if self.numC > 6 or self.numC < 1:
            raise ValueError('Number of comparison stars must be between 1 and 6')
        self.star_name = input('Source name:')
        filter = input('input the filter type (Clear, R, V, I, or B): ').upper()
        if filter == 'CLEAR':
            self.filter = 'CV'
        else:
            self.filter = filter
        if self.filter == 'I' or self.filter == 'R' or self.filter == 'CV':
            self.trans = 'NO'
        else:
            self.trans = input('Transforming the data? (YES or NO)').upper()
        self.day  = str(input("Date of observations (Ex: 04Feb2018):"))
        self.mtype = 'STD'
        self.group = 'na'
        self.chart = 'na'
        self.notes = 'na'
        

    def comp_stars(self):

        self.C2name=str(input('Name of first comparison star:'))
        self.C2mag = float(input('Magnitude of the first comparison star:'))        
        if self.numC == 1:
            self.C3mag = 'none'
            self.C3name = 'none'
        if self.numC  > 1:
            self.C3name =str(input('Name of second comparison star:'))
            self.C3mag = float(input('Magnitude of the second comparison star:'))
            
        if self.numC > 2:
            self.C4name = str(input('Name of third comparison star:'))
            self.C4mag = float(input('Magnitude of the third comparison star:'))
            
        if self.numC > 3:
            self.C5name =str(input("Name of fourth comparison star:"))
            self.C5mag = float(input('Magnitude of the fourth comparison star:'))
            
        if self.numC > 4:
            self.C6name = str(input("Name of the fifth comparison star"))
            self.C6mag = float(input('Magnitude of the fifth comparison star:'))
            
        if self.numC >5:
            self.C7name = str(input("Name of the sixth comparison star"))
            self.C7mag = float(input('Magnitude of the sixth comparison star:'))
            


    def read_file(self):

        filename = input('input the filename (csv):')+".csv"

        data = pd.read_csv(filename)

        self.time = data['J.D.-2400000']
        self.jd = data['JD_UTC']
        self.airmass = data['AIRMASS']
        self.Tcounts = data['Source-Sky_T1']
        self.Terror = data['Source_Error_T1']
        self.C2counts = data['Source-Sky_C2']
        self.C2error = data['Source_Error_C2']
        if self.numC > 1:
            self.C3counts = data['Source-Sky_C3']
            self.C3error = data['Source_Error_C3']
        if self.numC > 2:
            self.C4counts = data['Source-Sky_C4']
            self.C4error = data['Source_Error_C4']
        if self.numC > 3:
            self.C5counts = data['Source-Sky_C5']
            self.C5error = data['Source_Error_C5']
        if self.numC > 4:
            self.C6counts = data['Source-Sky_C6']
            self.C6error = data['Source_Error_C6']
        if self.numC > 5:
            self.C7counts = data['Source-Sky_C7']
            self.C7counts = data['Source_Error_C7']


    def mag_calc(self):

        mag_list = []
        error_list = []

        Tmag2, Terrmag2 = mag(self.Tcounts, self.Terror, self.C2counts, self.C2error, self.C2mag)
        print(Tmag2[0], Terrmag2[0])
        mag_list.append(Tmag2)
        error_list.append(Terrmag2)

        if self.numC > 1:
            Tmag3, Terrmag3 = mag(self.Tcounts, self.Terror, self.C3counts, self.C3error, self.C3mag)
            print(Tmag3[0],Terrmag3[0])
            mag_list.append(Tmag3)
            error_list.append(Terrmag3)

        if self.numC > 2:
            Tmag4, Terrmag4 = mag(self.Tcounts, self.Terror, self.C4counts, self.C4error, self.C4mag)
            print(Tmag4[0],Terrmag4[0])
            mag_list.append(Tmag4)
            error_list.append(Terrmag4)

        if self.numC > 3:
            Tmag5, Terrmag5 = mag(self.Tcounts, self.Terror, self.C5counts, self.C5error, self.C5mag)
            print(Tmag5[0],Terrmag5[0])
            mag_list.append(Tmag5)
            error_list.append(Terrmag5)

        if self.numC > 4:
            Tmag6, Terrmag6 = mag(self.Tcounts, self.Terror, self.C6counts, self.C6error, self.C6mag)
            print(Tmag6[0],Terrmag6[0])
            mag_list.append(Tmag6)
            error_list.append(Terrmag6)

        if self.numC > 5:
            Tmag7, Terrmag7 = mag(self.Tcounts, self.Terror, self.C7counts, self.C7error, self.C7mag)
            print(Tmag7[0],Terrmag7[0])
            mag_list.append(Tmag7)
            error_list.append(Terrmag7)

        return mag_list, error_list


    def make_plots(self, mag_list, error_list):
        fig = plt.figure(figsize=(12, 12))

        p1 = fig.add_subplot(221)
        p1.errorbar(self.time, mag_list[0], yerr=error_list[0], marker='d', markersize =5, linestyle='')
        plt.xlabel('Time (JD-2400000)', fontsize=16)
        plt.ylabel('Target Mag (relative to C2)', fontsize=16)
        plt.gca().invert_yaxis()

        if self.numC > 1:
            p2 = fig.add_subplot(222)
            p2.errorbar(self.time, mag_list[1], yerr=error_list[1], marker='d', markersize =5, linestyle='')
            plt.xlabel('Time (JD-2400000)', fontsize=16)
            plt.ylabel('Target Mag (relative to C3)', fontsize=16)
            plt.gca().invert_yaxis()

        if self.numC > 2:
            p3 = fig.add_subplot(223)
            p3.errorbar(self.time, mag_list[2], yerr=error_list[2], marker='d', markersize =5, linestyle='')
            plt.xlabel('Time (JD-2400000)', fontsize=16)
            plt.ylabel('Target Mag (relative to C4)', fontsize=16)
            plt.gca().invert_yaxis()

        if self.numC > 3:
            p4 = fig.add_subplot(224)
            p3.errorbar(self.time, mag_list[3], yerr=error_list[3], marker='d', markersize =5, linestyle='')
            plt.xlabel('Time (JD-2400000)', fontsize=16)
            plt.ylabel('Target Mag (relative to C5)', fontsize=16)
            plt.gca().invert_yaxis()


    def CBA_print(self, mag_list, error_list):
        self.comp_star = int(input('Which comparison star are you using to submit the data? (i.e. 1,2,3,4,5,or 6) '))
        comp_star = self.comp_star
        self.check_star = int(input('Which comparison star are you using as the check star? (i.e. 1,2,3,4,5,or 6) '))
        check_star = self.check_star
        if self.comp_star == self.check_star:
            raise ValueError('Comp star and check star cannot be the same')
        
        comp_star = self.comp_star
        check_star = self.check_star
        if self.comp_star == 1:
            cmag = self.C2mag
            cnam = self.C2name
        elif self.comp_star == 2:
            cmag = self.C3mag
            cnam = self.C3name
        elif self.comp_star == 3:
            cmag = self.C4mag
            cnam = self.C4name
        elif self.comp_star == 4:
            cmag = self.C5mag
            cnam = self.C5name
        elif self.comp_star == 5:
            cmag = self.C6mag
            cnam = self.C6name
        elif self.comp_star == 6:
            cmag = self.C7mag
            cnam = self.C7name

        if self.check_star == 1:
            kmag = self.C2mag
            knam = self.C2name
        elif self.check_star == 2:
            kmag = self.C3mag
            knam = self.C3name
        elif self.check_star == 3:
            kmag = self.C4mag
            knam = self.C4mag
        elif self.check_star == 4:
            kmag = self.C5mag
            knam = self.C5name
        elif self.check_star == 5:
            kmag = self.C6mag
            knam = self.C6name
        elif self.check_star == 6:
            kmag = self.C7mag
            knam = self.c7name
            
        
        exposure_time=str(input("input the exposure length (s):"))
        observers= input("input the people observing: ")
        notes= input("input any notes for the night (na if none):")
        
        print()
        print('# Variable: ' + self.star_name)
        print('# Date: ' + self.day )
        print('# Comp star: ' + cnam + '_AAVSO_(V=',cmag,')')
        print('# Check star: ' + knam + '_AAVSO_(V=',kmag,')')
        print('# Exp time (s): ' + exposure_time + ' s')
        print('# Filter: ' + self.filter + '')
        print('# Observatory: Michigan State Campus Observatory')
        print('# East Lansing, MI, USA')
        print('# Observers: ' + observers + '')
        print('# Comments: ' + notes + '')
        print('# JD              Var_Mag     Var_eMag   Airmass')
        for ii in range(10):
            print('{0:13.5f}    {1:7.3f}    {2:7.3f}    {3:7.3f}'. format(self.time[ii], mag_list[comp_star-1][ii],
            error_list[comp_star-1][ii], self.airmass[ii]))
        print("...         Data Points:", self.time.size)
        contiue=input("Write to file? (yes or no) ").lower()
        if contiue== 'yes':
            f = open('CBA_'+self.star_name+self.day+'_measurements.txt','w')
           
    
            print('# Variable: ' + self.star_name,file=f)
            print('# Date: ' + self.day, file=f )
            print('# Comp star: ' + cnam + '_AAVSO_(V=',cmag,')',file=f)
            print('# Check star: ' + knam + '_AAVSO_(V=',kmag,')',file=f)
            print('# Exp time (s): ' + exposure_time + ' s',file=f)
            print('# Filter: ' + self.filter + '',file=f)
            print('# Observatory: Michigan State Campus Observatory',file=f)
            print('# East Lansing, MI, USA',file=f)
            print('# Observers: ' + observers + '',file=f)
            print('# Comments: ' + notes + '',file=f)
            print('# JD              Var_Mag     Var_eMag   Airmass',file=f)
            fmt='{0:13.5f}    {1:7.3f}    {2:7.3f}    {3:7.3f}'
            for ii in range(0, self.time.size):
                print(fmt.format(self.time[ii], mag_list[comp_star-1][ii], error_list[comp_star-1][ii], self.airmass[ii]),file=f)  
            f.close()

    def AAVSO_print(self, mag_list, error_list):

        comp_star = self.comp_star
        check_star = self.check_star
        if self.comp_star == 1:
            cmag = self.C2mag
            cnam = self.C2name
        elif self.comp_star == 2:
            cmag = self.C3mag
            cnam = self.C3name
        elif self.comp_star == 3:
            cmag = self.C4mag
            cnam = self.C4name
        elif self.comp_star == 4:
            cmag = self.C5mag
            cnam = self.C5name
        elif self.comp_star == 5:
            cmag = self.C6mag
            cnam = self.C6name
        elif self.comp_star == 6:
            cmag = self.C7mag
            cnam = self.C7name

        if self.check_star == 1:
            kmag = self.C2mag
            knam = self.C2name
        elif self.check_star == 2:
            kmag = self.C3mag
            knam = self.C3name
        elif self.check_star == 3:
            kmag = self.C4mag
            knam = self.C4mag
        elif self.check_star == 4:
            kmag = self.C5mag
            knam = self.C5name
        elif self.check_star == 5:
            kmag = self.C6mag
            knam = self.C6name
        elif self.check_star == 6:
            kmag = self.C7mag
            knam = self.c7name

       
        obscode = input('input your AAVSO observer code: \n')

        print('#TYPE=Extended') #Do not change
        print('#OBSCODE=',obscode) #Your unique ID. Change one time, then let it be.
        print('#SOFTWARE=Custom Python Script') #Do not change
        print('#DELIM=,') #Do not change
        print('#DATE=JD') #Do not change
        print('#OBSTYPE=CCD') #Do not change
        print('#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES')
        for ii in range(10):
            print(self.star_name,',','{0:13.5f},{1:6.3f},{2:5.3f}'.format(self.time[ii]+2400000, mag_list[comp_star-1][ii],
            error_list[comp_star-1][ii]),',',self.filter,',',self.trans,',',self.mtype,',',cnam,',',cmag,
            ',',knam,',',kmag,',','{0:5.3f}'.format(self.airmass[ii]),',',self.group,',',self.chart,',',self.notes,sep='')
        print("...         Data Points:", self.time.size)
        contiue=input("Write to text file? (yes or no) ").lower()
        if contiue== 'yes':
            fp = open('AAVSO_'+self.star_name+self.day+'_measurements.txt','w')
            print('#TYPE=Extended',file=fp) #Do not change
            print('#OBSCODE=',obscode,file=fp) #Your unique ID. Change one time, then let it be.
            print('#SOFTWARE=Custom Python Script',file=fp) #Do not change
            print('#DELIM=,',file=fp) #Do not change
            print('#DATE=JD',file=fp) #Do not change
            print('#OBSTYPE=CCD',file=fp) #Do not change
            print('#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES',file=fp)
            for ii in range(0, self.time.size):
                print(self.star_name,',','{0:13.5f},{1:6.3f},{2:5.3f}'.format(self.time[ii]+2400000, mag_list[comp_star-1][ii],
                    error_list[comp_star-1][ii]),',',self.filter,',',self.trans,',',self.mtype,',',cnam,',',cmag,
                    ',',knam,',',kmag,',','{0:5.3f}'.format(self.airmass[ii]),',',self.group,',',self.chart,',',self.notes,sep='',file=fp)
            fp.close
       
