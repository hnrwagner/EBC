# --------------------------------------------------------------
# Python script for the analysis of FEA results
# 
#
#
# Dr.-Ing. Ronald Wagner
# Siemens Mobility GmbH
# 25.07.2019
# --------------------------------------------------------------

# Import the required libraries

#---------------------------------------------------------------

import pandas as pd
import xlsxwriter
import numpy as np
import matplotlib.pyplot as plt

 # define functions

def createDataColumn(row,col,data,data_modification):
    """
    This function creates excel columns from data
    Input:
        row - number of the row in excel
        col - number of the column in excel
        data - input data for the excel sheet
        data_modification - modifyer for the input data
    """
    for i in (data):
        worksheet.write_number  (row,col,i/data_modification)
        row += 1
        
def createInput(input_string):
    """
    This function read the input data of the FEA from a text file
    Input:
        input_string - name of the input data
    Output: 
        values - the input data 
    """
    data = np.loadtxt(input_string)
    values = data[:]
    return values
        

#------------------------------------------------------------------------------   
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Name of the Finite Element Analysis Modell    
#
#------------------------------------------------------------------------------   
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------ 


myName = ['S1_Z_715_angle_0']

# define axis title 

x_axis = ['Displacement u [mm]', 'Relative Stiffness S [kN/mm]', 'Relative Strain Energy E [kNmm/kN]']
y_axis = ['Axial Force F [kN]','Axial Force F [kN]','Axial Force F [kN]', 'Buckling Load N [kN]']

max_load = ['Numerical Collapse Load CL [kN]','Local Buckling Load LL [kN]']

# method which is evaluated

method = ['LRSM']

# axis of perturbation approach

pert_axis = ['LRSM radius-to-shell radius ratio, Rs/R']

# limits for the perturbation iteration

Pert_START =1
Pert_END = 11

for jc in range(0,1,1):
    filename = 'Design_'+str(myName[jc])+'_008.xlsx'
    workbook = xlsxwriter.Workbook(filename)
    workbook = xlsxwriter.Workbook(filename, {'nan_inf_to_errors': True})
    worksheet = []
    local_load = []
    local_load2 = []
    strain = []
    
    for ic in range(Pert_START,Pert_END,1):
        values_Int = []
        
        # define strings for excel sheet tabs
        
        myString = str(myName[jc])+'_EBC_Loop-'+str(ic)
        myString2 = str(myName[jc])+'-'+str(ic)
        worksheet =  workbook.add_worksheet(myString2)
        
        # read data
        
        values_N = createInput('Rf_'+str(myString)+'.txt')
        values_w = createInput('u_'+str(myString)+'.txt')
        values_N3y = createInput('v_'+str(myName[jc])+'.txt')
        
        # calculate derivative of axial force based on axial displacement
        # Relative Stiffness
        
        values_Diff = np.diff(values_N)/np.diff(values_w)
        values_Diff_2 = np.diff(values_Diff)/np.diff(values_w[0:len(values_w)-1])
#        aa = np.array(values_Diff)
#        min_max = np.where(aa == 0)[0]
        
        # calculate Integral of axial force based on axial displacement
        # Strain Energy
        
        for i in range(1,len(values_N),1):
            values_Int.append(np.trapz(values_N[0:i],values_w[0:i]))
            #print(values_N[i],i)
        
        # calculate derivative of Strain Energy based on axial force
        # Relative Strain Energy
        
        values_Int = np.asarray(values_Int)
        values_Diff2 = np.diff(values_Int)/np.diff(values_N[0:len(values_N)-1])
        
        # Numerical criterion for local buckling - 1 - Stiffness Criterion
        # Change of axial stiffness by more than 5 % equals local buckling (values based on experience by Dr. Wagner)
        
        for i in range(0,len(values_Diff)-1,1):
            a = values_Diff[i]/values_Diff[0]
            #print(a)
            if a < 0.95:
                local_load.append(values_N[i-1]/1000)
                break
            
        # Numerical criterion for local buckling - 2 - Energy Criterion
#        TH = values_N[np.argmax(values_N)]*0.2
#        
#        for i in range(0,len(values_N),1):
#            if values_N[i] >= TH:
#                index = i
#                break
#        
#        m,c = np.polyfit(values_w[0:index],values_N[0:index]/1000,1)
#        
#        y = m * values_Diff2 + c
        
        values_Diff3 = np.diff(values_N[0:len(values_N)-3])/np.diff(values_Diff2[0:len(values_N)-3])
        
        for i in range(0,len(values_Diff3),1):
            if values_Diff3[i] <= 0:
                index = i
                break
            else:
                index = np.argmax(values_N)
                
            
        local_load2.append(values_N[index]/1000)    



        
        bold = workbook.add_format({'bold': 1})        
        
        # Adjust the column width.
        
        worksheet.set_column(0, 0, 15)
        worksheet.set_column(1, 1, 20)
        worksheet.set_column(2, 2, 31)
        worksheet.set_column(3, 3, 25)
        worksheet.set_column(5, 5, 1)
        worksheet.set_column(4, 4, 27)
        worksheet.set_column(6, 6, 30)
        
        # define title for columns in excel
        
        worksheet.write('A1', y_axis[0], bold)
        worksheet.write('D1', 'Relative Stiffness [N/mm]', bold)
        worksheet.write('C1', 'Relative Strain Energy S [Nmm/kN]', bold)
        worksheet.write('B1', x_axis[0], bold)
        worksheet.write('G1', myString2, bold)
        worksheet.write('E1', 'Perfect Relative Strain Energy', bold)
        
        # create Datas in Excel (columns)
        
        createDataColumn(3,0,abs(values_N),1000)
        createDataColumn(3,1,abs(values_w),1)
        createDataColumn(3,2,values_Diff,1)
        createDataColumn(3,3,values_Diff_2,1)
        #createDataColumn(3,4,y,1)
        
        # define max. Load

        worksheet.write('M4', max_load[0], bold)
        worksheet.write('M5', np.max(values_N), bold)
        #worksheet.write('M7', max_load[1], bold)
        #worksheet.write('M8', local_load(ic), bold)
        ##--------------------------------------------------------------------------------
        
        # Insert First Diagramm for Reaction Load vs. Displacement
        
        ##--------------------------------------------------------------------------------
        
        chart = workbook.add_chart({'type': 'scatter'})
        
        # Limits for xy chart
        
        max_row=len(values_N)
        for i in range(0,1,1):
            col = i 
            chart.add_series({
                #'name':       [myString2, 0, 6],
                'categories': [myString2, 3, col+1, max_row+3, col+1],
                'values':     [myString2, 3, col, max_row+3, col],
                'marker':   {'type': 'none'},
                'line':     {'color': 'red', 'width': 0.5},})
    
        # Configure the chart axes.
        chart.set_y_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'name_font': {'size': 14, 'bold': True},'name': str(x_axis[0]),'min': 0.0,'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_y_axis({'name_font': {'size': 14, 'bold': True},'name': str(y_axis[0]),'min': 0.0,'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_legend({'none': True})

        
        # insert diagramm
        
        worksheet.insert_chart('G9', chart)
        
        
        ##--------------------------------------------------------------------------------
        
        # Insert Second Diagramm for Stiffness vs. Reaction Load
        
        ##--------------------------------------------------------------------------------
        
        
        chart = workbook.add_chart({'type': 'scatter'})
        
        # Limits for xy chart
        
        max_row=len(values_N)
        for i in range(0,1,1):
            col = i 
            chart.add_series({
                #'name':       [myString2, 0, 6],
                'categories': [myString2, 3, col+3, max_row+3, col+3],
                'values':     [myString2, 3, col, max_row+3, col],
                'marker':   {'type': 'none'},
                'line':     {'color': 'red', 'width': 0.5},})
    
        # Configure the chart axes.
        chart.set_y_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'name_font': {'size': 14, 'bold': True},'name': str(x_axis[1]),'min': 0.0,'max': values_Diff[0]*1.1,'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_y_axis({'name_font': {'size': 14, 'bold': True},'name': str(y_axis[1]),'min': 0.0,'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_legend({'none': True})
        
        
        # insert diagramm
        
        #worksheet.insert_chart('M24', chart)
        
        
        ##--------------------------------------------------------------------------------
        
        # Insert Third Diagramm for Reaction Load vs. Relative Strain Energy
        
        ##--------------------------------------------------------------------------------
        
        
        chart = workbook.add_chart({'type': 'scatter'})
        
        # Limits for xy chart
        
        max_row=len(values_N)
        for i in range(0,1,1):
            col = i 
            chart.add_series({
                #'name':       [myString2, 0, 6],
                'categories': [myString2, 3, col+1, max_row+3, col+1],
                'values':     [myString2, 3, col+2, max_row+3, col+2],
                'marker':   {'type': 'none'},
                'line':     {'color': 'red', 'width': 0.5},})
    
#            chart.add_series({
#                #'name':       [myString2, 0, 6],
#                'categories': [myString2, 3, col+2, max_row+3, col+2],
#                'values':     [myString2, 3, col+4, max_row+3, col+4],
#                'marker':   {'type': 'none'},
#                'line':     {'color': 'black', 'width': 0.5},})
    
    
        # Configure the chart axes.
        chart.set_y_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'num_font':  {'name': 'Times New Roman'}, })
        chart.set_x_axis({'name_font': {'size': 14, 'bold': True},'name': str(x_axis[2]),'min': 0.0,'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_y_axis({'name_font': {'size': 14, 'bold': True},'name': str(y_axis[2]),'min': np.min(values_Diff),'max': np.max(values_Diff),'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
        chart.set_legend({'none': True})
        
        
        # insert diagramm
        
        #worksheet.insert_chart('G24', chart)
        
        
        ##--------------------------------------------------------------------------------
        
        ##--------------------------------------------------------------------------------
        
        
        
        # define load vs. perturbation load diagramm
        
        if ic == Pert_END-1:
            myString3 = 'Global Buckling Load'
            myString4 = 'Local Buckling Load - for elastic-plastic buckling'
            myString6 = 'Local Buckling Load - for elastic buckling'
            myString5 = 'Lower-Bound Diagram'
            worksheet =  workbook.add_worksheet(myString5)
            worksheet.write('A1',myString , bold)  
            
            # define input for perturbation diagramm
            
            values_N1 = createInput('N_'+str(myName[jc])+'.txt')
            values_N2 = createInput('u_'+str(myName[jc])+'.txt')
            values_N3 = createInput('v_'+str(myName[jc])+'.txt')

            # define title for columns in excel
            
            worksheet.write('A1', str(myString3), bold)  
            worksheet.write('B1', str(x_axis[0]), bold)  
            worksheet.write('C1', str(pert_axis[0]), bold) 
            worksheet.write('D1', str(myString4), bold) 
            worksheet.write('E1', str(myString6), bold) 
            
            # create Datas in Excel (columns)
            
            createDataColumn(3,0,abs(values_N1),1000)
            createDataColumn(3,1,abs(values_N2),1)
            createDataColumn(3,2,abs(values_N3),1)
            createDataColumn(3,3,abs(np.asarray(local_load)),1)
            createDataColumn(3,4,abs(np.asarray(local_load2)),1)
            #createDataColumn(3,4,abs(np.asarray(local_pert)),1)
            
            # create xy - chart from Data
                
            chart = workbook.add_chart({'type': 'scatter'})
    
            # Configure the series of the chart from the dataframe data.
            max_row = len(values_N1)
            for i in range(0,1,1):
                col = i 
                chart.add_series({
                    'name':       [myString5, 0, 0],
                    'categories': [myString5, 3, col+2, max_row+3, col+2],
                    'values':     [myString5, 3, col, max_row+3, col],
                    'marker':     {'type': 'plus', 'size': 7, 'border': {'color': 'red'}},})
    
                chart.add_series({
                    'name':       [myString5, 0, 3],
                    'categories': [myString5, 3, col+2, max_row+3, col+2],
                    'values':     [myString5, 3, col+3, max_row+3, col+3],
                    'marker':     {'type': 'diamond', 'size': 7,'fill':   {'none': True}, 'border': {'color': 'black'}},})
    
            chart.set_plotarea({'fill':   {'color': 'white'}})
    
            # Configure the chart axes.
            chart.set_y_axis({'num_font':  {'name': 'Times New Roman'}})
            chart.set_x_axis({'num_font':  {'name': 'Times New Roman'}})
            chart.set_x_axis({'name_font': {'size': 14, 'bold': True},'name': str(pert_axis[0]), 'min': 0.0, 'max': max(values_N3),'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
            chart.set_y_axis({'name_font': {'size': 14, 'bold': True},'name': str(y_axis[3]), 'min': 0.0, 'max': max(values_N1/1000),'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
            #chart.set_legend({'none': True})
            chart.set_legend({'layout': {'x':      1.0,'y':      0.0,'width':  1.0,'height': 0.05,}})
            chart.set_plotarea({'layout': {'x':      0.20,'y':      0.20,'width':  0.8,'height': 0.6,}})
    
            # Insert the chart into the worksheet.
            worksheet.insert_chart('W4', chart)
            
            # create xy - chart from Data
#                
            chart = workbook.add_chart({'type': 'scatter'})
    
            # Configure the series of the chart from the dataframe data.
            max_row = len(values_N1)
            for i in range(0,1,1):
                col = i 
                chart.add_series({
                    'name':       [myString5, 0, 0],
                    'categories': [myString5, 3, col+2, max_row+3, col+2],
                    'values':     [myString5, 3, col, max_row+3, col],
                    'marker':     {'type': 'plus', 'size': 7, 'border': {'color': 'red'}},})
    
                chart.add_series({
                    'name':       [myString5, 0, 4],
                    'categories': [myString5, 3, col+2, max_row+3, col+2],
                    'values':     [myString5, 3, col+4, max_row+3, col+4],
                    'marker':     {'type': 'diamond', 'size': 7,'fill':   {'none': True}, 'border': {'color': 'black'}},})
    
            chart.set_plotarea({'fill':   {'color': 'white'}})
    
            # Configure the chart axes.
            chart.set_y_axis({'num_font':  {'name': 'Times New Roman'}})
            chart.set_x_axis({'num_font':  {'name': 'Times New Roman'}})
            chart.set_x_axis({'name_font': {'size': 14, 'bold': True},'name': str(pert_axis[0]), 'min': 0.0, 'max': max(values_N3),'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
            chart.set_y_axis({'name_font': {'size': 14, 'bold': True},'name': str(y_axis[3]), 'min': 0.0, 'max': max(values_N1/1000),'major_gridlines': {'visible': True,'line':{'color': '#CAE1FF'}}})
            #chart.set_legend({'none': True})
            chart.set_legend({'layout': {'x':      1.0,'y':      0.0,'width':  1.0,'height': 0.05,}})
            chart.set_plotarea({'layout': {'x':      0.20,'y':      0.20,'width':  0.8,'height': 0.6,}})
    
            # Insert the chart into the worksheet.
            worksheet.insert_chart('G4', chart)
           
    workbook.close()
        
            
        

