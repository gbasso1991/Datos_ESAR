#%%!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vacuity.py
Created on Tue Aug 30 15:46:23 2022

@author: Giuliano

Para levantar los archivos .csv de temperatura gnersados por el sensor Rugged
"""
import pandas as pd
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
import os
import fnmatch
import tkinter as tk
from tkinter import filedialog
from datetime import datetime
#from possessor import medida_cruda


#%% Seleccion Archivos muestra
root = tk.Tk()
root.withdraw()
texto_encabezado = "Seleccionar archivos con las medidas de muestra:"
path_m=filedialog.askopenfilenames(title=texto_encabezado,filetypes=(("Archivos .txt","*.txt"),("Archivos .dat","*.dat"),("Todos los archivos","*.*")))
directorio = path_m[0].rsplit('/',maxsplit=1)[0]

fnames_m = []
for item in path_m:    
    fnames_m.append(item.split('/')[-1])

print('Directorio de trabajo: '+ directorio +'\n')
print('Archivos de muestra: ')
for item in fnames_m:
    print(' >',item)

#%% Params del nombre
'''
Parámetros de la medida a partir de nombre del archivo de muestra: 'xxxkHz_yyydA_zzzMss_*.txt
'''
frec_nombre=[] #Frec del nombre del archivo. Luego comparo con frec ajustada
Idc = []       #Internal direct current en el generador de RF
delta_t = []   #Base temporal 
fecha_m = []   #fecha de creacion archivo, i.e., de la medida 

for i in range(len(fnames_m)):
    frec_nombre.append(float(fnames_m[i].split('_')[0][:-3])*1000)
    Idc.append(float(fnames_m[i].split('_')[1][:-2])/10)
    delta_t.append(1e-6/float(fnames_m[i].split('_')[2][:-3]))
    fecha_m.append(datetime.fromtimestamp(os.path.getctime(path_m[i])).strftime('%H:%M:%S'))

#%% Creo los dataframes 
for fecha in fecha_m:
    print(fecha)
#df_0 = medida_cruda(path_m[0],delta_t[0])
#df_1 = medida_cruda(path_m[1],delta_t[1])



#%% Archivpos de Temperatura 
#directorio= 
data = pd.read_csv('agua/220830_templog_01.csv',sep=';',header=5,names=('Timestamp','Temperatura'),usecols=(0,1),decimal=',',dtype = {'Timestamp': 'str', 'Temperatura': 'float'},parse_dates=['Timestamp'],engine='python') 
data.plot(x='Timestamp',y='Temperatura',kind='scatter')
data2 = pd.read_csv('paramagneto/220830_templog_00.csv',sep=';',header=5,names=('Timestamp','Temperatura'),usecols=(0,1),decimal=',',parse_dates=['Timestamp'],engine='python') 


#%%
filenames = os.listdir(directorio)  
for templog in fnmatch.filter(filenames,'*templog*'):
    print(templog)
#%%Timedelta
data['timedelta']=(data['Timestamp']-data['Timestamp'][0])
data2['timedelta']=data2['Timestamp']-data2['Timestamp'][0]

