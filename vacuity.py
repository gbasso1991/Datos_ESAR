#%%!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vacuity.py
Created on Tue Aug 30 15:46:23 2022

@author: Giuliano

Para levantar los archivos .csv de temperatura gnersados por el sensor Rugged

ATENCION: Getting file creation dates, on the other hand, is fiddly 
          and platform-dependent: 
    On Windows, a file's ctime  stores its creation date. 
    You can access this in Python through os.path.getctime() or 
    the .st_ctime attribute of the result of a call to os.stat(). 
    This won't work on Unix, where the ctime is the last time that 
    the file's attributes or content were changed.

    On Linux, this is currently impossible, at least without writing 
    a C extension for Python. Although some file systems commonly used 
    with Linux do store creation dates, the Linux kernel offers no way 
    of accessing them.

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
import os
import fnmatch
import tkinter as tk
from tkinter import filedialog
from datetime import datetime

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

'''
Parámetros de la medida a partir de nombre del archivo de muestra: 'xxxkHz_yyydA_zzzMss_*.txt
'''
frec_nombre=[]      #Frec del nombre del archivo. Luego comparo con frec ajustada
Idc = []            #Internal direct current en el generador de RF
samples_per_s = []  #Base temporal 
fecha_m = []        #fecha de creacion archivo, i.e., de la medida 

for i in range(len(fnames_m)):
    frec_nombre.append(float(fnames_m[i].split('_')[0][:-3])*1000)
    Idc.append(float(fnames_m[i].split('_')[1][:-2])/10)
    samples_per_s.append(1e-6/float(fnames_m[i].split('_')[2][:-3]))
    fecha_m.append(datetime.fromtimestamp(os.path.getmtime(path_m[i])))

delta_t_m = []
for elem in fecha_m:
    delta_t_m.append((elem - fecha_m[0]).total_seconds())

# #%% En windows, me conviene usar mtime, despues de todo, busco un timedelta
# datetime.fromtimestamp(os.stat(path_m[0]).st_ctime).strftime('%F %H:%M:%S')
# #%% 
# datetime.fromtimestamp(os.stat(path_m[0]).st_mtime).strftime('%F %H:%M:%S')

#%% Ahora levanto el log de temperaturas en .csv
#directorio= 
data = pd.read_csv('agua/220830_templog_01.csv',sep=';',header=5,names=('Timestamp','Temperatura'),usecols=(0,1),decimal=',',dtype = {'Timestamp': 'str', 'Temperatura': 'float'},parse_dates=['Timestamp'],engine='python') 
temperatura = pd.Series(data['Temperatura'])
temperatura = temperatura.to_numpy()

#= pd.read_csv('paramagneto/220830_templog_00.csv',sep=';',header=5,names=('Timestamp','Temperatura'),usecols=(0,1),decimal=',',parse_dates=['Timestamp'],engine='python') 

timestamp = pd.Series(data['Timestamp'])

timestamp = timestamp.to_numpy(dtype=datetime)

for e in timestamp:
    e = datetime.strptime(e[11:19],'%H:%M:%S')
    #print((e,)
#%%




#plt.xticks(rotation=45, ha='right')