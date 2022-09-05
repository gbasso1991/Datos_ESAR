#%%!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vacuity.py
Created on Tue Aug 30 2022

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
from scipy.interpolate import BSpline
#% Seleccion Archivos muestra
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
#%%
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

delta_t_m = [] #seg
for elem in fecha_m:
    delta_t_m.append(int((elem - fecha_m[0]).total_seconds()))

# #%% En windows, me conviene usar mtime, despues de todo, busco un timedelta
# datetime.fromtimestamp(os.stat(path_m[0]).st_ctime).strftime('%F %H:%M:%S')
# #%% 
# datetime.fromtimestamp(os.stat(path_m[0]).st_mtime).strftime('%F %H:%M:%S')

#% Ahora levanto el log de temperaturas en .csv

directorio= 'paramagneto/220830_templog_00.csv' #automatizar esto
#directorio= 'agua/220830_templog_01.csv' #automatizar esto

data = pd.read_csv(directorio,sep=';',header=5,
                    names=('Timestamp','Temperatura'),usecols=(0,1),
                    decimal=',',engine='python') 
temperatura = pd.Series(data['Temperatura']).to_numpy(dtype=float)

timestamp=[]
for time in pd.Series(data['Timestamp']):
    timestamp.append(time[11:19])
timestamp=np.array(timestamp,dtype=str)

#% Datetime de las medidas de muestra en funcion del horario del 1er dato

date_primer_dato = datetime(year=2022,month=8,day=30,
                            hour=10,minute=39,second=59) #queda automatizar esto 
#para agua, cambiar por 11:07:11 


time_m = [] 
for elem in delta_t_m:
    time_m.append((date_primer_dato + timedelta(0,elem)).strftime('%H:%M:%S'))
time_m = np.array(time_m) #H:M:S del registro de cada archivo muestra

#obtengo los indices de estos horarios en el timestamp

temp_m = []
indx_temp = []
for t in time_m:
    indx_temp.append(np.flatnonzero(timestamp==t)[0])
    temp_m.append(temperatura[np.flatnonzero(timestamp==t)[0]])
temp_m=np.array(temp_m)
#% Printeo lo obtenido
print('Archivos de muestra: ')
for i, item in enumerate(fnames_m):
    print(item[-8:-4],' ',str(temp_m[i]) + ' ºC')






#%% dif temporal en s entre el comienzo del registro y la primer medida de muestra
#delta_0 = (datetime.strptime(time_m[0],'%H:%M:%S') - datetime.strptime(timestamp[0],'%H:%M:%S')).total_seconds()

#%%
# timestamp = pd.to_datetime(pd.Series(data['Timestamp']).to_numpy())
# aux=[]
# for element in timestamp:
#     aux.append(element.strftime('%H:%M:%S'))

#for e in timestamp:
 #   e = datetime.strptime(e[11:19],'%H:%M:%S')
    #print((e,)
#plt.xticks(rotation=45, ha='right')
#%% Ploteo los ciclos y asocio la temperatura
#traigo los datos desde possessor

Ciclos_eje_H=Ciclos_eje_H
Ciclos_eje_M=Ciclos_eje_M
Ciclos_eje_H_cal=Ciclos_eje_H_cal
Ciclos_eje_M_cal_ua=Ciclos_eje_M_cal_ua
SAR=SAR

fig = plt.figure(figsize=(10,8),constrained_layout=True)
ax = fig.add_subplot(1,1,1)
axin = ax.inset_axes([0.60,0.08, 0.35,0.38])
axin.set_title('Calibración',loc='center')
#axin.yaxis.tick_right()
plt.setp(axin.get_yticklabels(),visible=True)
plt.setp(axin.get_xticklabels(),visible=True)
axin.yaxis.tick_right()
axin.grid()
axin.axhline(0,0,1,lw=0.9,c='k')
axin.axvline(0,0,1,lw=0.9,c='k')

for i in range(len(fnames_m)):      
    plt.plot(Ciclos_eje_H[i],Ciclos_eje_M[i],label=f'{fnames_m[i][-8:-4]}   {temp_m[i]}ºC')
    axin.plot(Ciclos_eje_H_cal[i], Ciclos_eje_M_cal_ua[i])
    axin.set_ylabel('M $(V\cdot s)$')
    axin.set_xlabel('H $(A/m)$')

plt.legend(loc='best',fancybox=True)
plt.grid()
plt.xlabel('Campo (A/m)',fontsize=15)
plt.ylabel('Magnetización (A/m)',fontsize=15)
plt.suptitle('Ciclos de histéresis en descongelamiento',fontsize=30)
#plt.savefig('Ciclos_histeresis_descong_bis.png',dpi=300,facecolor='w')

# %%
