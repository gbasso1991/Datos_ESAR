#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
god_of_thunder.py

Created on Tue Aug 30 14:14:59 2022

@author: giuliano

Para comparar fondos pre y post medidas en temperatura
"""
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import os
from datetime import datetime
from possessor import medida_cruda
#%% Seleccion
root = tk.Tk()
root.withdraw()
texto_encabezado = "Seleccionar archivos con las medidas de fondo a comparar:"
path_m=filedialog.askopenfilenames(title=texto_encabezado,filetypes=(("Archivos .txt","*.txt"),("Archivos .dat","*.dat"),("Todos los archivos","*.*")))
directorio = path_m[0].rsplit('/',maxsplit=1)[0]

fnames_m = []
for item in path_m:    
    fnames_m.append(item.split('/')[-1])

print('Directorio de trabajo: '+ directorio +'\n')
print('Archivos de fondo a comparar: ')
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
    fecha_m.append(datetime.fromtimestamp(os.path.getmtime(path_m[i])).strftime('%d-%m-%Y %H:%M'))

#%% Creo los dataframes 

df_0 = medida_cruda(path_m[0],delta_t[0])
df_1 = medida_cruda(path_m[1],delta_t[1])

resta = df_1['v']-df_0['v']


fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(20,8),constrained_layout=True)
#plt.figure(figsize=(10,8),constrained_layout=True)
ax1.plot(df_0['t'],df_0['v'],lw=0.9, label=fnames_m[0],alpha=1,zorder=0)
ax1.plot(df_1['t'],df_1['v'],lw=0.9, label=fnames_m[1],alpha=0.7,zorder=1)

ax1.legend(ncol=2)
ax1.grid()
ax1.set_ylabel('$\epsilon (V)$')
ax1.set_ylim(-0.03,0.03)


ax2.plot(df_0['t'],resta,lw=0.9,label='Resta',c='tab:green')
plt.xlim(0,max(df_0['t']))
ax2.grid()
plt.legend()
plt.xlabel('tiempo (s)')
ax2.set_ylabel('$\epsilon (V)$')
ax2.set_ylim(-0.03,0.03)
plt.suptitle('Resta de fondos pre y post medida',fontsize=20)
plt.savefig('Resta_fondos.png',dpi=300)
