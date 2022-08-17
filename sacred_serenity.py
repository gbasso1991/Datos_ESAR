#%% sacred_serenity.py
# https://www.youtube.com/watch?v=-as4LausEok&ab_channel=deathmetalchuck 
# Giuliano Basso 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from uncertainties import ufloat , unumpy
import pandas as  pd 
import seaborn as sns

#%%
#Tomo variables de possessor.py relativos a calibracion: 
                # pendiente_cal 
                # pendiente_cal_filtrada 
                # Campo_maximo_

pendiente_cal = Pendiente_cal
pendiente_cal_filtrada= Pendiente_cal_filtrada
Campo_kAm = Campo_maximo_c_kAm 
frecuencia_kHz = Frecuencia_c_kHz

#%% Resultados copiados para no correr innecesariamente possessor.py
pend_cal_aux = [2.1938891523014237e-14,
 2.66337420317448e-14,
 2.244892591284205e-14,
 2.5093936186902377e-14,
 2.2889081540309153e-14,
 2.035213826593563e-14,
 2.4802887528205646e-14,
 2.442687101041686e-14,
 2.4637249282335315e-14,
 2.3833582612029775e-14,
 2.4925024562133184e-14,
 2.2742834965158362e-14,
 2.0156486239419658e-14,
 2.130720596527412e-14,
 2.039952033243038e-14,
 2.083479358156623e-14,
 2.08259781660887e-14,
 2.1290206786958435e-14,
 1.9015664948463763e-14,
 2.383347577393824e-14,
 2.3660689673351158e-14,
 2.1468029246588e-14,
 1.9703556820356692e-14,
 1.7886554368547193e-14,
 2.0698230409970277e-14,
 2.0027969293701767e-14,
 1.8595471582505443e-14,
 1.9821678409493157e-14,
 1.9670103979176385e-14,
 1.8789227795546812e-14,
 2.2154513008749434e-14,
 2.3805172181294694e-14,
 2.1735039862426956e-14,
 2.0896137691287775e-14,
 1.6228986353358634e-14,
 2.0758206032199146e-14,
 1.5046265480912106e-14,
 2.222718611896522e-14,
 2.0695427889536545e-14,
 1.950270019671407e-14,
 2.4435526954187375e-14,
 2.139256092426729e-14,
 1.9703259553397456e-14,
 1.9721501577121833e-14,
 1.782460629237128e-14,
 1.9955444450676015e-14,
 2.2980683601724387e-14,
 1.9256291032974623e-14,
 2.73409814653112e-14,
 2.512655351008897e-14]
# ,
#  2.4594172350456962e-14,
#  2.4717353653770868e-14,
#  2.165195131639293e-14,
#  2.534677094200584e-14,
#  1.8756976823989546e-14,
#  1.6174394441555504e-14,
#  1.6868972930694622e-14,
#  2.0641287619492046e-14]
pendiente_cal_filtrada_aux = [2.230371025045133e-14,
 2.6633362165512437e-14,
 2.2711331346487755e-14,
 2.5093506345441616e-14,
 2.289233331881107e-14,
 2.001183327954973e-14,
 2.4791732688168448e-14,
 2.4420751504844277e-14,
 2.487755258437644e-14,
 2.3804643597248577e-14,
 2.4899067862560176e-14,
 2.273953752292939e-14,
 1.9874693801840522e-14,
 2.197583638682677e-14,
 2.0342895959248715e-14,
 2.0852044273321304e-14,
 2.111605023189915e-14,
 2.1245001713454813e-14,
 1.9207668542355723e-14,
 2.402516853278747e-14,
 2.353085378263747e-14,
 2.147154402219688e-14,
 1.9620924284512427e-14,
 1.8388243717878196e-14,
 2.0814762738618344e-14,
 2.0686412754746578e-14,
 1.8537103577675097e-14,
 1.986591060988049e-14,
 1.9669432898263832e-14,
 1.8546974678733467e-14,
 2.1263073943794346e-14,
 2.3709676381598414e-14,
 2.1760574215818546e-14,
 2.0919707482915188e-14,
 1.63431746395627e-14,
 2.1577421560274614e-14,
 1.5399828005459692e-14,
 2.2148546567211178e-14,
 1.987276469526923e-14,
 1.928004789861137e-14,
 2.4724941952520637e-14,
 2.1554033901730634e-14,
 2.013158098784757e-14,
 1.9956702372169692e-14,
 1.7219168062353743e-14,
 2.0025192524599326e-14,
 2.235354908666882e-14,
 1.9386978853227608e-14,
 2.7103684452884877e-14,
 2.5127883490905795e-14]
#  ,
#  2.4251845087922003e-14,
#  2.4546916039259263e-14,
#  2.1010191526637006e-14,
#  2.5270784344697774e-14,
#  1.879932349938216e-14,
#  1.6463238085697296e-14,
#  1.6894669378370844e-14,
#  2.0369234645462595e-14]
Campo_kAm_aux=[10.551177899999999,
 10.551177899999999,
 20.884583699999997,
 20.884583699999997,
 31.217989499999995,
 31.217989499999995,
 41.551395299999996,
 41.551395299999996,
 51.8848011,
 51.8848011,
 10.551177899999999,
 10.551177899999999,
 20.884583699999997,
 20.884583699999997,
 31.217989499999995,
 31.217989499999995,
 41.551395299999996,
 41.551395299999996,
 51.8848011,
 51.8848011,
 10.551177899999999,
 10.551177899999999,
 20.884583699999997,
 20.884583699999997,
 31.217989499999995,
 31.217989499999995,
 41.551395299999996,
 41.551395299999996,
 51.8848011,
 51.8848011,
 10.551177899999999,
 10.551177899999999,
 20.884583699999997,
 20.884583699999997,
 31.217989499999995,
 31.217989499999995,
 41.551395299999996,
 41.551395299999996,
 51.8848011,
 51.8848011,
 10.551177899999999,
 10.551177899999999,
 20.884583699999997,
 20.884583699999997,
 31.217989499999995,
 31.217989499999995,
 41.551395299999996,
 41.551395299999996,
 51.8848011,
 51.8848011]
#,
#  3.6622407,
#  3.6622407,
#  10.551177899999999,
#  10.551177899999999,
#  31.217989499999995,
#  31.217989499999995,
#  51.8848011,
#  51.8848011]

frecuencia_kHz= [97.81548032872013,
 97.86157830340441,
 97.81571976773141,
 97.79541459894554,
 97.86107860796157,
 97.83964393876292,
 97.88995631658761,
 97.85816284144794,
 97.92872850627569,
 97.93724659313501,
 134.43895143307802,
 134.49424250879383,
 134.52302937349327,
 134.47875279790742,
 134.56904541381456,
 134.5570263340215,
 134.70855133099474,
 134.61145625750152,
 134.83682079113555,
 134.82120400852617,
 208.37883321279506,
 208.4016690936,
 208.27016841877617,
 208.26864246670078,
 208.28309606674873,
 208.2381935312378,
 208.35553122117014,
 208.30444901444224,
 208.39728839105535,
 208.42873742262074,
 234.2910606600161,
 234.3125114726098,
 234.11897229981182,
 234.1098966134046,
 234.10561243219527,
 234.00686354735663,
 234.0003112973212,
 233.9629415828658,
 234.0648754652537,
 234.12562713468336,
 262.7039636189564,
 262.66681030159407,
 262.07152634641454,
 262.0267261238117,
 262.00382234736213,
 261.96478094256594,
 262.1989541727353,
 262.1722292987853,
 262.4509069701014,
 262.4703176515558]
#  ,
#  300.1981712521084,
#  300.18782260803874,
#  299.9391866193163,
#  299.88368556178483,
#  299.042497258614,
#  299.0677094242394,
#  299.8151925759934,
#  299.8096066689678]

#%% Realizo estadistica sobre freq, campo, y pendiente

#Frec: promedio en 10 valores 
frecuencia=[]
frecuencia_aux = []
for i in range(0,len(frecuencia_kHz),10):
    frecuencia.append(ufloat(np.mean(frecuencia_kHz[i:i+9]),np.std(frecuencia_kHz[i:i+9])))
#Campo Max
campo_aux = []
for j in range(0,len(Campo_kAm_aux),2):
    campo_aux.append(round(np.mean(Campo_kAm_aux[j:j+1]),2))
campo = np.unique(campo_aux)    
#Pendiente
m=[]
m_incert=[]
for k in range(0,len(pend_cal_aux),2):
    m.append(np.mean(pend_cal_aux[k:k+9]))
    m_incert.append(ufloat(np.mean(pend_cal_aux[k:k+9]), np.std(pend_cal_aux[k:k+9])))

m_matriz = np.asarray(m)
m_matriz= m_matriz.reshape(5,5)

df_pendiente = pd.DataFrame(m_matriz,columns=campo,index=frecuencia )
#df_pendiente= df_pendiente.reindex(index=df_pendiente.index[::-1])
texto_m = (np.asarray(["{:.2e}".format(value) for value in m_incert]).reshape(5, 5))
  
#Pendiente filtrada 
m_filtrada = []
m_filtrada_incert=[]
for l in range(0,len(pendiente_cal_filtrada_aux),2):
    m_filtrada.append(np.mean(pend_cal_aux[k:k+9]))
    m_filtrada_incert.append(ufloat(np.mean(pendiente_cal_filtrada_aux[l:l+9]),np.std(pendiente_cal_filtrada_aux[l:l+9])))
m_filtr_new = np.asarray(m_filtrada)
m_filtr_new= m_filtr_new.reshape(5,5)    
df_pendiente_filtrada = pd.DataFrame(m_filtr_new,columns=campo,index=frecuencia )

texto_m_filtrada= (np.asarray(["{:.2e}".format(value) for value in m_filtrada_incert]).reshape(5, 5))



#%%  

#%% Dataframes 

#df_pendiente_filtrada= df_pendiente_filtrada.reindex(index=df_pendiente_filtrada.index[::-1])
#%%

#%%
sns_plot = sns.heatmap(df_pendiente,annot=texto_m, fmt="",linewidths=0.25)
sns.set(rc={"figure.figsize":(10,8)})
fig = sns_plot.get_figure()

plt.title('Pendiente de calibracion',fontsize=15)
plt.xlabel('Campo (kA/m)')
plt.ylabel('Frecuencia (kHz)')
plt.tight_layout()
#fig.savefig("heatmap_pendiente_calibracion.png",dpi=100)
#%%

sns_plot = sns.heatmap(df_pendiente_filtrada,annot=texto_m_filtrada, fmt="",linewidths=0.25)
sns.set(rc={"figure.figsize":(10,8)})
fig = sns_plot.get_figure()


plt.title('Pendiente de calibracion (filtrada)',fontsize=15)
plt.xlabel('Campo (kA/m)')
plt.ylabel('Frecuencia (kHz)')
plt.tight_layout()
fig.savefig("heatmap_pendiente_calibracion_filtrada.png",dpi=100)
#%%
sns.set()
flights = df_pendiente
flights = flights.pivot("campo", "frecuencia", "pendiente")
ax = sns.heatmap(flights)
plt.title("Heatmap Flight Data")
plt.show()
#%%
fig, ax = plt.subplots()
im = ax.imshow(df_pendiente)
cbar = ax.figure.colorbar(im, ax=ax)
ax.set_xticks(np.arange(len(campo)), label=campo)
ax.set_yticks(np.arange(len(frecuencia)), label=frecuencia)

ax.set_title('Pendiente')


#%%
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
#%%
# Import libraries

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
 
# Define Data

x = frecuencia_kHz
y = Campo_kAm
z = pendiente_cal_filtrada

# Create Figure

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
 
# Create Plot

ax.scatter3D(x, y, z)
 
# Show plot

plt.show()
#%%
ax = plt.figure(figsize=(10,8)).add_subplot(projection='3d')

# setting up a parametric curve
t = np.arange(0, 2*np.pi+.1, 0.01)
x = frecuencia_kHz
y = Campo_kAm
z = pendiente_cal_filtrada
estep = 15

ax.errorbar(x, y, z, 0.2)


ax.set_xlabel("X label")
ax.set_ylabel("Y label")
ax.set_zlabel("Z label")

plt.show()
#%
# 
# %

