import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Data.csv')

time = df['Time'].values.tolist()
cIV = df['Cocaine Intravenous (ng/mL)'].values.tolist()
bIV = df['Benzoylecgonine Intravenous (ng/mL)'].values.tolist()
cIN = df['Cocaine Intranasal (ng/mL)'].values.tolist()
bIN = df['Benzoylecgonine Intranasal (ng/mL)'].values.tolist()
cm = df['Cocaine Smoked (ng/mL)'].values.tolist()
bm = df['Benzoylecgonine Smoked (ng/mL)'].values.tolist()

#INTRAVENOUS COCAINE

cIVp1, cIV_up=plt.subplots()
cIVp2, cIV_sol= plt.subplots()
diagonal = np.linspace(0,250,399)

#plt.plot(time,cIV,label='Intravenous Cocaine Graph')
cIVa = cIV[::-1]   
cIV2 = cIVa[:-1]

#line of best fit
a,b = np.polyfit(cIVa[:-1],cIVa[1:],1)
line1 = f"Line of best fit: y = {round(a, 2)}x + {round(b, 2)}"

# bubble sorting algorithm to find the closest value to diagonal (equilibrium point)
for i in range(0,len(diagonal)):
    index_store = 398
    if (abs(diagonal[i]-cIV2[i])) < abs(diagonal[index_store]-cIV2[index_store]):
        index_store = i 

#writing updating function with initial concentration of 0

def cIV_updating(x):
    return (x-b)/a
values = [220]
terminal = np.linspace(0,12,400)

for i in range(len(terminal)):
    values.append(cIV_updating(values[i]))


# ALL GRAPHING
cIV_up.plot(diagonal,diagonal,label='Diagonal')
cIV_up.scatter(cIVa[:-1],cIVa[1:],c='black',s=10,label='Intravenous Cocaine Raw Data Updating Function')
cIV_up.plot(cIVa[:-1],a*np.int64(cIVa[:-1])+b,label=str(line1),color='k',linestyle='dashed')

cIV_sol.plot(terminal,values[:-1],c='k',label='CIV Solution')
cIV_sol.scatter(time,cIV,label='raw data',c='k',s=50)

cIV_sol.legend()
cIV_up.legend()
plt.show()

#Riemann Sum (AUC analysis)



