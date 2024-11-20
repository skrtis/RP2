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

# First updating function graph
#plt.plot(time,cIV,label='Intravenous Cocaine Graph')
diagonal = np.linspace(0,250,399)
plt.scatter(cIV[:-1],cIV[1:],label='Intravenous Cocaine Updating Function')
a,b = np.polyfit(cIV[:-1],cIV[1:],1)
plt.plot(cIV[:-1],a*np.int64(cIV[:-1])+b)
#plt.plot(diagonal,diagonal,label='Diagonal')
plt.show()

cIV2 = cIV[:-1]
# bubble sorting algorithm to find the closest value to diagonal (equilibrium point)
for i in range(0,len(diagonal)):
    index_store = 398
    if (abs(diagonal[i]-cIV2[i])) < abs(diagonal[index_store]-cIV2[index_store]):
        index_store = i 
print(cIV2[i])

