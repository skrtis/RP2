import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Data.csv')

time = np.array(df['Time'].values.tolist())
cIV = np.array(df['Cocaine Intravenous (ng/mL)'].values.tolist())
bIV = np.array(df['Benzoylecgonine Intravenous (ng/mL)'].values.tolist())
cIN = np.array(df['Cocaine Intranasal (ng/mL)'].values.tolist())
bIN = np.array(df['Benzoylecgonine Intranasal (ng/mL)'].values.tolist())
cm = np.array(df['Cocaine Smoked (ng/mL)'].values.tolist())
bm = np.array(df['Benzoylecgonine Smoked (ng/mL)'].values.tolist())

# GAMMA DISTRIBUTION TESTING
from scipy.stats import gamma

x = time
y = cIN

# Step 1: Fit a gamma distribution to the concentration data (y)
shape, loc, scale = gamma.fit(y)

# Step 2: Generate the gamma distribution curve based on the time points (x)
# Use the time points as the range for the gamma curve
pdf_fitted = gamma.pdf(x, shape, loc, scale)

# Step 3: Plot the original data and the fitted gamma curve
plt.plot(x, y, 'bo-', label='Original Data')  # Original time vs concentration
plt.plot(x, pdf_fitted, 'r--', label=f'Gamma fit\nShape={shape:.2f}, Loc={loc:.2f}, Scale={scale:.2f}')  # Fitted gamma curve
plt.legend()
plt.xlabel('Time')
plt.ylabel('Plasma Concentration (ng/mL)')
plt.title('Gamma Distribution Fit to Plasma Concentration Data')
plt.show()
