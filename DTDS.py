import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Data.csv')

time = np.array(df['Time'].values.tolist())
cIV = np.array(df['Cocaine Intravenous (ng/mL)'].values.tolist())
bIV = np.array(df['Benzoylecgonine Intravenous (ng/mL)'].values.tolist())
cIN = np.array(df['Cocaine Intranasal (ng/mL)'].values.tolist())
bIN = np.array(df['Benzoylecgonine Intranasal (ng/mL)'].values.tolist())
cSM = np.array(df['Cocaine Smoked (ng/mL)'].values.tolist())
bSM = np.array(df['Benzoylecgonine Smoked (ng/mL)'].values.tolist())

cIN_adjusted = cIN[np.argmax(cIN):]
time_adjusted = time[np.argmax(cIN):]

def cocaine_analysis(time, cocaine_data, route_name):
    # Reversal and Pre-Processing 
    cocaine_rev = cocaine_data[::-1]
    cocaine2 = cocaine_rev[:-1]

    # Line of best fit
    a, b = np.polyfit(cocaine_rev[:-1], cocaine_rev[1:], 1)
    line1 = f"Line of best fit: y = {round(a, 2)}x + {round(b, 2)}"

    # Updating Function based on Best Fit 
    def updating(x):
        return (x - b) / a

    # Updating concentrations with an initial value
    values = [max(cocaine_data)]  
    terminal = np.linspace(0, 12, 400)
    
    for i in range(len(terminal)):
        values.append(updating(values[i]))

    # Graphs for analysis
    diagonal = np.linspace(0, 250, 399)
    
    # Plotting figures
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax_rsum = plt.subplots(1, 3, figsize=(15, 5))

    # Plot diagonal and raw data with line of best fit
    ax1.set_title("Updating Function and Line of Best Fit")
    ax1.plot(diagonal, diagonal, label='Diagonal')
    ax1.scatter(cocaine_rev[:-1], cocaine_rev[1:], c='black', s=10, label=f'{route_name} Cocaine Raw Data Updating Function')
    ax1.plot(cocaine_rev[:-1], a * np.int64(cocaine_rev[:-1]) + b, label=str(line1), color='k', linestyle='dashed')
    
    # Plot the solution from updating function
    ax2.set_title(f'{route_name} Solution')
    ax2.plot(terminal, values[:-1], c='k', label=f'{route_name} Solution')
    ax2.scatter(time, cocaine_data, label='Raw Data', c='k', s=50)
    ax2.legend()
    ax1.legend()

    # Trapezoidal Riemann Sum (AUC analysis)
    fig3.suptitle(f'{route_name} AUC Analysis')
    
    # Left sum
    ax_rsum[0].set_title('Left Sum')
    time_left = time[:-1]
    cocaine_left = cocaine_data[:-1]
    ax_rsum[0].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[0].plot(time_left, cocaine_left, c='r')
    ax_rsum[0].bar(time_left, cocaine_left, width=0.03, alpha=0.5, align='edge', edgecolor='k')

    # Middle sum
    ax_rsum[1].set_title('Middle Sum')
    time_mid = (time[:-1] + time[1:]) / 2
    cocaine_mid = (cocaine_data[:-1] + cocaine_data[1:]) / 2
    ax_rsum[1].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[1].plot(time_mid, cocaine_mid, 'b.', markersize=10)
    ax_rsum[1].bar(time_mid, cocaine_mid, width=0.03, alpha=0.5, edgecolor='k')

    # Right sum
    ax_rsum[2].set_title('Right Sum')
    time_right = time[1:]
    cocaine_right = cocaine_data[1:]
    ax_rsum[2].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[2].plot(time_right, cocaine_right, c='r')
    ax_rsum[2].bar(time_right, cocaine_right, width=0.03, alpha=0.5, align='edge', edgecolor='k')

    plt.show()


cocaine_analysis(time_adjusted,cIN_adjusted,'Cut Intranasal')