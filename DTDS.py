import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Data.csv')
time_data = np.array(df['Time (h)'].values.tolist())
cIV = np.array(df['Cocaine Intravenous (ng/mL)'].values.tolist())   
cIN = np.array(df['Cocaine Intranasal (ng/mL)'].values.tolist())
cSM = np.array(df['Cocaine Smoked (ng/mL)'].values.tolist())

def differentiate(time, cocaine_data, route_name):
    # Calculate the rate of change (first derivative)
    rate_of_change = np.diff(cocaine_data) / np.diff(time)

    # Adjust time array to match the length of the rate_of_change array
    time_diff = time[:-1] + np.diff(time) / 2

    # Plot the rate of change
    fig, ax = plt.subplots()
    ax.set_title(f'Rate of Elimination of Raw {route_name} Data')
    ax.plot(time_diff, rate_of_change, label=f'{route_name} Rate of Elimination', color='k')
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Rate of Change (ng/mL per hour)')
    ax.legend()

    plt.show()

def cocaine_exp_analysis(time, cocaine_data, route_name):

    a,b = np.polyfit(cocaine_data[:-1], cocaine_data[1:], 1)
    line1 = f"Line of best fit: $c_{{t+0.03}} = {round(a, 2)}c_t + {round(b, 2)}$"
    
    predicted = a * cocaine_data[:-1] + b
    residual = np.sum((cocaine_data[1:] - predicted) ** 2)
    squares = np.sum((cocaine_data[1:] - np.mean(cocaine_data[1:])) ** 2)
    r_squared = 1 - (residual / squares)
    line1 += f", $R^2 = {round(r_squared, 4)}$"
    print(r_squared)
   
    def updating(x):
         return a * x +b

    # Updating concentrations with an initial value
    values = [max(cocaine_data)]  
    terminal = np.linspace(0, 12, 401)
    repeated_values =[0]  
    
    for i in range(len(terminal)):
            values.append(updating(values[i]))
        
    # REPEATED DOSES    
    for i in range(len(terminal)):
        interval = 33
        dose = 400
        if i%interval==0: #interval for repeated doses
            repeated_values.append(updating(repeated_values[i]+dose)) #amount of repeated dose
        else: 
            repeated_values.append(updating(repeated_values[i]))

    iconc1 = [250]
    iconc2 = [200]
    iconc3 = [100]
    vstore = 0
    for v in range(len(terminal)):
        if updating(iconc1[v])<20:
            vstore=v
            print('Time to reach 20 ng/mL for 250 ng/mL initial concentration:',round(terminal[v],2))
            break
        iconc1.append(updating(iconc1[v]))
        iconc2.append(updating(iconc2[v]))
        iconc3.append(updating(iconc3[v]))

    

    print(round(250-iconc1[-1],2),round(200-iconc2[-1],2),round(100-iconc3[-1],2)) 

    # Finding the Equilibria of the Repeated Doses
    mins=[]
    newtime=[]
    for i in range(len(repeated_values)):
        if i%interval==0:
            newtime.append(i*0.03)
            mins.append(repeated_values[i])


    rate_of_change = np.diff(mins) / np.diff(newtime)
    equilibria_pos = 0
    for i in rate_of_change:
        if i < 1e-3:
            equilibria_pos = newtime[np.where(rate_of_change == i)[0][0]]
            break
    print(equilibria_pos)

    # Adjust the size and intervals of diagonal to match cocaine_data[1:]
    diagonal = np.linspace(min(cocaine_data[1:]), max(cocaine_data[1:]), len(cocaine_data[1:]))

    # Find the index where cocaine_data[1:] is the closest to the corresponding value on the diagonal
    min_diff = np.inf
    closest_index = -1
    for i in range(len(diagonal)):
        diff = abs(cocaine_data[i + 1] - diagonal[i])
        if diff < min_diff:
            min_diff = diff
            closest_index = i

    #print(f'Equilibrium Point Using Raw Data: {cocaine_data[closest_index + 1]} ng/mL')
    intersection_x = (b / (1 - a))
    intersection_y = intersection_x

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    #fig4, ax4 = plt.subplots()

    # Plot diagonal and raw data with line of best fit
    ax1.set_title(f"{route_name} Updating Function, Diagonal and Line of Best Fit")
    ax1.plot(diagonal, diagonal, label='Diagonal')
    ax1.plot(cocaine_data[closest_index],cocaine_data[closest_index+1],'ro', label='Equilibrium Point from Raw Data')
    ax1.plot(intersection_x, intersection_y, 'bo', label='Equilibrium Point from Updating Function')
    ax1.scatter(cocaine_data[:-1], cocaine_data[1:], c='black', s=10, label=f'{route_name} Cocaine Raw Data Updating Function')
    ax1.set_ylabel('Plasma Cocaine Concentration $c_{{t+1}}$ (ng/mL)')
    ax1.set_xlabel('Plasma Cocaine Concentration $c_t$ (ng/mL)')
    ax1.legend()

    # Plot the solution from updating function
    #ax2.set_title(f'{route_name} Raw Data and Solution Determined by Updating Function')
    ax2.set_title(f'{route_name} Solution From 250 ng/mL to 20 ng/mL')
    #ax2.plot(terminal, values[:-1], c='k', label=f'{route_name} Solution')
    ax2.plot(terminal[:vstore+1],iconc1[:vstore+1],'ro',label='250 ng/mL Initial Concentration')
    #ax2.scatter(time, cocaine_data, label='Raw Data', c='k', s=50)
    ax2.set_xlabel('Time (h)')
    ax2.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax2.legend()

    # Plot Solutions
    ax3.set_title(f'{route_name} Solution With Repeated Doses of {dose} ng/mL at {round(interval*0.03,2)} h Intervals')
    ax3.plot(terminal, repeated_values[:-1], c='k', label=f'{route_name} Solution for Repeated Doses')
    ax3.plot(newtime, mins, label='Min Values', c='k')
    # the above line is used to check if the mins are plotted correctly
    ax3.set_xlabel('Time (h)')
    ax3.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax3.legend()

    '''
    ax4.set_title(f'{route_name} Solution Determined by Inverse Updating Function')
    value_rev = values[::-1]
    ax4.plot(terminal, value_rev[:-1], c='k', label=f'{route_name} Inverse Solution')
    ax4.set_xlabel('Time (h)')
    ax4.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax4.legend()
    '''

    plt.show()
    tov = {'Time (h)': terminal, 'Repeated Values (ng/mL)': repeated_values[:-1]}
    df_table = pd.DataFrame(tov)
    print(df_table.to_string())
 
# computing the Riemann Sum
def rsum_trapezoid(time, cocaine_data,route_name):
    return np.sum(np.diff(time) * ((cocaine_data[:-1] + cocaine_data[1:]) / 2)) 


def riemann_sum_analysis(time, cocaine_data, route_name):
    fig3, ax_rsum = plt.subplots(1, 3, figsize=(15, 5))
    fig3.suptitle(f'{route_name} AUC Analysis')
    
    print(f'Riemann Sum for {route_name} Route of Administration: '+str(np.sum(np.diff(time) * ((cocaine_data[:-1] + cocaine_data[1:]) / 2))))

    # Left sum
    ax_rsum[0].set_title('Left Sum')
    time_left = time[:-1]
    cocaine_left = cocaine_data[:-1]
    ax_rsum[0].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[0].plot(time_left, cocaine_left, c='r')
    ax_rsum[0].bar(time_left, cocaine_left, width=0.03, alpha=0.5, align='edge', edgecolor='k')
    ax_rsum[0].set_xlabel('Time (h)')
    ax_rsum[0].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    # Middle sum
    ax_rsum[1].set_title('Middle Sum')
    time_mid = (time[:-1] + time[1:]) / 2
    cocaine_mid = (cocaine_data[:-1] + cocaine_data[1:]) / 2
    ax_rsum[1].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[1].plot(time_mid, cocaine_mid, 'b.', markersize=10)
    ax_rsum[1].bar(time_mid, cocaine_mid, width=0.03, alpha=0.5, edgecolor='k')
    fig3.text(0.5,0.01,f'Riemann Sum for {route_name} Route of Administration: '+str(np.sum(np.diff(time) * ((cocaine_data[:-1] + cocaine_data[1:]) / 2))), ha='center')
    ax_rsum[1].set_xlabel('Time (h)')
    ax_rsum[1].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    # Right sum
    ax_rsum[2].set_title('Right Sum')
    time_right = time[1:]
    cocaine_right = cocaine_data[1:]
    ax_rsum[2].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[2].plot(time_right, cocaine_right, c='r')
    ax_rsum[2].bar(time_right, cocaine_right, width=0.03, alpha=0.5, align='edge', edgecolor='k')
    ax_rsum[2].set_xlabel('Time (h)')
    ax_rsum[2].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    plt.show()


def run(time,list,route_name):
    #differentiate(time[np.argmax(list):],list[np.argmax(list):],route_name)
    cocaine_exp_analysis(time[np.argmax(list):],list[np.argmax(list):],route_name)
    plt.show()


#run(time_data,cIV,'Intravenous')
riemann_sum_analysis(time_data,cIV,'Intravenous')
riemann_sum_analysis(time_data,cIN,'Intranasal')
riemann_sum_analysis(time_data,cSM,'Smoked')