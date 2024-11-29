import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Data.csv')
time_data = np.array(df['Time (h)'].values.tolist()) #general time dataset
cIV = np.array(df['Cocaine Intravenous (ng/mL)'].values.tolist()) # Intravenous cocaine data
cIN = np.array(df['Cocaine Intranasal (ng/mL)'].values.tolist()) # Intranasal cocaine data
cSM = np.array(df['Cocaine Smoked (ng/mL)'].values.tolist()) # Smoked cocaine data

def differentiate(time, cocaine_data, route_name): #this function is primarily used for differentiating the raw data for the elimination rate
    # Calculate the rate of change (first derivative)
    rate_of_change = np.diff(cocaine_data) / np.diff(time) #dividing the differences of the cocaine data by the differences of the time data

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

def cocaine_exp_analysis(time, cocaine_data, route_name): #this is the general analysis function, which includes the updating function, solution, repeated doses, and the inverse updating function

    a,b = np.polyfit(cocaine_data[:-1], cocaine_data[1:], 1) #polyfit is used to fit the data to a linear model
    line1 = f"Line of best fit: $c_{{t+0.03}} = {round(a, 2)}c_t + {round(b, 2)}$" #this is a string with the equation of the line of best fit
    
    predicted = a * cocaine_data[:-1] + b #this creates predicted values for the data based on the updating function
    residual = np.sum((cocaine_data[1:] - predicted) ** 2) #this calculates the residuals
    squares = np.sum((cocaine_data[1:] - np.mean(cocaine_data[1:])) ** 2) #this calculates the sum of squares
    r_squared = 1 - (residual / squares) #this calculates the r squared value
    line1 += f", $R^2 = {round(r_squared, 4)}$" #this adds the r squared value to the string
    print("The r-squared value is: "+str(r_squared)) #this prints the r squared value
   
    def updating(x):
         return a * x + b #this is the updating function, which uses a,b determined by polyfit

    # UPDATING FUNCTION (generates based on updating function)
      # the lists for the repeated function is "terminal" and "values"
    values = [max(cocaine_data)] 
    terminal = np.linspace(0, 12, 401)
    for i in range(len(terminal)):
            values.append(updating(values[i]))
    
    diagonal = np.linspace(min(cocaine_data[1:]), max(cocaine_data[1:]), len(cocaine_data[1:]))

    # Find the index where cocaine_data[1:] is the closest to the corresponding value on the diagonal
    min_diff = np.inf
    closest_index = -1
    for i in range(len(diagonal)):
        diff = abs(cocaine_data[i + 1] - diagonal[i])
        if diff < min_diff:
            min_diff = diff
            closest_index = i

    print(f'Equilibrium Point Using Raw Data: {cocaine_data[closest_index + 1]} ng/mL')
    intersection_x = (b / (1 - a))
    intersection_y = intersection_x

    # REPEATED DOSES UPDATING FUNCTION
        # the lists for the repeated function is "terminal" and "repeated_values"
    repeated_values =[0] #this is the list for the repeated doses
    for i in range(len(terminal)):
        interval = 33
        dose = 250
        if i%interval==0: #interval for repeated doses
            repeated_values.append(updating(repeated_values[i]+dose)) #amount of repeated dose
        else: 
            repeated_values.append(updating(repeated_values[i]))

    # Plotting the Moving Average for Repeated Doses
        # the lists for the repeated function is "moving_time" and "moving_average"
    moving_average = [0]
    counter = []
    moving_time=[0]
    for i in repeated_values:
        if len(counter) < 33:
            counter.append(i)
        else:
            moving_average.append(np.mean(counter))
            moving_time.append(repeated_values.index(i)*0.03)
            counter.pop(0)
            counter.append(i)

    # these are the calculations for three levels of dosages: 250, 200, and 100 ng/mL, and how much is eliminated in a certain time (in this case 20)
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


    # These are the calculations for the equilibria of the repeated doses, one differentiating minima and one using moving average
    rate_of_change = np.diff(mins) / np.diff(newtime)
    equilibria_pos = 0
    for i in rate_of_change:
        if i < 1e-3:
            equilibria_pos = newtime[np.where(rate_of_change == i)[0][0]]
            break
    print("Mins:" + str(equilibria_pos))

    roc2 = np.diff(moving_average) / np.diff(moving_time)
    equilibria_pos2 = 0 
    for i in roc2:
        if i < 1e-3:
            equilibria_pos2 = moving_time[np.where(roc2 == i)[0][0]]
            break
    print("Moving Average:" + str(equilibria_pos2))

    # All Figure Plots
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()

    # Intravenous Updating Plot with Raw Updating Data, Line of Best Fit, Diagonal and Equilibrium Points    
    ax1.set_title(f"{route_name} Updating Function, Diagonal and Line of Best Fit")
    ax1.plot(diagonal, diagonal, label='Diagonal')
    ax1.plot(cocaine_data[closest_index],cocaine_data[closest_index+1],'ro', label='Equilibrium Point from Raw Data')
    ax1.plot(intersection_x, intersection_y, 'bo', label='Equilibrium Point from Updating Function')
    ax1.scatter(cocaine_data[:-1], cocaine_data[1:], c='black', s=10, label=f'{route_name} Cocaine Raw Data Updating Function')
    ax1.set_ylabel('Plasma Cocaine Concentration $c_{{t+1}}$ (ng/mL)')
    ax1.set_xlabel('Plasma Cocaine Concentration $c_t$ (ng/mL)')
    ax1.legend()

    # Intravenous Solutions: Raw Data and Generated from Updating Function
    ax2.set_title(f'{route_name} Raw Data and Solution Determined by Updating Function')
    ax2.set_title(f'{route_name} Solution From 250 ng/mL to 20 ng/mL')
    ax2.plot(terminal, values[:-1], c='k', label=f'{route_name} Solution')
    ax2.plot(terminal[:vstore+1],iconc1[:vstore+1],'ro',label='250 ng/mL Initial Concentration')
    ax2.scatter(time, cocaine_data, label='Raw Data', c='k', s=50)
    ax2.set_xlabel('Time (h)')
    ax2.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax2.legend()

    # Intravenous Solution for Repeated Doses: Minima and Moving Averages Plotted
    ax3.set_title(f'{route_name} Solution With Repeated Doses of {dose} ng/mL at {round(interval*0.03,2)} h Intervals')
    ax3.plot(terminal, repeated_values[:-1], c='k', label=f'{route_name} Solution for Repeated Doses')
    ax3.plot(newtime, mins, label='Min Values', c='k')
    ax3.plot(moving_time, moving_average, label='Moving Average', c='r')
    ax3.set_xlabel('Time (h)')
    ax3.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax3.legend()
    
    # Backwards DTDS Solution from Inverse Updating Function
    ax4.set_title(f'{route_name} Solution Determined by Inverse Updating Function')
    values_rev = values[::-1]
    ax4.plot(terminal, values_rev[:-1], c='k', label=f'{route_name} Inverse Solution')
    ax4.set_xlabel('Time (h)')
    ax4.set_ylabel('Plasma Cocaine Concentration (ng/mL)')
    ax4.legend()

    plt.show()

    # Table of Values for Repeated Doses
    tov = {'Time (h)': terminal, 'Repeated Values (ng/mL)': repeated_values[:-1]}
    df_table = pd.DataFrame(tov)
    print(df_table.to_string())
 
def rsum_trapezoid(time, cocaine_data,route_name): # Function for Computing the Trapezoidal Riemann Sum
    return np.sum(np.diff(time) * ((cocaine_data[:-1] + cocaine_data[1:]) / 2)) 

def riemann_sum_analysis(time, cocaine_data, route_name): #Riemann Sum Analysis
    fig3, ax_rsum = plt.subplots(1, 3, figsize=(15, 5))
    fig3.suptitle(f'{route_name} AUC Analysis')

    rsum = rsum_trapezoid(time, cocaine_data,route_name)
    print(f'Riemann Sum for {route_name} Route of Administration: '+str(rsum))

    #Plot of Left sum
    ax_rsum[0].set_title('Left Sum')
    time_left = time[:-1]
    cocaine_left = cocaine_data[:-1]
    ax_rsum[0].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[0].plot(time_left, cocaine_left, c='r')
    ax_rsum[0].bar(time_left, cocaine_left, width=0.03, alpha=0.5, align='edge', edgecolor='k')
    ax_rsum[0].set_xlabel('Time (h)')
    ax_rsum[0].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    # Plot of Middle sum
    ax_rsum[1].set_title('Middle Sum')
    time_mid = (time[:-1] + time[1:]) / 2
    cocaine_mid = (cocaine_data[:-1] + cocaine_data[1:]) / 2
    ax_rsum[1].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[1].plot(time_mid, cocaine_mid, 'b.', markersize=10)
    ax_rsum[1].bar(time_mid, cocaine_mid, width=0.03, alpha=0.5, edgecolor='k')
    fig3.text(0.5,0.01,f'Riemann Sum for {route_name} Route of Administration: '+str(np.sum(np.diff(time) * ((cocaine_data[:-1] + cocaine_data[1:]) / 2))), ha='center')
    ax_rsum[1].set_xlabel('Time (h)')
    ax_rsum[1].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    # Plot of Right sum
    ax_rsum[2].set_title('Right Sum')
    time_right = time[1:]
    cocaine_right = cocaine_data[1:]
    ax_rsum[2].scatter(time, cocaine_data, label='Raw Data', c='k')
    ax_rsum[2].plot(time_right, cocaine_right, c='r')
    ax_rsum[2].bar(time_right, cocaine_right, width=0.03, alpha=0.5, align='edge', edgecolor='k')
    ax_rsum[2].set_xlabel('Time (h)')
    ax_rsum[2].set_ylabel('Plasma Cocaine Concentration (ng/mL)')

    plt.show()


def run(time,list,route_name): # Run Function that Helps Adjust the Data for the Differentiation and Analysis Functions
    differentiate(time[np.argmax(list):],list[np.argmax(list):],route_name)
    cocaine_exp_analysis(time[np.argmax(list):],list[np.argmax(list):],route_name)
    plt.show()


run(time_data,cIV,'Intravenous') #Intraveous Analysis
#Riemann Sum Analysis for Intravenous, Intranasal and Smoked Routes of Administration
riemann_sum_analysis(time_data,cIV,'Intravenous') 
riemann_sum_analysis(time_data,cIN,'Intranasal')
riemann_sum_analysis(time_data,cSM,'Smoked')