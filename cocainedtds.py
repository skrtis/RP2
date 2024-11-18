import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict

r = 4
k = 2
terminal_time=10

def logistic_update(x,R=r,K=10):
    return R*(1-x/K)+x

#cobwebbing
def plot_cobweb(f, m_0, steps): #function to cobweb
    
    cobwebs_x=[m_0] #starting x-value
    cobwebs_y=[f(m_0)] #starting y-value
    
    for i in range(1,steps):
        if i%2==1: #if the step i is odd
            cobwebs_x.append(cobwebs_y[i-1]) #append the same values along the diagonal y=x
            cobwebs_y.append(cobwebs_y[i-1])
        elif i%2==0: #if the step i is even
            cobwebs_x.append(cobwebs_y[i-1]) #append the same x-value, but "move up" to the updating curve. 
            cobwebs_y.append(f(cobwebs_x[i]))
            
            
    xvalues=np.linspace(0,max(cobwebs_x), 500) #linspace for plotting
    plt.title("Cobweb Plot")
    plt.plot(xvalues,xvalues, label="y=x") #plot the diagonal y=x
    plt.plot(xvalues, f(xvalues), label="updating function") #plot the updating function
            
    plt.plot(cobwebs_x,cobwebs_y) #plot the cobwebs
    for x,y in zip(list(np.around(np.array(cobwebs_x),2)),list(np.around(np.array(cobwebs_y),2))): #add some labels
        label = f"({x},{y})"
        #uncomment the following code to see point labels. It gets too messy with a lot of steps
        #plt.annotate(label, # this is the text
                    #(x,y), # these are the coordinates to position the label
                     #textcoords="offset points", # how to position the text
                     #xytext=(0,10), # distance from text to points (x,y)
                     #ha='center') # horizontal alignment can be left, right or center
    plt.legend()
    plt.draw()  
    
    plt.figure() 
    
    

    
    distinct_solutions=list(OrderedDict.fromkeys(cobwebs_x))#get the x-components of the cobwebs to plot the solution. This orders them and deletes repeats
    plt.title("Plot of the Cobweb x-coordinates")
    plt.plot(list(range(len(distinct_solutions))),distinct_solutions) #plots the solution set. 
    plt.draw()

    plt.show()

plot_cobweb(logistic_update,1,10)


df = pd.read_excel('raw.xlsx')

