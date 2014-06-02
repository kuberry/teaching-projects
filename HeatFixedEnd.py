import numpy
import scipy
import matplotlib
import matplotlib.pyplot as plt
import time

# This code performs a 1-d simulation of the distribution
# of heat over time via the finite difference method using
# an explicit, forward, first-order accurate time
# scheme and a centered, second-order spatial
# discretization.
# 
# Author: Paul Kuberry, 2014


def initialCondition(x):
    return numpy.maximum(0,-numpy.absolute(10*x-5)+.5)

# Spatial discretization
numPoints = 101;

# Visual Parameters
# How long the movie should play over
displayTime = 3.0;
framesPerSecond = 30;

# Thermal diffusivity
alpha = 0.1;
finalTime = 0.1;

# Domain Dimensions
xMin = 0;
xMax = 1;

# Calculate mesh size and time step size 
h = float(xMax-xMin)/(numPoints-1);
# Time step needed to be stable
dt = 2*(h**2)/(alpha*numpy.pi**2); 
numTimeSteps = float(finalTime)/dt;
numTimeSteps = numpy.ceil(numTimeSteps); 
# Recompute time steps so that we come out evenly at final time
dt = float(finalTime)/numTimeSteps;

# Initialize vectors used to store solutions
x = numpy.linspace(xMin,xMax,numPoints);
# Solution at time step n
y_n = initialCondition(x);
# Solution at time step n - 1
y_nm1 = numpy.empty(x.size);
# The simulation begins at rests
y_nm1[:] = y_n;

# Prepare Matplotlib for interactive plotting
matplotlib.use('TkAgg');
plt.ion();
fig = plt.figure();
# Plot the initial condition
ax, = plt.plot(x,y_n,lw=2);
plt.xlabel('x-axis');
plt.ylabel('y-axis');
plt.title('Heat Simulation');
# Reset the y-axis limits based on the assumption that the wave will travel
plt.ylim((-numpy.max(numpy.absolute(y_n)),numpy.max(numpy.absolute(y_n))))
plt.show();
# Determine how often to plot the solution so that the movie lasts $displayTime seconds
# Note: If the spatial discretization is fine, this may take longer that #displayTime
# seconds to show, and you may wish to comment out the time.sleep() call.
numPlotUpdates = max(1,numpy.floor(numTimeSteps/(displayTime*framesPerSecond)));
plotTimeDelay = displayTime/(numTimeSteps/numPlotUpdates);


for i in xrange(0,int(numTimeSteps),1):
    # current_time = (i+1)*dt;

    # Discretize in time with u_tt ~= (u_n-2u_nm1+u_nm2)/dt^2
    #                and with u_xx ~= (u_nm1(x+h)-2u_nm1(x)+u_nm1(x-h))/h^2
    #
    # On the end points we need a second order accurate difference
    # formula     u_xx ~= (-u_nm1(x+2h)+4u_nm1(x+h)-3u_nm1(x))/2h
    # where x is left end point, and
    #
    #                 u_tt = c * u_xx 
    #
    
    for j in xrange(1,int(numPoints)-1,1):
        y_n[j] = y_nm1[j] + dt*alpha*(1./(h**2))*(y_nm1[j+1]-2*y_nm1[j]+y_nm1[j-1]);

    # Update previous two time step solutions    
    y_nm1[:] = y_n;

    # Plotting logic 
    if numpy.mod(i,numpy.ceil(numTimeSteps/(displayTime*framesPerSecond)))==0:
        ax.set_ydata(y_n);
        fig.canvas.draw()
        #time.sleep(plotTimeDelay) 
        #time.sleep(numpy.max(1./(displayTime*framesPerSecond),1./numTimeSteps))
        
    
    
        
    
    

    























# bx = plt.gca();
# bx.relim();
# bx.autoscale();
# Trigger emacs to run this script using the "compile" command
# ;;; Local Variables: ***
# ;;; compile-command: "python HeatFixedEnd.py" ***
# ;;; end: ***
