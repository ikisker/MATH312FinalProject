from numpy import *
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
import IPython
import wavio
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt 
from scipy.stats import linregress
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

def WanderersPath():
    def f(t,y):
        a,b,c,d,e = y
        ap = b*c+d*e-2*b*e
        bp = c*d+e*a-2*c*a
        cp = d*e+a*b-2*d*b
        dp = e*a+b*c-2*e*c
        ep = a*b+c*d-2*a*d
        return array([ap,bp,cp,dp,ep])
    ic = array([2,-2,0,-1,4]) / 5
    t = linspace(0,200,3000)
    sol = solve_ivp(f,[0,t.max()],ic,t_eval = t,atol=1e-9,rtol=1e-9)
    fig,ax = plt.subplots()
    ax.plot(sol.y[0,:],sol.y[1,:])
    ax.plot(cos(arange(629)/100),sin(arange(629)/100),'k-')
    ax.grid()
    ax.set_aspect("equal")
    ax.set_title("A wanderer's path");w
    plt.show()

def Chaos():
    def Lorentz(t,a): 
        u = a[0] # unpack x
        v = a[1] # unpack y 
        w = a[2] # unpack z
        dotu = 16*(v-u)
        dotv = u*(45.6-20*w)-v
        dotw = 5*u*v - 4*w
        return array([dotu,dotv,dotw]) # return vector of derivatives

    ic = array([-2.3,-1.9,2.2]) # I made this up, with no intuition
    tend = 40 # simulation end time
    t = linspace(0,tend,1000) # times when we want to know solutions
    sol = solve_ivp(Lorentz,[0,tend],ic,t_eval=t,atol=1e-9,rtol=1e-9) # :) :) :) 
    assert sol.success 
    u = sol.y[0,:] # extract x(t)
    v = sol.y[1,:] # extract y(t)
    w = sol.y[2,:] # extract z(t)

    fig,ax = plt.subplots(figsize=(10,4))
    ax.plot(t,u,label="u")
    ax.plot(t,v,label="v")
    ax.plot(t,w,label="w")
    ax.set_xlabel("Time")
    ax.set_ylabel("Lorentz components")
    ax.set_title("Chaos, beautiful chaos")
    ax.legend()
    ax.grid();
    plt.show()


def Chaos2():
    def Lorentz(t,a): 
        u = a[0] # unpack x
        v = a[1] # unpack y 
        w = a[2] # unpack z
        dotu = 16*(v-u)
        dotv = u*(45.6-20*w)-v
        dotw = 5*u*v - 4*w
        return array([dotu,dotv,dotw]) # return vector of derivatives

    ic = array([-2.3,-1.9,2.2]) # I made this up, with no intuition
    tend = 20000 # simulation end time
    t = linspace(0,tend,12*44100) # times when we want to know solutions
    sol = solve_ivp(Lorentz,[0,tend],ic,t_eval=t,atol=1e-9,rtol=1e-9) # :) :) :) 
    assert sol.success 
    u = sol.y[0,:] # extract u(t)
    fig,ax = plt.subplots(figsize=(12,4))
    ax.plot(arange(30000)/44100,u[:30000],label="u")
    ax.set_xlabel("Time")
    ax.set_ylabel("Amplitude")
    ax.set_title("An audio signal?")
    ax.legend()
    ax.grid();
    # u = sol.y[0,:] # extract x(t)
    # v = sol.y[1,:] # extract y(t)
    # w = sol.y[2,:] # extract z(t)

    # fig,ax = plt.subplots(figsize=(10,4))
    # ax.plot(t,u,label="u")
    # ax.plot(t,v,label="v")
    # ax.plot(t,w,label="w")
    # ax.set_xlabel("Time")
    # ax.set_ylabel("Lorentz components")
    # ax.set_title("Chaos, beautiful chaos")
    # ax.legend()
    # ax.grid();
    plt.show()



Chaos2()
