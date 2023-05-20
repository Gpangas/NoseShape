#--------------------------------Wet Surface Calculation-------------------------------
#The Current program was developed to calcute the surface wetted area for four diferent
#nose cone shape (Power series,Elliptical, Van Karman and LV-Haack)
#--------------------------------------------------------------------------------------

#Library-------------------------------------------------------------------------------
import numpy as np
import sympy as smp
import scipy as sp
from scipy.integrate import quad
import matplotlib.pyplot as plt

#Input Data-----------------------------------------------------------------------------
R = 70 #Radius (mm) 
F = np.arange(0.1,5,0.1) #Fineness

A_w=np.empty((len(F)))
A_w_parabolla=np.empty((len(F)))
A_w_elliptical=np.empty((len(F)))
A_w_VonKarman=np.empty((len(F)))
A_w_LV_Haack=np.empty((len(F)))


#Determine the wetted area of a given nose cone shape for a range of fineness
#It's used the formula A_w = 2*pi*f(x)*sqrt(1+(f(x)')^2), so:
#f_ - the coordinate function
#dif_f_ - the derivative of the coordinate function
#int _f_ - the function to integrate

for i in range (len(F)):
    L = F[i]*2*R
    #Power Series(Parabolla)-----------------------------------------------------------
    f_parabolla = lambda x: R*np.sqrt(x/L)
    dif_f_parabolla = lambda x: R/(2*np.sqrt(x*L))
    int_parabolla = lambda x: f_parabolla(x)*np.sqrt(1+pow(dif_f_parabolla(x),2))
    A_w_parabolla[i] = 2*np.pi/(10**6)*quad(int_parabolla,0,L)[0]

    #Elliptical------------------------------------------------------------------------
    f_elliptical = lambda x: R*np.sqrt(1-(x/L)**2)
    dif_f_elliptical = lambda x: -R*x/(L*np.sqrt(L**2-x**2))
    int_elliptical = lambda x: f_elliptical(x)*np.sqrt(1+pow(dif_f_elliptical(x),2))
    A_w_elliptical[i] = 2*np.pi/(10**6)*quad(int_elliptical,0,L)[0]

    #Von Karman------------------------------------------------------------------------
    theta = lambda x: np.arccos(1-2*x/L)
    f_VonKarman = lambda x: R/np.sqrt(np.pi)*np.sqrt(theta(x)-np.sin(2*theta(x))/2)
    dif_f_VonKarman = lambda x: (R/np.sqrt(np.pi))*(L+2*x)/(4*np.sqrt(L*x*(L-x))*np.sqrt(L*theta(x)-np.sqrt(x*(L-x))))
    int_VonKarman = lambda x: f_VonKarman(x)*np.sqrt(1+pow(dif_f_VonKarman(x),2))
    A_w_VonKarman[i] = 2*np.pi/(10**6)*quad(f_VonKarman,0,L)[0]

    #LV-Haack------------------------------------------------------------------------
    theta = lambda x: np.arccos(1-2*x/L)
    f_LV_Haack = lambda x: R/np.sqrt(np.pi)*np.sqrt(theta(x)-np.sin(2*theta(x))/2+(1/3)*(np.sin(theta(x)))**2)
    dif_f_LV_Haack = lambda x: (R/np.sqrt(np.pi))/(2*np.sqrt(theta(x)+np.sin(2*theta(x))/2+(1/3)*np.tan(theta(x))**2)*((1-theta(x))/np.sqrt(x*(L-x))+4*L**2*np.sqrt(x*(L-x))/(L-2*x)**4))
    int_LV_Haack = lambda x: f_LV_Haack(x)*np.sqrt(1+pow(dif_f_LV_Haack(x),2))
    A_w_LV_Haack[i] = 2*np.pi/(10**6)*quad(f_LV_Haack,0,L)[0]

    #Confirmation code--------------------------------------------------------------
    #y = smp.symbols('y', dtype=int)
    #f = R*smp.sqrt(y/L)
    #dif_f = smp.diff(f,y)
    #int = f*smp.sqrt(1+pow(dif_f,2))
    #A_w[i] = 2*np.pi*smp.integrate(int,(y,0,L))
    #print( A_w[i], A_w_parabolla[i])

    #Output data 
    print('For ',L,' the wet area is ',A_w_parabolla[i],A_w_elliptical[i],A_w_VonKarman[i],A_w_LV_Haack[i])


#Plot Graph Wet Area vs Fineness---------------------------------------------------
plt.plot((F), A_w_parabolla, label = "Parabolla")
plt.plot((F), A_w_elliptical, label = "Elliptical")
plt.plot((F), A_w_VonKarman, label = "Von Karman")
plt.plot((F), A_w_LV_Haack, label = "LV-Haack")
  
plt.xlabel('Finenesse')
plt.ylabel('Wet area, [m^2]')
plt.title('Graph Wet Area vs Fineness')
  
plt.legend()
plt.show()


#Plot the nose shapes for R=1 an F=3-----------------------------------------------

n=1000
R_i = 1
F_i = 3

x_i=np.empty(n)
y_parabolla=np.empty(n)
y_elliptical=np.empty(n)
y_VonKarman=np.empty(n)
y_LV_Haack=np.empty(n)

theta_i = lambda x: np.arccos(1-2*x/F_i)

for i in range (n):
    x_i[i] = ((i)/n)*F_i
    y_parabolla[i] = R_i*np.sqrt(x_i[i]/F_i)
    y_elliptical[i] = R_i*np.sqrt(1-((F_i-x_i[i])/F_i)**2)
    y_VonKarman[i] = R_i/np.sqrt(np.pi)*np.sqrt(theta_i(x_i[i])-np.sin(2*theta_i(x_i[i]))/2)
    y_LV_Haack[i] = R_i/np.sqrt(np.pi)*np.sqrt(theta_i(x_i[i])-np.sin(2*theta_i(x_i[i]))/2+(1/3)*(np.sin(theta_i(x_i[i])))**2)

plt.plot(x_i, y_parabolla, label = "Parabolla")
plt.plot(x_i, y_elliptical, label = "Elliptical")
plt.plot(x_i, y_VonKarman, label = "Von Karman")
plt.plot(x_i, y_LV_Haack, label = "LV-Haack")
  
plt.xlabel('Length')
plt.ylabel('Radius')
plt.title('Graph Radius vs Length')
  
plt.legend()
plt.show()
