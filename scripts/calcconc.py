import math
import re
import csv
import numpy as np
from getRateCoefficients import getRateCoefficients

def calcconc(ratecoefficients, o2init, fuelinit,time):
    """returns values of concentration for HO2, H2O2 and formaldehyde at the specified time.  
    all inputs shoudl be in mol, m3 and s and the outputs are the same units. 
    """ 
    #initiate constants
    kho2self=ratecoefficients[0][0]
    kh2o2uni=ratecoefficients[1][0]        
    kho2fuel=(ratecoefficients[2][0])*fuelinit
    ko2fuel =ratecoefficients[3][0]*o2init*fuelinit
    kho2formal=ratecoefficients[4][0]
    exp_eigval=[ (-kh2o2uni/2-(kh2o2uni*(kh2o2uni+8*kho2fuel))**(1./2)/2),
                            ((kh2o2uni*(kh2o2uni+8*kho2fuel))**(1./2)/2-kh2o2uni/2)]
    exp_eigvec=[[(exp_eigval[0]/2/kh2o2uni),(exp_eigval[1]/2/kh2o2uni)],
                    [1,1]]

    #phase 1
    time_knee=np.log((kh2o2uni/4/kho2self*kho2fuel)/ko2fuel)/exp_eigval[1]
    print time_knee
    HO2_value_exp=2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]/exp_eigval[0]-exp_eigvec[0][0]/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*time)/exp_eigval[0]+exp_eigvec[0][0]*np.exp(exp_eigval[1]*time)/exp_eigval[1])
    H2O2_value_exp=((ko2fuel*math.exp(-(time*(kh2o2uni + (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2)*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) - 2*ko2fuel*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) + kh2o2uni*ko2fuel*math.exp(-(time*(kh2o2uni - (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2) + ko2fuel*math.exp(-(time*(kh2o2uni - (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2)*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) - kh2o2uni*ko2fuel*math.exp(-(time*(kh2o2uni + (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2))/(2*kh2o2uni*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))
    formal_value_exp=time*ko2fuel+kho2fuel*(2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]*time/exp_eigval[0]-exp_eigvec[0][0]*time/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*time)/exp_eigval[0]**2+exp_eigvec[0][0]*np.exp(exp_eigval[1]*time)/exp_eigval[1]**2))+2*kh2o2uni*(2*ko2fuel/(1/exp_eigvec[0][1]-1/exp_eigvec[0][0])*(time/exp_eigval[0]-time/exp_eigval[1]-np.exp(exp_eigval[0]*time)/exp_eigval[0]**2+np.exp(exp_eigval[1]*time)/exp_eigval[1]**2))
    formal_value_intercept=0*ko2fuel+kho2fuel*(2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]*0/exp_eigval[0]-exp_eigvec[0][0]*0/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*0)/exp_eigval[0]**2+exp_eigvec[0][0]*np.exp(exp_eigval[1]*0)/exp_eigval[1]**2))+2*kh2o2uni*(2*ko2fuel/(1/exp_eigvec[0][1]-1/exp_eigvec[0][0])*(0/exp_eigval[0]-0/exp_eigval[1]-np.exp(exp_eigval[0]*0)/exp_eigval[0]**2+np.exp(exp_eigval[1]*0)/exp_eigval[1]**2))
    formal_value_exp=formal_value_exp-formal_value_intercept
#     print str(formal_value_intercept)+'    '+str(formal_value_exp)+'      '+str(time)
    if time_knee>time:
        HO2_value_exp=2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]/exp_eigval[0]-exp_eigvec[0][0]/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*time)/exp_eigval[0]+exp_eigvec[0][0]*np.exp(exp_eigval[1]*time)/exp_eigval[1])      
#         print str(HO2_value_exp)+'function values'
        return [HO2_value_exp,H2O2_value_exp,formal_value_exp,HO2_value_exp,H2O2_value_exp,formal_value_exp]
    
    #determining constants for phase 2
    HO2_value_knee=2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]/exp_eigval[0]-exp_eigvec[0][0]/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*time_knee)/exp_eigval[0]+exp_eigvec[0][0]*np.exp(exp_eigval[1]*time_knee)/exp_eigval[1])
    H2O2_value_knee=((ko2fuel*math.exp(-(time_knee*(kh2o2uni + (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2)*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) - 2*ko2fuel*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) + kh2o2uni*ko2fuel*math.exp(-(time_knee*(kh2o2uni - (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2) + ko2fuel*math.exp(-(time_knee*(kh2o2uni - (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2)*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2) - kh2o2uni*ko2fuel*math.exp(-(time_knee*(kh2o2uni + (kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))/2))/(2*kh2o2uni*(kh2o2uni*(kh2o2uni + 8*kho2fuel))**(1./2)))
    formal_value_knee=time_knee*ko2fuel+kho2fuel*(2*ko2fuel/(exp_eigvec[0][0]-exp_eigvec[0][1])*(exp_eigvec[0][1]*time_knee/exp_eigval[0]-exp_eigvec[0][0]*time_knee/exp_eigval[1]-exp_eigvec[0][1]*np.exp(exp_eigval[0]*time_knee)/exp_eigval[0]**2+exp_eigvec[0][0]*np.exp(exp_eigval[1]*time_knee)/exp_eigval[1]**2))+kh2o2uni*(2*ko2fuel/(1/exp_eigvec[0][1]-1/exp_eigvec[0][0])*(time_knee/exp_eigval[0]-time_knee/exp_eigval[1]-np.exp(exp_eigval[0]*time_knee)/exp_eigval[0]**2+np.exp(exp_eigval[1]*time_knee)/exp_eigval[1]**2))
    
    HO2_slope_lin=(kho2fuel)/2*kh2o2uni/kho2self
    H2O2_slope_lin=HO2_slope_lin*np.sqrt(kho2self/kh2o2uni)
    constant_HO2_knee=HO2_value_knee-(HO2_slope_lin)*time_knee
    constant_H2O2_knee=constant_HO2_knee*np.sqrt(kho2self/kh2o2uni)
    formal_value_lin_indefinite_int=(ko2fuel+constant_H2O2_knee**2*2*kh2o2uni+kho2fuel*constant_HO2_knee)*time_knee+(kh2o2uni/4/kho2self*kho2fuel**2+2*kh2o2uni*constant_H2O2_knee*0.5*kho2fuel*np.sqrt(kh2o2uni/kho2self))*time_knee**2+2/3*kh2o2uni*1/4*kho2fuel**2*kh2o2uni/kho2self*time_knee**3
    constant_formal_knee=formal_value_knee-formal_value_lin_indefinite_int
    
    #phase 2 calculations
    H2O2_value_lin=(H2O2_slope_lin*time+constant_H2O2_knee)**2
    HO2_value_lin=HO2_slope_lin*time+constant_HO2_knee
    formal_value_lin=(ko2fuel+constant_H2O2_knee**2*2*kh2o2uni+kho2fuel*constant_HO2_knee)*time+(kh2o2uni/4/kho2self*kho2fuel**2+2*kh2o2uni*constant_H2O2_knee*0.5*kho2fuel*np.sqrt(kh2o2uni/kho2self))*time**2+2/3*kh2o2uni*1/4*kho2fuel**2*kh2o2uni/kho2self*time**3+constant_formal_knee    

    #using tangent approximation
    constant_tan=(np.arctan(kho2formal*H2O2_value_knee/kho2fuel))**2-.5*np.sqrt(kh2o2uni*kho2fuel*kho2formal/kho2self)*time_knee
    H2O2_value_tan=kho2fuel/kho2formal*(np.tan(.5*np.sqrt(kh2o2uni*kho2fuel*kho2formal/kho2self)*time+constant_tan))**2
    HO2_value_tan=np.sqrt(kh2o2uni/kho2self*H2O2_value_tan)
    formal_value_tan=H2O2_value_tan
#    return [HO2_value_lin,H2O2_value_lin,formal_value_lin,HO2_value_tan,H2O2_value_tan,formal_value_tan]
    return [HO2_value_tan,H2O2_value_tan,formal_value_tan]
