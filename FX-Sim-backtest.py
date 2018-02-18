# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 17:31:48 2018

@author: Malek
"""
import numpy as np
from matplotlib import pyplot as plt

#### backtesting riskfactors
# [FXSims,FXCcyList,dfDates]
def FXSimBacktestingFull(Sims,CcyList,PlotFlag):
    """ short function to count the number of breaches for all fx pair simulations"""
    
    FXSimCount = []
    for i in range(0,len(CcyList)):
        FX98 = np.percentile(Sims[0,i,:,:],98,axis=0,interpolation='nearest')
        FX02 = np.percentile(Sims[0,i,:,:],2,axis=0,interpolation='nearest')
        FXSpot = Sims[:,i,0,0]
        if PlotFlag:
            plt.clf()
            plt.plot(FX98)
            plt.plot(FX02)
            plt.plot(FXSpot)
            plt.show()
        
        Counter = 0
        for t in range(0,len(FXSpot)):
            if FXSpot[t] > FX98[t] or FXSpot[t] < FX02[t]:
                Counter += 1
        
        FXSimCount.append(Counter)
        
    return FXSimCount

def FXSimBacktestingRolling(Sims,CcyList,PlotFlag):
    """ short function to count the number of breaches for all fx pair simulations"""
    
    FXSimCount = []
    for i in range(4,len(CcyList)):
        for t in range(0,len(Sims[:,i,0,0])-4):
            FX98 = np.percentile(Sims[t,i,:,4],98,axis=0,interpolation='nearest')
            FX02 = np.percentile(Sims[t,i,:,4],2,axis=0,interpolation='nearest')    
            FXSpot = Sims[t+4,i,0,0]
            
            Counter = 0
            if FXSpot > FX98 or FXSpot < FX02:
                Counter += 1
        
        FXSimCount.append(Counter)
        
    return FXSimCount


def TradeBacktestingFull(trade):
    """ function to count number of breaches for a trade """
    
    return None

def TradeBacktestingRolling(trade):
    """ function to count number of breaches for a trade """
    
    return None

def PortfolioBacktestingFull(portfolio):
    """ function to count number of breaches for a portfolio """
    
    return None

def PortfolioBacktestingRolling(portfolio):
    """ function to count number of breaches for a portfolio """
    
    return None

