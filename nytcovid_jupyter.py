
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scipy.optimize as optz

import pdb

nyt = pd.read_csv('./nytimes_apr8.csv')

census = pd.read_csv('./uscensus_2010_county_pop_density.csv')

states = pd.read_csv('./us_states.csv', index_col='state')

# Manual changes:
#nyt.loc[nyt.county=='New York City', 'county']='New York'
census.loc[census.shape[0], 'area_name']='New York City, NY'
nyboroughs = ['New York, NY', 'Kings, NY', 'Bronx, NY', 'Queens, NY', 'Richmond, NY']
census.loc[census.area_name=='New York City, NY', 'pop'] = census.loc[census.area_name.isin(nyboroughs), 'pop'].sum()
census.loc[census.area_name=='New York City, NY', 'density'] = census.loc[census.area_name=='New York City, NY', 'pop'] / ( census.loc[census.area_name.isin(nyboroughs), 'pop'] / census.loc[census.area_name.isin(nyboroughs), 'density'] ).sum()

# create new column containing County Name, State abreviation
nyt.loc[:,'abrev'] = list(states.abrev[nyt['state']])
nyt['ctst'] = nyt.county + ', ' + nyt.abrev

# convert dates to date data type:
nyt.date = pd.to_datetime( nyt.date )



class County :
    
    
    def __init__(self, name, date, cases, deaths):
        
        self.name = name
        self.date = [date]
        self.cases= [cases]
        self.deaths=[deaths]
               
    def __repr__(self) :
        return self.name + ' ' + str(self.date[-1]) + ' ' + str(self.cases[-1]) + ' ' + str(self.deaths[-1]) + ' '
    
    # define dt as the time since the 5th case
    def deltatime(self) :
        for i0 in range(len(self.cases)) : 
            if self.cases[i0] >= 5 : break
        
        self.dt = [ (dd - self.date[i0]).days for dd in self.date ]
    
    # fit the number of cases vs time to 
    # 2 exponential functions with different
    # growth rates which meet at time = tbend
    def fitcases(self) :
        
        # uncertainty = Poisson errors
        self.dcases = [cc**0.5 for cc in self.cases]
        
        # fitting
        pars, cov = optz.curve_fit(expfunc, self.dt, self.cases, p0=[5, 3, 12, 5], sigma=self.dcases)
        
        # save the best-fit parameters to the County instance
        self.aa = pars[0]
        self.t0 = pars[1]
        self.tbend = pars[2]
        self.t1 = pars[3]
        
        # save best-fit parameters to the pandas dataset
        census.loc[census.area_name==self.name, 't0'] = float(self.t0)
        census.loc[census.area_name==self.name, 't1'] = float(self.t1)
        census.loc[census.area_name==self.name, 'tbend'] = float(self.tbend)

    # fit the number of deaths vs time to 
    # 2 exponential functions with different
    # growth rates which meet at time = tbend       
    def fitdeaths(self) :
        
        # remove zero points (usually at the start)
        self.nonzerodeaths = np.array(self.deaths)
        self.nonzerodeaths = self.nonzerodeaths[np.where(self.nonzerodeaths != 0)]
        
        # Poisson errors
        self.ddeaths = self.nonzerodeaths ** 0.5
        
        # Define new time axis that exclude zero points
        self.dt2 = np.array(self.dt)[np.where(self.nonzerodeaths != 0)]
        
        pars, cov = optz.curve_fit(expfunc, self.dt2, self.nonzerodeaths, p0=[0.1, 3, 12, 5], sigma=self.ddeaths)
        
        self.t0d = pars[1]
        
        census.loc[census.area_name==self.name, 't0d'] = float(self.t0d)
        
        


# create new dataset (counties) if it doesn't exist
try :
    
    ll=len(counties)
    
except NameError :

    counties = {}

    ii=0
    while ii < nyt.shape[0] :
    
    
        if nyt.loc[ii, 'county'] != 'Unknown' :
        
            if nyt.loc[ii, 'ctst'] not in counties.keys() : 
                
                counties[nyt.loc[ii, 'ctst']] = County( nyt.loc[ii, 'ctst'], nyt.loc[ii, 'date'], nyt.loc[ii, 'cases'], nyt.loc[ii, 'deaths'] )
                
            else :
                
                counties[nyt.loc[ii, 'ctst']].date.append( nyt.loc[ii, 'date'] )
                counties[nyt.loc[ii, 'ctst']].cases.append( nyt.loc[ii, 'cases'] )
                counties[nyt.loc[ii, 'ctst']].deaths.append( nyt.loc[ii, 'deaths'] )
                
        ii = ii + 1
   
    for name in counties : 
        counties[name].deltatime()     



def expfunc(tt, aa, t0, tbend, t1) :
        
    nexp = []
    
    for ttt in tt : 
        if ttt < tbend : nexp.append(aa * np.exp(ttt/t0) )
        else : nexp.append( aa * np.exp(tbend/t0) * np.exp((ttt-tbend)/t1) )

    return nexp



def casesvsdensity() :

    for name in counties : 
    
        county=counties[name]
    
        if len(county.dt) >= 10 :
            
            # fit the cases vs time for each county if we have 10 or more
            # data points
            county.fitcases()

    # plot of exponential constant vs density
    pl.figure(figsize=(12, 8))
    pl.loglog(census.density[census.t0.notna() & (census.t0 < 1e6)], 1/census.t0[census.t0.notna() & (census.t0 < 1e6)], '.', alpha=0.2)
    
    def powlaw(xx, aa, pp) :
        return aa * xx**pp
    
    pars, cov = optz.curve_fit(powlaw, census.density[census.t0.notna() & (census.t0 < 1e6)], 1/census.t0[census.t0.notna() & (census.t0 < 1e6)], p0=[0.1, 1])
    
    
    xx = np.linspace(1, 3e4, 2)
    pl.loglog(xx, powlaw(xx, pars[0], pars[1]))
    pl.xlabel('Population density (person/sq mile)')
    pl.ylabel('Exponential rate (per day)')
    
    pl.show()

    return pars, cov

def deathsvsdensity() :

    for name in counties : 
    
        county=counties[name]

        if sum([dd!=0 for dd in county.deaths]) >= 10 :
            
            county.fitdeaths()

    pl.figure(figsize=(12, 8))
    pl.loglog(census.density[census.t0d.notna() & (census.t0d < 1e3)], 1/census.t0d[census.t0d.notna() & (census.t0d < 1e3)], '.', alpha=0.2)
    
    def powlaw(xx, aa, pp) :
        return aa * xx**pp
    
    pars, cov = optz.curve_fit(powlaw, census.density[census.t0d.notna() & (census.t0d < 1e3)], 1/census.t0d[census.t0d.notna() & (census.t0d < 1e3)], p0=[0.1, 1])
    
    xx = np.linspace(1, 3e4, 2)
    pl.loglog(xx, powlaw(xx, pars[0], pars[1]))


def checkfit() :
    
    #pl.figure()
    """pl.ion()
    pl.show()"""
    
    
    ii=0
    for name in counties : 
    
        county=counties[name]


#        try :
#            t0 = county.t0 
#            
#            if ii % 2 == 0 : 
#                fig, axs = pl.subplots(1, 2, figsize=(12,8))
#                
#            # plot the data
#            axs[ii % 2].semilogy(county.dt, county.cases, 'o')
#            
#            # plot the best-fit model
#            xx = np.linspace(min(county.dt), max(county.dt), 20)
#            axs[ii % 2].plot(xx, expfunc(xx, county.aa, county.t0, county.tbend, county.t1))
#
#            axs[ii % 2].xlabel('Days after 5th case')
#            axs[ii % 2].ylabel('Number of cases')
#            axs[ii % 2].title(name)
#
#            if ii % 2 == 1: 
#
#                pl.show()
#                
#                _ = input('Enter "q" to exit loop, anything else to continue: ')
#                if _ == 'q' : break
#                #elif _ == '' : pass
#                #else : exec(_)
#            
#        except : pass


#        try :
#            t0 = county.t0 
#            
#            if ii % 2 == 0 : 
#                print('create figs')
#                
#            print('plot ', ii%2)
#
#            if ii % 2 == 1: 
#
#                print('show both figs')
#                
#                _ = input('Enter "q" to exit loop, anything else to continue: ')
#                if _ == 'q' : break
#                #elif _ == '' : pass
#                #else : exec(_)
#            
#        except : pass
#        ii += 1

        try :
            t0 = county.t0 
            
            if ii % 2 == 0 : 
                pl.figure(figsize=(12,8))
            pl.subplot(1,2, (ii % 2)+1 )
            
            # plot the data
            pl.semilogy(county.dt, county.cases, 'o')
            
            # plot the best-fit model
            xx = np.linspace(min(county.dt), max(county.dt), 20)
            pl.plot(xx, expfunc(xx, county.aa, county.t0, county.tbend, county.t1))

            pl.xlabel('Days after 5th case')
            pl.ylabel('Number of cases')
            pl.title(name)
            
            if ii % 2 == 1: 

                pl.show()
                
                _ = input('Enter "q" to exit loop, anything else to continue: ')
                if _ == 'q' : break
            
        except : pass
        ii += 1 
    
def plotall() : 
    
    pl.figure(figsize=(12, 8))

    for name in counties : 
        county=counties[name]
        pl.semilogy(county.dt, county.cases, 'k', alpha=0.2)

    pl.xlabel('Days after 5th case')
    pl.ylabel('Number of cases')
    pl.show()
        
def plotallscaled() : 
    
    #pl.figure()

    for name in counties : 
        county=counties[name]
        #pdb.set_trace()
        pl.semilogy(county.dt, 
                    np.array(county.cases)/
                    float(census.loc[census.area_name==name, 'pop']), 
                    'k', alpha=0.2)
        
    pl.xlabel('Days after 5th case')
    pl.ylabel('Number of cases / total population')
    pl.show()
    