#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 20:33:06 2018

@author: DewiGould

REFERENCES:
    Log-binning code used belongs to Imperial College London (supplied as part of Project)
"""

import networkx as nx
from random import randint, choice, sample, random,seed
import seaborn as sns #make plots look 'nicer'
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import time
from scipy.stats import linregress,ks_2samp
from scipy.interpolate import interp1d
from log_bin_CN_2016 import log_bin
import math
import itertools
sns.set(font_scale=2.5, style="ticks") #set styling preferences for matplotlib plots

class Network:
    """
    Main class 'Network'.
    
    Initialised with the following inputs:
        N = total number of nodes required in your network.
    Class contains several methods required throughout the project.
    
    Various different models (BA, random attachment, mixed attachment) can be implemented throughout
    the below methods by assigning either 'BA', 'RA' or 'MA' respectively to the value of the 'model'variable in the function inputs.
    
    Functions allow for tests of multiple values of m and/or N, for any of the modes described above, through
    the input parameters of the individual methods (see functions for individual descriptions)
    """
    
    def __init__(self,N):
        self.N = N              #total number of nodes in the network
        self.t = 0              #time
        self.my_list = []       #add node to this list every time you add an edge to network (for preferential attachment)
        self.G = nx.Graph()     #This is the graph instance - contains all edges and nodes through NetworkX module
        
    def empty(self): 
        """Empties all data fields to build a new graph.
        """
        self.__init__(self.N)
       
    def to_gephi(self):
        """
        Write networkx graph instance as .gexf file to be visualised in Gephi software package.
        """
        nx.write_gexf(self.G,"graph_visual.gexf")
        
    def maximum_degree(self):
        """Returns the maximum degree of the current network contained within self.G object
        """
        degrees_raw = nx.degree(self.G).values() #list of degrees of each node
        return max(degrees_raw)
        
    def degree_distribution(self):
        """ Returns a tuple containing: 
            [dictionary of observed degrees and their frequency]
            [list of all raw degrees as observed on the network]
        """
        degrees_raw = nx.degree(self.G).values() #list of degrees of each node
        degrees = Counter(degrees_raw)            
        return degrees,degrees_raw #return raw degrees for log-binning
        

    def p_inf_predicted(self,k,m,model):
        """
        Function returns the predicted value of the time stationary degree probability p_infinity(k)
        Inputs:
            k = degree
            m = number of nodes added per time step
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
        Outputs:
            Value of predicted degree probability
        """
        if model == "BA":
            k = float(k)
            A = 2*m*(m+1)
            return float(A)/float(k*(k+1)*(k+2))
        if model == "RA":
            frac = (float(m)/float(m+1))**(k-m)
            return frac/float(m+1)
        if model == "MA":
            m = float(m)
            k = float(k)
            return (12*m*((3*m)+1)*((3*m)+2)*((3*m)+3))/((k+(2*m))*(k+(2*m)+1)*(k+(2*m)+2)*(k+(2*m)+3)*(k+(2*m)*4))
      
    def k1_predicted(self,m,N,model):
        """
        Function returns the predicted value for the maximum degree: k1
        Inputs:
            m = number of nodes added per time step
            N = total number of nodes in the model being tested
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
        Outputs:
            Value of predicted maximum degree, k1
        """
        if model == "BA":
            return (-1 + np.sqrt(1+(4*N*m*(m+1))))/2.0
        if model == "RA":
            return m - ((math.log(N))/((math.log(m))-(math.log(m+1))))

    def algo_check(self,model):
        """
        Function to perform algorithm checks for the probability implementations for each of the three models
        Inputs:
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
        Outputs:
            Graph of node choice frequency versus degree
        """
        if model == "BA":
            results = []
            seed(0)                                         #set random seed
            s = 0
            for p in range(10):
                vals = []
                s +=1 
                seed(s)                                     #update random seed
                degrees = np.linspace(1,100,100)
                my_list = []
                for d in degrees:                           #create uniform degree distribution
                    for k in range(int(d)):
                        my_list.append(d)            
                degs = Counter(my_list)
                for j in range(10000):                    #Run node choice algorithm 10^6 times and store results
                    choose = choice(my_list)                #BA algorithm
                    vals.append(degs[choose])               #Find the degree of the node chosen
                results.append(vals)
            plt.figure()
            for p in results:
                res = Counter(p)
                plt.plot(res.keys(),res.values())
            plt.xlabel("Vertex Degree")
            plt.ylabel("Node Choice Frequency")
            
        if model == "RA":                                   #Same process by analogy to BA case
            results = []
            seed(0)
            s = 0
            for p in range(10):
                vals = []
                s +=1
                seed(s)
                nodes = np.linspace(1,100,100)
                my_list = []
                for d in nodes:
                    for k in range(int(d)):
                        my_list.append(d)            
                degs = Counter(my_list)
                for j in range(10000):
                    choose = choice(nodes)                  #RA algorithm
                    vals.append(degs[choose])
                results.append(vals)
            plt.figure()
            for p in results:
                res = Counter(p)
                plt.plot(res.keys(),res.values())
            plt.xlabel("Vertex Degree")
            plt.ylabel("Node Choice Frequency")  
        
        if model == "MA":                                   #Same process by analogy to BA case
            plt.figure()
            results = []
            seed(0)
            s = 0
            for p in range(10):
                vals = []
                s +=1
                seed(s)
                degrees = np.linspace(1,100,100)
                my_list = []
                for d in degrees:
                    for k in range(int(d)):
                        my_list.append(d)            
                degs = Counter(my_list)
                for j in range(10000):
                    choose = choice(my_list)                #MA algorithm
                    vals.append(degs[choose])
                results.append(vals)
            plt.figure()
            res = Counter(results[0])
            plt.plot(res.keys(),res.values(),color="green",label="q=1")
            for p in results[1:]:
                res = Counter(p)
                plt.plot(res.keys(),res.values(),color="green")
            results = []
            seed(0)
            s = 0
            for p in range(10):
                vals = []
                s +=1
                seed(s)
                nodes = np.linspace(1,100,100)
                my_list = []
                for d in nodes:
                    for k in range(int(d)):
                        my_list.append(d)            
                degs = Counter(my_list)
                for j in range(10000):
                    choose = choice(nodes)
                    vals.append(degs[choose])
                results.append(vals)
            res = Counter(results[0])
            plt.plot(res.keys(),res.values(),color="blue",label="q=0")
            for p in results[1:]:
                res = Counter(p)
                plt.plot(res.keys(),res.values(),color="blue")
            plt.xlabel("Vertex Degree")
            plt.ylabel("Node Choice Frequency") 
            plt.legend(fontsize=30)


            
    def varying_m(self,plot,checking,m_range,tests,model,testing):
        """
        Function to perform first set of tasks in the Project Notes - examine behaviour for multiple 'm' values.
        Inputs:
            plot = Boolean: True = see plots of raw degree distribution
            checking = Boolean: True = see plots of log-binned degree distribution
            m_range = list of m values to be tested.
            tests = Boolean: True - perform statistical tests on the data and theoretical prediction, output is graph.
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
            testing = Boolean: set to True to run SMALLER number of repetitions to quickly check operations
        Ouputs:
            For all three boolean inputs, the outputs are graphs used in the final project report.   
        
        Statistical Testing:
            A Kolmogorov-Smirnov test was performed to compare the predicted degree distributions with the obtained data.
        To examine the performance of the models, the KS test was performed over an incrementally increasing range of data points
        I.e. starting with a very small window (k=m,k=m+1), and expanding this up to the entire data set.        
        """
        results = []
        for m in m_range:
            start = time.time()
            if m == 1:
                n = int(0)
            else:
                n= int(float(m)**(1.0/3.0))
            repetitions = 100*(3**(4-n)) #Repetitions are maximum 81000, down to 100 at m=81
            if testing:
                repetitions = 100
            print repetitions
            values = []
            count = 0
            for j in range(repetitions): #build graphs depending on model type.
                if model == "BA":
                    model_test = BAmodel(self.N,2*m,m,"new")
                elif model == "RA":
                    model_test = RAmodel(self.N,2*m,m)
                elif model == "MA":
                    model_test = MAmodel(self.N,2*m,m,0.5)
                kdist,degrees_raw =  model_test.degree_distribution() #get degree distribution
                values.extend(degrees_raw)
                self.empty()
                if count in [1,2,3,4,5,10,20,100,200,300,400,500,1000,2000,3000,4000,5000,6000,10000]:
                    print count
                    
                count+=1
            kdist = Counter(values) #collect data over all repetitions
            results.append((m,kdist,values,int(repetitions*self.N))) #collect full data for this value of m
            print "Finished for m = ", m
            print "time taken: ", time.time()-start
            
        for_plot, for_check,for_tests = [[],[],[]] #empty lists to be filled with data to be plotted
        tests_2 = []
        for values in results:
            if plot:
                N = len(values[2])
                degrees_filtered = [i for i in values[1].keys() if i>=values[0]] #remove k<m
                p_k = [float(values[1][i])/float(N) for i in degrees_filtered] #find theoretical values for these k values
                for_plot.append((degrees_filtered,p_k,values[0]))
            if checking:
                degrees = [i for i in values[2] if i>=values[0]] #filter k<m
                centres, counts = log_bin(degrees,min(degrees),1.0,1.1) #log bin the degree distribution data.
                for_check.append((centres,counts,values[0]))
            if tests:
                degrees = [i for i in values[2] if i>=values[0]] #filter k<m
                results_two = []
                c = Counter(degrees)
                degrees_obs = c.keys()
                counts = c.values()
                for i in np.arange(values[0],max(degrees_obs)+1):  #define variable window to perform tests over.
                    true = counts[degrees_obs.index(values[0]):i]
                    degs = degrees_obs[degrees_obs.index(values[0]):i]
                    if len(true) != 0:
                        pred = [self.p_inf_predicted(l,values[0],model) for l in degs]
                        pred = [l*(float(sum(true))/(float(sum(pred)))) for l in pred] #get predicted counts with SAME sum as obtained results
                        if len(true) == len(pred): #to account for slight over production
                            results_two.append((i,ks_2samp(true,pred)[1])) #perform KS tests
                for_tests.append((values[0],results_two))
            if tests: #Do same KS tests but on CUMULATIVE DISTRIBUTION  to deal with fat tail CDF Data
                degrees = [i for i in values[2] if i>=values[0]] #filter k<m
                ress_2 = []
                centres = Counter(degrees).keys()
                counts = np.cumsum(Counter(degrees).values())
                for i in centres:
                    if i != min(centres):
                        true = counts[centres.index(min(centres)):centres.index(i)]
                        true = [int(values[3]*a) for a in true] #make them a count not a probability.
                        degs = centres[centres.index(min(centres)):centres.index(i)]
                        if len(true) !=0:
                            pred = [self.p_inf_predicted(l,values[0],model) for l in degs]
                            pred = np.cumsum(pred).tolist()
                            pred = [l*(float(sum(true))/(float(sum(pred)))) for l in pred]
                            if len(true) == len(pred):
                                ress_2.append((i,ks_2samp(true,pred)[1]))
                tests_2.append((values[0],ress_2))                              
        if plot:
            plt.figure()
            for j in for_plot:
                plt.scatter(j[0],j[1],marker='o',edgecolors = 'black')
                plt.plot(np.linspace(1,max(j[0]),100),[self.p_inf_predicted(i,j[2],model) for i in np.linspace(1,max(j[0]),100)],'--',label=r'$p_{{\infty}}(k), m={0}$'.format(j[2]))
            plt.xlabel("k")
            plt.ylabel("p(k)")
            plt.legend(fontsize=30)
            plt.xscale("log")
            plt.yscale("log")
            plt.show()  
        if checking:
            plt.figure()
            for j in for_check:
                 plt.scatter(j[0],j[1],marker='o',edgecolors = 'black')
                 plt.plot(np.linspace(j[2],max(j[0]),100),[self.p_inf_predicted(i,j[2],model) for i in np.linspace(j[2],max(j[0]),100)],'--',label=r'$p_{{\infty}}(k), m={0}$'.format(j[2]))
            plt.xlabel("k")
            plt.ylabel(r"$\tilde{p}(k)$")
            plt.legend(fontsize=30)
            plt.xscale("log")
            plt.yscale("log")
            plt.show()
            
            plt.figure()
            for j in for_plot:
                cdf = np.cumsum(j[1])
                plt.scatter(j[0],cdf,marker='o',edgecolors = 'black')
                plt.plot(j[0],np.cumsum([self.p_inf_predicted(i,j[2],model) for i in j[0]]),'--',label=r'$p_{{\infty}}(k), m={0}$'.format(j[2]))
            plt.xlabel("k")
            plt.ylabel(r"$\tilde{p}(k)$")
            plt.legend(fontsize=30)
            plt.xscale("log")
            plt.yscale("log")
            plt.show()
                
        if tests:
            plt.figure()
            for res in for_tests:
                plt.plot([j[0] for j in res[1]],[j[1] for j in res[1]],label="m=%s"%res[0])
            plt.xlabel(r"$k_{max}$")
            plt.ylabel("K-S Test, p-val")
            plt.legend(fontsize=30)
            plt.show()
            
            plt.figure()
            for res in tests_2:
                plt.plot([j[0] for j in res[1]],[j[1] for j in res[1]],label="m=%s"%res[0])
            plt.xlabel(r"$k_{max}$")
            plt.ylabel("K-S Test, p-val")
            plt.legend(fontsize=30)
            plt.show()
                        

    def max_k(self,N_range,m_range,model,testing):
        """
        Function to perform analysis on the maximum degree of the model.
        Inputs:
            N_range = list of N values to be tested (N = total number of nodes in network)
            m_range = list of m values to be tested (m = number of nodes added per time step)
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
            testing = Boolean: set to True to run SMALLER number of repetitions to quickly check operations

        Output:
            Graph of theory vs. data measurements for maximum degree of network.
            Any other graphs used in report during this section.
            
        Function runs simulations, over a number of repetitions, to find the average numerical maximum degree of the network.
        This value is compared against the theoretical prediction of the maximum degree for that given m and N.
        
        A comparison of the theory and data is made graphically by plotting the two with associated errors.
        """
        results = []
        r2 = []
        for m in m_range:
            inter = []
            int2 = []
            for n in N_range:
                start = time.time()
                repetitions = 100000000/n #1million nodes in total (10^7-n)
                if testing:
                    repetitions = 100
                print repetitions
                values = []
                count = 0
                for j in range(repetitions): #create graphs according to user-defined algorithm of choice.
                    if model == "BA":
                        model_test = BAmodel(int(n),2*m,m,"new") 
                    elif model == "RA":
                        model_test = RAmodel(int(n),2*m,m)
                    elif model == "MA":
                        model_test = MAmodel(int(n),2*m,m,0.5)
                    kdist, degrees = model_test.degree_distribution()
                    values.append(max(degrees)) #measure maximum degree
                    if count in [1,2,3,4,5,10,20,30,40,50,100,1000,2000,3000,5000,6000]:
                        print count
                    count +=1
                print "Done in time: ", time.time()-start
                inter.append((n,np.mean(values),np.std(values)))  #take simple average and find standard deviation
                int2.append((n,values,repetitions))
            results.append((m,inter))
            r2.append((m,int2))
        
        for values in r2: #Plot distribution of measurements of max_k
            m, points = values
            plt.figure()
            palette = itertools.cycle(sns.color_palette()) #set colours of plot
            for p in points:               
                a = Counter(p[1])
                y = [float(i)/float(p[2]) for i in a.values()]
                plt.plot(a.keys(),y,zorder=1,marker='o',color=next(palette), label="N=%s"%p[0])
            plt.xlabel(r"$k_{1}$")
            plt.ylabel(r"$p(k_{1})$")
            plt.legend(fontsize=30)
        for values in r2: #Plot distribution of measurements of max_k
            m, points = values
            plt.figure()
            palette = itertools.cycle(sns.color_palette()) #set colours of plot
            for p in points:               
                a = Counter(p[1])
                y = [float(i)/float(p[2]) for i in a.values()]
                plt.scatter(a.keys(),y,zorder=1,color=next(palette), label="N=%s"%p[0],edgecolors = 'black')
            plt.xlabel(r"$k_{1}$")
            plt.ylabel(r"$p(k_{1})$")
            plt.legend(fontsize=30)          
        new_results = []
        for values in r2: #perform data collapse using hypothesis of MB distribution
            m, points = values
            plt.figure()
            palette = itertools.cycle(sns.color_palette()) #set colours of plot
            vals = []
            for p in points:
                a = Counter(p[1])
                x1 = a.keys()
                y = a.values()
                vp = float(x1[y.index(max(y))]) #most probably value
                y = [float(i)/float(p[2]) for i in y]
                y = [float(i)*float(vp) for i in y]
                x1 = [float(i)/vp for i in x1]
                pairs = [(i,y[x1.index(i)]) for i in x1]
                pairs = sorted(pairs, key=lambda a: a[0])
                xplot = [l[0] for l in pairs]
                yplot = [l[1] for l in pairs]
                f = interp1d(xplot,yplot,"cubic")
                xplot = np.linspace(min(xplot),max(xplot),50)
                yplot = f(xplot)
                plt.plot(xplot,yplot,marker='o', color=next(palette), label="N=%s"%p[0])  
                mean = (2.0/np.pi)*vp
                err = (2.0/np.pi)*np.std(p[1])
                err2 = err
                if err > vp:
                    err2 = 0.999*vp
                vals.append((p[0],mean,err,err2))
            new_results.append((m,vals))
            plt.xlabel(r"$\frac{k_{1}}{k_{p}}$")
            plt.ylabel(r"$k_{p}p(k_{1})$")
            plt.legend(fontsize=30)
             
        for values in r2: #perform data collapse using hypothesis of MB distribution
            m, points = values
            plt.figure()
            palette = itertools.cycle(sns.color_palette()) #set colours of plot
            for p in points:
                a = Counter(p[1])
                x1 = a.keys()
                y = a.values()
                vp = float(x1[y.index(max(y))]) #most probably value
                y = [float(i)/float(p[2]) for i in y]
                y = [float(i)*float(vp) for i in y]
                x1 = [float(i)/vp for i in x1]
                pairs = [(i,y[x1.index(i)]) for i in x1]
                pairs = sorted(pairs, key=lambda a: a[0])
                xplot = [l[0] for l in pairs]
                yplot = [l[1] for l in pairs]
                f = interp1d(xplot,yplot,"cubic")
                xplot = np.linspace(min(xplot),max(xplot),50)
                yplot = f(xplot)
                plt.scatter(xplot,yplot,marker='o', color=next(palette), label="N=%s"%p[0],edgecolors='black')  
            plt.xlabel(r"$\frac{k_{1}}{k_{p}}$")
            plt.ylabel(r"$k_{p}p(k_{1})$")
            plt.legend(fontsize=30) 
       
        plt.figure()
        palette = itertools.cycle(sns.color_palette()) #set colours of plot
        for values in results:
            m, points = values
            color = next(palette)
            if model == "BA" or model == "MA": #find slope of data line using linear regression
                slope, intercept, r_value, p_value, std_err = linregress([np.log(i[0]) for i in points],[np.log(i[1]) for i in points])
                print "Measured slope: ", slope, "with error ", std_err            
                plt.errorbar([i[0] for i in points],[i[1] for i in points],yerr=[i[2] for i in points],color=color,marker='o',capsize=5,label = r"$m={0}, slope = {1:.3f} \pm{2:.3f}$".format(m,slope,std_err))
            else:
                plt.errorbar([i[0] for i in points],[i[1] for i in points],yerr=[i[2] for i in points],color=color,marker='o',capsize=5,label = "m={0}".format(m))
            x_range = np.linspace(np.log10(min([i[0] for i in points])),np.log10(max([i[0] for i in points])),100)
            x_range = [10**x for x in x_range]
            if model != "MA":
                plt.plot(x_range,[self.k1_predicted(m,i,model) for i in x_range],color=color,linestyle='dashed',label=r"$ Theoretical Prediction: \langle k_{1} \rangle$")
        plt.xlabel("N")
        plt.ylabel(r"$k_{1}$")
        plt.legend(fontsize=30)
        plt.xscale("log")
        plt.yscale("log")
        plt.show() 
        
        plt.figure() #with updated errors and means
        palette = itertools.cycle(sns.color_palette()) #set colours of plot
        for values in new_results:
            m, points = values
            color = next(palette)
            if model == "BA" or model == "MA": #find slope of data line using linear regression
                slope, intercept, r_value, p_value, std_err = linregress([np.log(i[0]) for i in points],[np.log(i[1]) for i in points])
                print "Measured slope: ", slope, "with error ", std_err            
                plt.errorbar([i[0] for i in points],[i[1] for i in points],yerr=[i[2] for i in points],color=color,marker='o',capsize=5,label = r"$m={0}, slope = {1:.3f} \pm{2:.3f}$".format(m,slope,std_err))
            else:
                plt.errorbar([i[0] for i in points],[i[1] for i in points],yerr=[i[2] for i in points],color=color,marker='o',capsize=5,label = "m={0}".format(m))
            x_range = np.linspace(np.log10(min([i[0] for i in points])),np.log10(max([i[0] for i in points])),100)
            x_range = [10**x for x in x_range]
            if model != "MA":
                plt.plot(x_range,[self.k1_predicted(m,i,model) for i in x_range],color=color,linestyle='dashed',label=r"$ Theoretical Prediction: \langle k_{1} \rangle$")
        plt.xlabel("N")
        plt.ylabel(r"$k_{1}$")
        plt.legend(fontsize=30)
        plt.xscale("log")
        plt.yscale('log', nonposy='clip')
        plt.show()         
        
    def varying_N(self,N_range,m,model,plot,collapse,deviation,testing):
        """
        Function to perform tasks associated with the 'Varying N' section of the Project.
        Inputs:
            N_range = list of N values to be tested (N = number of nodes in network)
            m = number of nodes added per time step.
            model = 'BA', 'RA', 'MA': for either Barabasi-Albert, Random Attachment or Mixed Attachment.
            plot = Boolean: True - to plot obtained data (log-binned)
            collapse = Boolean: True - to perform data collapse and observe resulting graph
            deviation = Boolean: True - analyse the deviation from theory of the data.
            testing = Boolean: set to True to run SMALLER number of repetitions to quickly check operations

        Outputs:
            Outputs of the various methods are all graphs included in this section of the report for visual and statistical analysis.
        
        Data Collapse performed according to procedure laid out in Project Report.
        
        """
        results = []
        for n in N_range:
            start =time.time()
            repetitions= int(100000000/float(n)) #repetitions to be averaged over
            if testing:
                repetitions = 100
            print repetitions
            values = []
            count = 0
            for j in range(repetitions): #set-up the networks
                if model == "BA":
                    model_test = BAmodel(int(n),2*m,m,"new")
                elif model == "RA":
                    model_test = RAmodel(int(n),2*m,m)
                elif model == "MA":
                    model_test = MAmodel(int(n),2*m,m,0.5)
                kdist,degrees_raw =  model_test.degree_distribution()
                values.extend(degrees_raw)
                self.empty()
                if count in [1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000,10000,50000,75000,100000]:
                    print count
                count +=1
            print "Time taken: ", time.time()-start
            kdist = Counter(values)
            results.append((int(n),kdist,values))
                    
        to_plot, to_collapse, to_deviate = [[],[],[]] #value to be filled for plotting
        for values in results:
            degrees = values[2]
            degrees = [i for i in values[2] if i>=m]
            if min(degrees) != 0:
                centres, counts = log_bin(degrees,min(degrees),1.0,1.1) #log-bin the degree distribution
            else:
                centres, counts = log_bin(degrees,m,1.0,1.1) #log-bin the degree distribution
            counts = counts.tolist()
            if plot:
                to_plot.append((centres,counts,values[0],degrees))
            if collapse:
                if model == "MA":
                    print "Can't do this, no analytic expression found for k1"
                else:
                    p_inf = [self.p_inf_predicted(l,m,model) for l in centres] #calculated expected probabilities
                    k1 = [self.k1_predicted(m,values[0],model) for l in centres]
                    collapsed = [float(l)/float(p_inf[counts.index(l)]) for l in counts] #collapsed
                    kcoll = [float(l)/float(k1[centres.index(l)]) for l in centres]      
                    to_collapse.append((kcoll,collapsed,values[0]))
            if deviation:
                if model == "MA":
                    print "Can't do this - don't have an analytic expression for k1"
                else:
                    if collapse != True:
                        p_inf = [self.p_inf_predicted(l,m,model) for l in centres]
                        collapsed = [float(l)/float(p_inf[counts.index(l)]) for l in counts]
                    to_deviate.append((centres, collapsed,values[0]))      
        if plot:
            plt.figure()
            for j in to_plot:
                plt.scatter(j[0],j[1],marker='o',edgecolors = 'black',label='N=%s'%j[2])
            plt.plot(np.linspace(1,max(j[3]),100),[self.p_inf_predicted(i,m,model) for i in np.linspace(1,max(j[3]),100)],'--',label=r'$p_{{\infty}}(k), m={0}$'.format(m))
            plt.xlabel("k")
            plt.ylabel(r"$\tilde{p}(k)$")
            plt.legend(fontsize=30)
            plt.xscale("log")
            plt.yscale("log")
            plt.show() 
            
        if collapse and model != "MA":
            if len(to_collapse) != 0:
                plt.figure()
                for j in to_collapse:
                    plt.scatter(j[0],j[1],marker='o',edgecolors='black',label="N=%s"%j[2])
                plt.xlabel(r"$\frac{k}{k_{1}}$")
                plt.ylabel(r"$\frac{p_{observed}}{p_{\infty}}$")
                plt.xscale("log")
                plt.yscale("log")
                plt.legend(fontsize=30)
                plt.show()
                
        if deviation and model != "MA":
            plt.figure()
            for j in to_deviate:
                plt.scatter(j[0],j[1],marker='o',edgecolors='black',label="N=%s"%j[2])
            plt.plot(np.linspace(1,max(j[0]),100),[1]*100,'--',label=r"$\frac{p_{observed}}{p_{\infty}} = 1$")
            plt.xlabel("k")
            plt.ylabel(r"$\frac{p_{observed}}{p_{\infty}}$")
            plt.xscale("log")
            plt.yscale("log")
            plt.legend(fontsize=30)
            plt.show()   
                          
            
class BAmodel(Network):
    """Type of Network: Barabasi-Albert Model.
    Nodes added using preferntial attachment algorithm
    
    Inputs:
        N = total number of nodes
        N_initial = number of nodes in the initial configuration
        m = number of nodes added per time step.
        method = "one", "two", "three", "new" - define the different algorithms to generate preferential attachment
        
    Outputs:
        self.G will be fully generated BA graph according to the above inputs.
    """
    
    def __init__(self,N,N_initial,m,method):
        Network.__init__(self,N)
        self.m = m
        self.N_initial = N_initial #initial number of nodes at t=0 (must be >m)
        self.method = "None"
        #Initialise graph at t=0
        for i in range(N_initial):
            self.G.add_node(i)         #Add nodes        
        
        #Choose to create a random Graph with E and N fixed - edges added completely randomly.
        N_edges_initial = int(m)                 #Choose E(0) = 10m
        for i in range(N_edges_initial):
            current_num_edges = self.G.number_of_edges()
            while self.G.number_of_edges() == current_num_edges:    #make sure to add a unique edge
                vertex1 = randint(0,self.N_initial) #randomly choose a vertex to add edge to
                vertex2 =  vertex1
                while vertex2 == vertex1:                    #avoid self-edges
                    vertex2 = randint(0,self.N_initial)     #randomly choose an edge to create
                self.G.add_edge(vertex1,vertex2)
            self.my_list.append(vertex1)
            self.my_list.append(vertex2)
        self.build(method) #Build the graph via preferential attachment
     
    def build(self,method): #add nodes with preferential attachment
        self.method = method
        T = self.N - self.N_initial #period
        count = 0
        for i in range(self.N_initial,T+self.N_initial):
            self.t+=1 #increment time
            self.G.add_node(i) #Add vertex
            if method == "one": #Third method in the notes
                for j in range(self.m): #Add m edges
                    current_num_edges = self.G.number_of_edges()
                    while self.G.number_of_edges() == current_num_edges: #make sure the edge you add is a new
                        vertex2 = i
                        while vertex2 == i: #avoid self-edges
                            vertex2 = choice(self.my_list) #choose node propto to its degree by randomly choosing from this list
                        count +=1
                        self.G.add_edge(i, vertex2)
                    self.my_list.append(i)
                    self.my_list.append(vertex2)
            if method == "new": #quicker version of Method "One" - without checks for unique edges and self-loops (faster)
                for j in range(self.m):
                    v2 = choice(self.my_list) #choose node propto its degree.
                    self.G.add_edge(i,v2)
                    self.my_list.append(i)
                    self.my_list.append(v2)                                 
            if method == "two": #First method in the notes
                for j in range(self.m):
                    current_num_edges = self.G.number_of_edges()
                    while self.G.number_of_edges() == current_num_edges: #avoid adding edge that already exists
                        vertex2 = i 
                        while vertex2 == i: #avoid self-edges
                            edge = sample(self.G.edges(), 1)[0] #choose edge at random
                            vertex2 =  edge[randint(0,1)] #choose end of this edge at random
                        count +=1
                        self.G.add_edge(i,vertex2) 

class RAmodel(Network):
    """Type of Network: Random Attachment.
    Nodes added via random attachment algorithm
    
    Inputs:
        N = number of nodes
        N_initial = number of nodes in initial configuration
        m = number of nodes added per time step.
        
    Outputs:
        fully generated network contained within self.G
    """
    
    def __init__(self,N,N_initial,m):
        Network.__init__(self,N)
        self.m = m
        self.N_initial = N_initial #initial number of nodes at t=0 (must be >m)
        self.method = "None"
        #Create initial graph
        for i in range(N_initial):
            self.G.add_node(i) #add nodes, choose completely disconnected graph
        self.build_random()    #Build the graph, i.e. add the nodes via random attachment up until N nodes reached
        
    def build_random(self):
        T = self.N - self.N_initial #Period
        for i in range(self.N_initial,T+self.N_initial):
            self.t +=1 #increment time
            self.G.add_node(i) #add new node
            for j in range(self.m):
                v2 = randint(0,i-1) #randomly attach m new edges within this loop
                self.G.add_edge(i,v2)
   
class MAmodel(Network):
    """Type of Network: Mixed Attachment
    Nodes attached via mixed attachment algorithm.
    
    Inputs:
        N = number of nodes
        N_initial = number of nodes in initial configuration
        m = number of nodes added per time step.
        
    Ouputs:
        Fully generated network contained within self.G
    """
    
    def __init__(self,N,N_initial,m,q):
        Network.__init__(self,N)
        self.m = m
        self.N_initial = N_initial #initial number of nodes at t=0 (must be >m)
        self.q = q
        #Create initial graph
        for i in range(N_initial):
            self.G.add_node(i) #add nodes
        #Choose a Random Graph; E(0) = m10
        N_edges_initial = int(m) 
        for i in range(N_edges_initial):
            current_num_edges = self.G.number_of_edges()
            while self.G.number_of_edges() == current_num_edges: #make sure edges added are unique
                vertex1 = randint(0,self.N_initial)
                vertex2 =  vertex1
                while vertex2 == vertex1: #avoid self-edges
                    vertex2 = randint(0,self.N_initial)
                self.G.add_edge(vertex1,vertex2)
            self.my_list.append(vertex1)
            self.my_list.append(vertex2)
        self.build_mixed()        #Build the graph, i.e. add the nodes via mixed attachment up until N nodes reached
    
    def build_mixed(self):
        T = self.N - self.N_initial
        for i in range(self.N_initial,T+self.N_initial):
            self.t +=1
            self.G.add_node(i)
            for j in range(self.m):
                r = random()
                if r <= self.q: #preferential attachment
                    v2 = choice(self.my_list) #choose node propto its degree.
                else: #random attachment
                    v2 = randint(0,i-1) #choose node randomly from list of all available nodes.
                self.G.add_edge(i,v2)
                self.my_list.append(i)
                self.my_list.append(v2)