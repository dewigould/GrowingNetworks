# GrowingNetworks

Codes submitted:
NetworksCode.py

#—————————————————————————#

- Networks.py

Code contains all necessary scripts to obtain figures and numerical results used in report. 
The below snippets can be run to obtain sample results of everything used in report.


***In all cases below, if the parameter “BA” appears - this can simply by swapped for “RA” or “MA”
to test the same thing for the different models (Barabasi and Albert (BA), Random Attachment (RA) and
mixed attachment (MA)

Algorithm Testing (generate graphs used to check algorithms were working properly):

import NetworksCode as n
model = n.Network(10000) #10000 is number of nodes in graph to be built
model.algo_check(“BA”) 

***Testing parameters have been modified inline to make this run quickly to just give a rough picture of what was obtained.


	•	Varying m, fixed N (generate graphs used to observe graph attributes for varying m, but fixed N values)
import NetworksCode as n
model = n.Network(10000) #10000 is number of nodes in graph to be built
model.varying_m(True,True,[1,3,9],True,”BA”,True)

List [1,3,9] is list of mranges to be tested over
The three Boolean values can be modified to print/ not print certain figures (all set True here to see everything)
Fourth Boolean value is set to True to run a smaller number of repetitions so you can visualise some basic results quickly.


	•	Maximum k (generate graphs used to observe distribution of k_1 values and how it changes with N)
import NetworksCode as n
model = n.Network(10000) #10000 is number of nodes in graph to be built
model.max_k([100,1000,10000],[3],”BA”,True)

List [100,1000,10000] is list of different N to be tested
3 is the value of fixed m to be used
Fourth Boolean value is to set number of repetitions lower - for quick visualisation purposes for code checking.


	•	Varying N, fixed m (generate graphs used to observe graph attributes for varying N, but fixed m values)
import NetworksCode as n
model = n.Network(10000) #10000 is number of nodes in graph to be built
model.varying_N([100,1000,10000],3,”BA”,True,True,True,True)

list [100,1000,10000] is list of N values to be tested
3 is value of fixed m to be tested
Boolean parameters switch on and off the printing of values
Fourth Boolean value is to set number of repetitions lower - for quick visualisation purposes for code checking.

