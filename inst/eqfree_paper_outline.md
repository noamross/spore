# Outline For Equation-Free Optimization Paper

## Introduction

-   Control of disease as an economic problem: minimizing disease impacts in the face of limited resources for control.  [@Brandeau2003] Wildlife [@Horan2005; Bicknell1999; Fenichel2010]   (TODO See additional citations in OPENQUALS for variety of methods applied to invasives)
-   Importance of individual-based processes (network contact structure, host and pathogen variation) in disease dynamics.
-   IBMs work where population-level models do not capture appropriate dynamics or answer questions. [@Grimm2005]
### Other approaches

-  The general approach to optimization of IBMs has been to apply control techniques
   to reduced-form models.
   - Dynamic optimal control described in @Clark1990, using @Pontryagin1967
   
-  Applying optimal control techniques to a simplified mean-field model [@Federico2013]
   can work for homogenous models, but yields non-optimal results under heterogeneity.
-  Simplified models can be adjusted via non-mechanistic modifiers to closer match IBM
behavior not captured under mean-field conditions.  These adjusted models can yield
better results than mean-field models. However, the applicability of these adjusted models
is limited the the range over which they were parameterized, and model parameters
are difficult to interprety [Oremland 2015].
-  If reduced form models have low enough dimensionality, discretization can allow for
   approximate optimal control through  for heuristic methods that search the entire space of possible controls.[@Oremland2015].
-  In addition to optimal control, IBMs can be useful in comparison of a small set of management
scenarios [LOTS OF EXAMPLES].

### Introduction to Equation-Free Modeling

-   IBM dynamics can be difficult to capture as ODEs [@Brown2004]

-   Equation-free (EF), or multi-scale modeling [@Armaou2004, but see @Kevrekidis2009a for review]
is a method for capturing and analyzing population-scale dynamics of individual-based models while
bypassing the deriving population-scale equations, which exist conceptually but often
can not be derived.
-   Originally developed for applications in physical chemistry, it has been applied to biological  [@Erban2006], ecological [@Raghib2010a]
and epidemiological systems [@Cisternas2004; @Williams2015].
    
-  The EF approach allows computationally efficient simulations of the *expected value* of
   population-level summaries of an individual based model.
   
### Outline/Preview

-   Describe an optimization problem on host-parasite model system
-   Describe an EF-based method of dynamic optimization
-   Solve the problem with the EF methood
-   Demonstrate that EF-method yields better results than optimizing on the population-level model.

## Methods

### Model system

-   The IBM is an individual-based representation of the @Anderson1978 model of
    host-parisite dynamics.
    
    -   Continuous time, discrete population.
    -   Full distribution of infections represented
    -   Demographic stochasticity via Gillespie process
    -   External source of disease particles entering the population

-   ODE representation of the sytem makes assumptions of
    -   large population size (no stochasticity)
    -   constant shape of distribution of infections amongst individuals
    
-   ODE representation useful for equilibrium solutions, but not transients [@Adler1992; @Ross2015]
-   ODE system is reduced in the following ways:
    -   It is an approximation of an infinite system of differential equations
    -   It lacks demographic stochasticity
    
-   This IBM is relatively simple.  It lacks several features that are trivial to incorporate into an IBM but challenging to represent
in closed-form solutions, such as individual variation in traits or heterogenity in the host contact structure.  

-   Economic control problem:
    -   Control mechanism is reduction of rate of new disease particles entering the sytem.
    -   Maximize profit,  measured as time-integrated ecosystem service value of hosts, minus cost of control mechanism
    

-   [EQUATIONS]

### EF Control Approach

-   @Pontryagin1967's maximum principle.  

-   Use EF simulation to simulate dynamics of system as well as estimate derivatives via finite differencing
-   Derive adjoint equations as functions of numerical estimates of host and and pathogen dynamics
-   Simulate optimal path by solving for Hamiltonian maximum at each time step

-  Lifting and restriction functions translate between full distribution of infections and summary statistics
   of distribution
   -   Use Conway-Maxwell-Poisson to allow for both under- and over-dispersed distributions with changing
       aggregation.

## Results

-   IBM and ODE versions of the model have different behavior
-   EF version of IBM correctly produces results similar to overall IBM:
    -   Plots: Dynamics of ODE, IBM (average of many runs), and EF versions of model
        under conditions of similar behavior (large population, equilibrium), and 
        different behavior (small population, transients)
        
-   IBM and ODE have different responses to the the same controls.
    -   Characteristics of the control path as solved for the ODE
        -   Plot: ODE dynamics under control
        -   Plot: IBM dynamics (average and several example runs) under ODE-derived control
        
-   EF optimization yields better results for the IBM than ODE-derived control
    -   Plot: IBM dynamics under EF optimization
    -   Plot: Compare net present value of
        -   IBM and ODE models without control
        -   IBM and ODE models with ODE-derived control
        -   IBM under EF-derived control.

## Discussion

### Implications for disease control

-   [Short section, this is really a separate, less method-intense paper]
-   External control only has limitations; best used to reduce impacts, not prevent outbreak
-   Quantitative implications of stochasticity

### Methodological issues

-  Importance of lifting, restriction functions
-  Sensitivity to number of simulations and micro-step size, integration algorithms
-  Boundary conditions.  ODE can provide guidance in values. 


### Broad implications

-  EF-based optimization has applications on a wide variety of models and situations.
  -  Complement to discrete approaches, ensure agreement with theory.

-  Still a model reduction approach, but reduction takes place at a different scale
-  Advantages
    -   Continuous in time and space. 
    -   Interpretability.  Model-fitting approaches to reduced models leave us with uninterpretable values.
    -   Also applies to black-box processes
-  Challenges
    -   Computational cost (though many methods have this, and Hamiltonian approach avoids much of it)
    -   Lifting and restriction functions must be chosen carefulling.  Model reduction occurring at the micro-scale.
    -   Software
        - Some details of R package
   
### Future work

-   Comparison of multiple objectives.


Journals: Royal Society Interface? Ecological Economics?
Methods in Ecol and Evo.  AmNat for sepearate disease paper.  
