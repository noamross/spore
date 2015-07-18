# Optimal Control in Individual-Based Models using Multi-Scale Approaches

## Introduction

A central challenge in disease ecology is developing effective management strategies to control harmful diseases and limit their damage to wildlife, human health, and ecosystem services. Using models of the dynamics of hosts and disease, this can be formulated as an economic optimization problem: how to maximize the benefits of a strategy, or how to minimize the costs of meeting a disease control goal.  

Disease management may be considered an economic problem: Given limited resources, how to best manage disease to minimize impacts? Or, conversely, given management targets, how to best limit the cost of disease control? Research addressing such questions - economic epidemiology - is a growing field [@Brandeau2003; @REFS], though the literature to on the economics of wildlife disease is more limited. [@Horan2005; Bicknell1999; Fenichel2010].  The bioeconomics of invasive species management is more broadly explored and has many close parallels with the economics of disease management [@REF (Epanchin-Niell, Liebhold].

  [TODO See additional citations in OPENQUALS for variety of methods applied to invasives]

Individual-based models (commonly called agent-based models) simulate  of systems based on the behavior and properties of individuals, rather than the population as a whole.  These can include network contact structure, spatial heterogeneity, individual variation, all of which affect disease transmission.  Frequently such models can not be expressed as a closed-form set of equations at the population level without making strong assumptions about the distribution of individuals.

Importance of individual-based processes (network contact structure, host and pathogen variation) in disease dynamics.
-   IBMs work where population-level models do not capture appropriate dynamics or answer questions. [@Grimm2005]

Individual-based models have have been used in the management of disease in a "scenario" framework.

Agent-based models my include optimizing agents, rather than be optimized for management, in order to simulate the behavior of individuals.  

### Other approaches

One approach to optimization of IBMs has been to apply control techniques
to reduced-form models.  The optimization problem can be solved on these simplified ODEs.  When the reduced-form model represents a mean-field approximation, it may perform poorly under when the underlying system is heterogenous [@Federico2013].  Another approach is to generate population-level models with non-mechanistic equations that approximate IBM behavior.  Such methods can yield better results than mean-field models. However, the applicability of these adjusted models is limited the the range over which they were parameterized, and model parameters are difficult to interpret [Oremland 2015].

Optimal control, as practices by 

@Pontryagin1967 described the necessary conditions to determine the optimal control path for a dynamical system.  The application of these methods to dynamical ecological systems is detailed in @Clark1990.  In brief, the optimal control path can be derived by solving maximization problem:

$$\max_{h >= 0} \int_{t=0}^T \pi(x, h, t) \, dt\, \text{ s.t. } \dot N, \dot P$$.

Various numerical methods can be used to find the optimal control path.  Where $f(t)$ has a closed-form, 

-  Applying optimal control techniques to a simplified mean-field model [@Federico2013]
   can work for homogenous models, but yields non-optimal results under heterogeneity.
-  
-  If reduced form models have low enough dimensionality, discretization can allow for
   approximate optimal control through  for heuristic methods that search the entire space of possible controls.[@Oremland2015].
-  In addition to optimal control, IBMs can be useful in comparison of a small set of management
scenarios [LOTS OF EXAMPLES].

### Introduction to Equation-Free Modeling

-   IBM dynamics can be difficult to capture as ODEs [@Brown2004]

-   Equation-free (EF), or multi-scale, modeling [@Armaou2004, but see @Kevrekidis2009a for review] is a method for capturing and analyzing population-scale dynamics of individual-based models while bypassing the deriving population-scale equations, which exist conceptually but often can not be derived.  It has been used in a variety of applications, including physical chemistry, movement ecology and epidemiology [@REFS, Levin, @Williams2015]

In equation-free modeling, the population-scale system is simulated using repeated short-repeated bursts of the individual-based model.  At each time step, a "lifting" operator maps the population-scale state to an ensemble of possible individual-scale states.  Each of these instances of the individual-based model are simulated for a short burst of time, then are each reduced back to a population-level state with a "restricting" function.  The results from all simulations are average to generate a mean result.  

The system dynamics may be simulated by the lift-simulate-constrict cycle alone.  However, it is computationally advantageous to use in conjunction with a projection step.  In this approach, rather than simulating for an entire time-step, the simulation is run for a fraction of a time step to obtain an estimate the time derivative of the population-scale system.  The population-scale state is then projected forward using this value in a differential equation solver using methods such forward Euler or Runge-Kutta.  In addition these approximate derivatives may be used 

For a stochastic system, the EF framework simulates the expected value of the simulation.  Each micro-step consists of many simulation of 


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

The model is a host-pathogen system similar to the macroparasite model of @Anderson1978, where individuals can host multiple pathogens of the same type, and suffer increased mortality with greater infection load.  The primary differences are the separation of birth and death processes and influx of infectious particles from outside the system.

In the mean-field version of the system, the host population $N$ increases via density-dependent reproduction ($rN(1-N/K)$) and decreases via a constant intrinsic mortality rate ($d$), and and additional mortality for each infection in a host ($\alpha$). 

$$\frac{dN}{dt} &= r N (1 - N/K) - \alpha P - d N$$

The pathogen population $P$ grows with new infections that are density-dependent $\lambdaPN$, and decreases via pathogen mortality ($\mu P$), background host mortality ($dP$), and mortality due to the disease, which affects the most-infected hosts most $alpha(P + P/N)$.  The form of of this term depends on the assumption of a random (Poisson) distribution of pathogens among hosts.  Finally, additional parasites enter the system via external propagule pressure $\lambda_{ex}$. $\lambda_{ex}$ may be reduced by control effort $h(t)$, which reduces propugule pressure by a factor of $e^{h(t)}$. 

$$\frac{dP}{dt} &= \lambda P N - \mu P - d P - \alpha P - \alpha P^2/N + e^{-h} N \lambda_{ex}$$

It assumption of (1) a large population size adequately represented by continuous variables, (2) a constant distribution of parasites among individuals, and (3) no stochastic processes.

The individual-based version of this represents the population as discrete counts of hosts and parasites. In the IBM, the state of the system is represented as a vector of individual hosts, each with a discrete number of infections.  Births, deaths, new infections and loss of infections (recovery) in each individual $i$ with number of infections $j$, occur stochastically according  to stochastic rates $r$:

$$\begin{aligned}
r_{i, birth} &= r(1-N/K) \\
r_{i, death} &= \alpha j + d\\
r_{i, infection} &= \lambda * \sum_{i=1}^N j_i + lambda_{ex} e^{h(t)},
r_{i, recovery} &= \mu j
\end{aligned}$$

This stochastic process occurs in continuous time, implemented via Gillespie's stochastic simulation algorithm (SSA) [@Gillespie1976].

Multi-scale simulation requires "lifting" and "restriction" functions, which map macro-scale system states to the micro-scale and vice-versa.  In this case, I generate micro-level states from the macro state (lifting) by random sampling (with replacement) of $P$ individuals to distribute the individuals (a Poisson process).   For the reverse (restriction), I simply sum the total host and parasite populations.

### Control problem

The management problem is to maximize the host population provides an ecosystem service of value $v$ per unit time, proportional to the host population size.  The control $h$, which consists of reducing the arrival rate of new infections, has a cost $c$.   Thus, the optimization problem is to maximize the net benefits over the course of a time period $T$, subject to the dynamics of the system.

$$\max_{h >= 0} \int_{t=0}^T  \pi(x,h,t), dt\, \text{ s.t. } \frac{dx}{xt} = f(x, h, t) $$.

From this I derive the Hamiltonian equation, which must be optimized for all time points:

$$\begin{align}
\mathcal{H} = &vN - ch + \\
&\nu_1 \left( rN(1 - N/K) - \alpha P - d N \right) + \\
&\nu_2 \left(\lambda P N - \mu P - d P - \alpha P - \alpha (P^2)/N + e^{-h}  N \lambda_{ex} \right)
\end{align}$$

Using the adjoint principal, I can derice the paths of the shadow values $\nu_1$ and $\nu_2$:

$$\begin{align}
\frac{d\nu_1}{dt} &= -v - \nu_1 \left(r - d - 2r \frac{N}{K}\right) - \nu_2 \left(\lambda P + \alpha \frac{P^2}{N^2} + e^{-h} \lambda_{ex} \right) + \delta \nu_1 \\
\frac{d\nu_2}{dt} &= \nu_1 \alpha - \nu_2 \left(\lambda N - \mu - d - \alpha - \frac{2 \alpha P}{N} \right) + \delta \nu_2
\end{align}$$




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

For the equation-free simulation approach, the system dynamics are a "black box" which may be simulated but under which the derivatives $\dot N$ and $\dot P$ are unknown.

At each time step, I simulate the system using the lift-simulate-restrict cycle repeatedly in order optimize

The derivative of the shadow states can be derived via the *adjoint principle*:

$$\frac{d\nu}{dt} = - \frac{d\mathcal H}{dx} = -\frac{d\pi}{dx} - \nu \frac{df}{dx}$$

The profit function $\pi(\dot)$ is known, so $\frac{d\pi}{dx}$ can be determined analytically.  To calculate $\frac{df}{dx}$, 


While the equation-free framework can be used on "black box" simulators, taking advantage of some properties of the individual-based model can improve computational efficiency and performance of the method.  In this case, I do so in two ways.  First, as the IBM uses the Gillespie SSA, I use a single variable-time "Gillespie" micro-step rather than a fixed-period micro-step.  This allows me to estimate $dx/dt$ without any micro-step offset.  Second, I take advantage of the mapping between continuous population states at the macro-level and discrete counts at the micro-level.  When lifting a continuous macro-level to create an ensemble of random micro-level states, some will have populations rounded down from the continuous level, and some rounded up.  The difference in $x'$ between these populations is used to determine $x'/x$ without additional simulations.




-   Use EF simulation to simulate dynamics of system as well as estimate derivatives via finite differencing

For a single variable system $df / dx$, can be approximated if the second derivative of the state variable is calculated by the integrator.  That is $$\frac{df(\dot)}{dx}  / dx \approx \frac{\Delta f(\dot) / \Delta t}{\Delta x / \Delta t}$.  However, for higher-dimensional systems, this approximation fails. Instead the full matrix of $$\frac{df(\dot)}{dx}$$ (the Jacobian matrix) is approximated by perturbing the simulation and measuring the difference in outcomes. 


-   Derive adjoint equations as functions of numerical estimates of host and and pathogen dynamics
-   Simulate optimal path by solving for Hamiltonian maximum at each time step

-  Lifting and restriction functions translate between full distribution of infections and summary statistics
   of distribution
   -   Use Conway-Maxwell-Poisson to allow for both under- and over-dispersed distributions with changing
       aggregation.

## Results

In the absence of control, the the disease invades the system rapidly.  This suppresses the population via and increased overall mortality rate for hosts, and the system reaches a stable equilibrium.

In the ODE system, the internal equilibrium can be calculate as [TODO].

The individual-based, stochastic model results in figure.  Shown are the mean path and the 95% confidence intervals for the system state through time, 

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
    -   Lifting and restriction functions must be chosen carefully.  Model reduction occurring at the micro-scale.
    -   Software
        - Some details of R package
   
-  Stochasticity:  The equation-free model produces an estimate of the IBM's *mean* behavior, and thus an optimal control path for the mean conditions. The primary utility of the method is not in stochastic control.  [See a host of other methods for this: Talk to Carl/Jim].  Rather, its value lies in optimizing where the IBM behaviors are not easily expressed in closed form equations.  Stochastic behavior makes the method computationally difficult.  The variation in outcomes is amplified in the estimate of the derivatives and even more so in the estimation of of the Jacobian matrix. 
### Future work

-   Comparison of multiple objectives.


Journals: Royal Society Interface? Ecological Economics?
Methods in Ecol and Evo.  AmNat for sepearate disease paper.  
