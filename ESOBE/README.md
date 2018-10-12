# ESOBE 2018

Hello! In this repository you will find the complementary material to my work _"Redistributive effect of monetary policy: evidence from Brazil"_, which is essentially my (future) master thesis.

## File description

* [Paper](https://github.com/aishameriane/Presentations/blob/master/ESOBE/Schmidt_Moura23072018.pdf) - The version submitted to the event.
* [Poster](https://github.com/aishameriane/Presentations/blob/master/ESOBE/P%C3%B4ster%20ESOBE.pdf) - PDF file of the poster.
* [New version] - I'm working on a newer version, but still needs some adjustmets in the text and results. If you keep scrolling on this readme, you will find some explanations on what I am trying to do now.
* [Codes] - In the paper I used a combination of Matlab and R code. It still has some parts of the comments in Portuguese, so as soon as I finish fixing this, will put in here. My intention is to migrate everything to R, but for now the estimation is being computed on Matlab, while de data download, cleaning and post processing is made on R. If you want the code, please send me an email.

## What was done so far

This paper is about evaluating the impact of conventional monetary policy on capital-labor ratio using data from Brazil in a TVP-VAR with Wishart stochastic volatility. I wanted to do something inspired on [Mumtaz and Theophilopoulou](//www.sciencedirect.com/science/article/pii/S0014292117301332) (2017), but using data from Brazil. At the same time, I wanted to avoid using Primiceri (2005) paper due to the "ordering problem" inherent to this particular model. That's when my thesis became "extending Uhlig's (1997) model to a TVP-VAR (the stochastic volatility was already included in the original model). Since Uhlig's propositions do not allow for this generalization in a straighforward form, I am borrowing Windle and Carvalho (2014) theorems to implement the filtering, smoothing and prediction steps for the covariance matrices and I combine this with a Carter and Kohn (1994) algorithm in order to estimate the TVP coefficients. My argument to use this model is that it allows for changes in volatility (as has been observed in Brazil macroeconomic aggregates) and the model coefficients, which I would expect to vary as well - at least some of them, for example, the inflation parameters in the interest rate equation.

The use of the capital-labor ratio is motivated (or at least I tried tomotivate) in the paper. Essentially this quantity can be seen as the quocient between the share of the capital and share of labor incomes with respect to the GDP. If everyone had the same proportion of income from capital and income from labor, changes in capital-labor ratio could be seen like some sort of household portfolio (income) rebalancig. But this is not the case, i.e., richer people tend to have more capitals income than other people, so changes in capital income ratio would imply a redistribution of income, in some sense.

## What I am going to do next

### In the data

* I used an automatic method do determine the number of lags in the model. Something nice to see is the effect of including or removing lags, because it is possible that the capital labor ratio takes more time to respond to a monetary shock (it would be the case if the capital income have a large portion coming from physical capital);
* I've made some transformations in the data: per capita GDP is annual variation of the Per capita GDP and exchange rate is the monthly variation of the ex. rate. I did this to have stationary series, but I'm not sure if I'm not loosing information in here and some more investigation would be helpful.
* This is a sugestion from prof. Meurer: if I could get some data similar to the capital labor ratio, but discriminating by types of capital income sources, I could estimate (probably a simpler) model to see if my results hold;

## In the econometric model

* I could try to use Mateus (my classmate) model (adapted for my covariance speficication) to obtain residuals - But this would imply adapting his theorems to my situation. Mateus' model is almost the same TVP-VAR as mine, except from the volatility in the state transition equation that is a little bit more general. 
* Results so far indicates three types of coefficients: 1) constant and zero; 2) constant and non zero; 3) TVP. Shrinkage methods would be nice to use in this cada
* Right now I'm using Korobilis code to compute impulse response functions, but I want to make something similar to the IRF in the bvarsv models, where you can use the time, impulse and response functions as arguments;
* Something further that I want to study is this SVAR with TVP from [Bognanni](https://www.clevelandfed.org/en/newsroom-and-events/publications/working-papers/2018-working-papers/wp-1811-a-class-of-time-varying-parameter-structural-vars.aspx) that may be interesting in this particular problem;
* Explore stochastic volatility on the state transition equation;
* Some measures of goodness-of-fit

## Other

* Improve the analysis because right now it is lacking some explanations and graphs. For example, it may be of interest to compute the IRF of inflation on K/L;
* I need to work on the efficiency of the algorithm;
* Finish migrating the entire algorithm to R;
* Identify which coefficients are TVP or not.
