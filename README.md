# Overview

Welcome to my project from my research as an undergraduate student. This project highlights skills in R, such as graphing and simulating data, and has some heavy statistical theory.

Here is the full paper I have written about this topic and presented in front of a math commitee: [Hunter_Proseminar_Final.pdf](Hunter_Proseminar_Final.pdf)


Here is a link to the raw code that I have compiled: [01_Code.md](01_Code.md)

This section is just going to go over the main idea behind this project and the results we got from it. All graphs and results can be found within the paper above.


# The Objective

 This paper aims to explore the implications of analyzing overdispersed count data under the
 Poisson distribution, which breaks the equidispersion assumption. Using negative binomial as a comparison,
 our goal is to understand the potential consequences of model misspecification and to guide researchers
 toward more robust statistical practices.


# Methods

We conducted a simulation study to compare Poisson and negative binomial. We wrote a function
that would repeatedly generate datasets while also allowing a way to change the relationship between µ and
σ2. We did this by simulating our outcome variable as y ∼ NegBin(µ,θ), where θ is the dispersion parameter
that changes the relationship between µ and σ2. After the data was simulated, each dataset was used to
construct a model under both Poisson and negative binomial distributions. The simulated raw data and,
separately, the analysis results from both Poisson and negative binomial regressions were stored for further
analysis. The regression analysis results were examined with respect to bias, mean square error (MSE), type
I error, and model selection using AIC.

# The Analysis

## Simulation Study 

 To explore the differences in these models, a simulation study was conducted. A simulation study is the
 process of creating synthetic, or made up, data. This was done using R with the packages tidyverse,
 MASS, and broom.
 Our goal was to make a function that would repeatedly generate datasets while also allowing a way to
 change the relationship between µ and σ2. The predictor, x, was simulated first, with x ∼ N(0,1). Then
 the mean was defined as µ = exp{β0 +β1x,}, leading to simulating the count outcome, y ∼ NegBin(µ,θ).
 Note that θ is the dispersion parameter and allows us to change the relationship between µ and σ2. After
 the data was simulated, each dataset was used to construct a model under both the Poisson and negative
 binomial distributions. Regardless of Poisson or negative binomial, we modeled
 ln(y) = β0 +β1x.
 The simulated raw data and, separately, the analysis results from both Poisson and negative binomial
 regressions were stored for further analysis. The regression analysis results were examined with respect to
 bias, mean square error (MSE), type I error, and model selection using AIC.

## Definition of Parameters
 For the simulation, we held the intercept and slope constant at β0 = 1.5 and β1 = 0.25, respectfully. Then, we
 considered three sample sizes, n ∈ {25,100,500}, and five dispersion parameters, θ ∈ {1,10,50,500,2000}.
 The dispersion parameters correspond to the following scenarios:
 • θ =2000: µ =σ2
 • θ =500: µ≈σ2
 • θ =50: µ<σ2
 • θ =10: µ<<σ2
 • θ =1: µ<<<σ



# Conclusion 
 In our simulation study, our results agreed with statistical theory. The bias and MSE results showed that
 there wasn’t a major difference between the two models, which was expected as the slopes themselves
 aren’t affected when the data is overdispersed. In the standard error results, we saw that negative binomial
 had a higher SE than Poisson for θ ∈ {1,10,50}, where the data was overdispersed. This visualizes the
 underestimation of the standard error.
 Wesawsimilar results in the rejection rate graphs as well, which showed when the data was overdispersed,
 the Poisson model would reject the null, H0 : β = 0, 30-40% of the time. In the vast majority of the simulated
 datasets, the AIC suggested that negative binomial was the better fit when the data was overdispersed. Based
 on these results, we can conclude that when the data is overdispersed and Poisson is used, our conclusions
 may not be accurate.
 One thing to highlight is that in the θ = 500 case, where µ ≈ σ2, Poisson and negative binomial are
 behaving similarly. Because there is not a clear advantage in estimating an additional parameter with the
 negative binomial, we recommend using Poisson regression for simplicity

 ## Further Study

  To further this research, this study can be tailored to a specific scientific application using real world data
 and parameters. To better examine what is happening between θ = 50 and θ = 500, additional simulations
 can be performed with 50 ≤ θ ≤ 500. We can also look into overdispersion in other models, such as the
 quasi-Poisson, and zero-inflated models.