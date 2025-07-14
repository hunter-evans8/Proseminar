# Code
This section will be the raw code, broken up into steps. For insights and shortened breakdown look at the README section.



## Simulation Code
```r
doit <- function(n, beta0, beta1, theta, iterations){
  data_out <- tibble(y = numeric(), 
                     x = numeric(),
                     n = numeric(),
                     beta0 = numeric(),
                     beta1 = numeric(),
                     theta = numeric())

  poisson_out <- tibble(term = character(),
                        estimate = numeric(),
                        std.err = numeric(),
                        statistic = numeric(),
                        p.value = numeric(),
                        iteration = integer(),
                        n = numeric(),
                        beta0 = numeric(),
                        beta1 = numeric(),
                        theta = numeric())
  
  negbin_out <- tibble(term = character(),
                       estimate = numeric(),
                       std.err = numeric(),
                       statistic = numeric(),
                       p.value = numeric(),
                       iteration = integer(),
                       n = numeric(),
                       beta0 = numeric(),
                       beta1 = numeric(),
                       theta = numeric())  
  
  for (k in 1:iterations){

  x <- rnorm(n = n, mean = 0, sd = 1)
  mu <- exp(beta0 + beta1 * x)
  y <- MASS::rnegbin(n, mu = mean(mu), theta = theta)
  
  data <- tibble(y, x) %>%
    mutate(iteration = k,
           n = n,
           beta0 = beta0,
           beta1 = beta1,
           theta = theta)
  
  data_out <- rbind(data_out, data)
    
# Fit Poisson model
  poisson_model <- tidy(glm(y ~ x, data = data, family = "poisson")) %>%
    mutate(iteration = k,
           n = n,
           beta0 = beta0,
           beta1 = beta1,
           theta = theta)

# Fit Negative Binomial model
  negbin_model <- tidy(MASS::glm.nb(y ~ x, data = data)) %>%
    mutate(iteration = k,
           n = n,
           beta0 = beta0,
           beta1 = beta1,
           theta = theta)

  poisson_out <- rbind(poisson_out, poisson_model)
  negbin_out <- rbind(negbin_out, negbin_model)
  }
  
data_out <<- data_out  
poisson_out <<- poisson_out
negbin_out <<- negbin_out
}

set.seed(68340) 
doit(n = 100, beta0 = 1.5, beta1 = 0.25, theta = 999999, iterations = 100)
write.csv(data_out, " hunterevans/data_n_100_beta0_15_beta1_0_25_theta_999999.csv")
write.csv(poisson_out, "hunterevans/1-data/poisson output/poisson_n_100_beta0_15_beta1_0_25_theta_999999.csv")
write.csv(negbin_out, "hunterevans/1-data/negbin output/negbin_n_100_beta0_15_beta1_0_25_theta_999999.csv")

```


## Analysis Code 
```r
library(tidyverse)

data <- read_csv("hunterevans/1-data/simulated data/Copy of data_n_75_theta_10.csv")

summary <- data %>% 
  group_by(iteration) %>%
  summarize(mean(y), var(y))

summary %>% ggplot(aes(x = `mean(y)`, 
                       y = `var(y)`)) +
  geom_point() + 
  geom_abline(intercept = 0) +
  theme_bw() +
  ylim(0, 20) + xlim(0, 20) + 
  labs(x = "mean", 
       y = "variance", 
       title = "n=75, theta=10")

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(broom) # For tidying model output

# Group by dataset, nest, and apply the Poisson regression
poisson_models <- data %>%
  group_by(iteration) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(y ~ x, family = poisson, data = .x)),
    model_summary = map(model, broom::tidy) # Optional: to extract model summaries
  )

# To view results, unnest or access summaries
# For example, view the model summaries for each dataset:
poisson_summaries <- poisson_models %>% 
  unnest(model_summary) %>%
  filter(term == "x") %>%
  mutate(reject = if_else(p.value < 0.05, 1, 0)) %>%
  ungroup()

poisson_summaries %>% count(reject)

979100/(9791+209)

# Group by dataset, nest, and apply the negative binomial regression
nb_models <- data %>%
  group_by(iteration) %>%
  nest() %>%
  mutate(
    model = map(data, possibly(~ MASS::glm.nb(y ~ x, data = .x))),
    model_summary = map(model, broom::tidy), # Optional: to extract model summaries
    aic = map_dbl(model, ~ if(!is.null(.x)) AIC(.x) else NA_real_),  # Calculate AIC
    bic = map_dbl(model, ~ if(!is.null(.x)) BIC(.x) else NA_real_)   # Calculate BIC
  )

# To view results, unnest or access summaries
# For example, view the model summaries for each dataset:
nb_summaries <- nb_models %>% 
  unnest(model_summary) %>%
  filter(term == "x") %>%
  mutate(reject = if_else(p.value < 0.05, 1, 0)) %>%
  ungroup()

nb_summaries %>% count(reject)

979100/(9791+209)
956900/(9569+431)






# Define a "safe" version of glm.nb using possibly, with NULL as the fallback
safe_glm_nb <- possibly(
  ~ glm.nb(y ~ x, data = .x, control = glm.control(maxit = 5000)),
  otherwise = NULL
)

# Filter out groups with low variance in `y` before fitting
nb_models <- data %>%
  group_by(iteration) %>%
  nest() %>%
  # Filter datasets with sufficient variation in y
  filter(map_lgl(data, ~ var(.x$y, na.rm = TRUE) > 0)) %>%
  mutate(
    model = map(data, safe_glm_nb),
    model_summary = map(model, ~ if(!is.null(.x)) broom::tidy(.x) else NULL)
  )

# View successful models
nb_summaries <- nb_models %>%
  filter(!is.null(model_summary)) %>%
  unnest(model_summary)

# Identify failed groups
failed_groups <- nb_models %>%
  filter(is.null(model_summary)) %>%
  select(iteration)


```

## Graphing Code

### Poisson PMF Graph

```r
lambda <- 10
x_vals <- 0:25 

# Compute the PMF values
pmf_vals <- dpois(x_vals, lambda)

# Plot the PMF
barplot(pmf_vals,
        names.arg = x_vals,
        col = "steelblue",
        main = paste("PMF of Poisson Distribution (λ =", lambda, ")"),
        xlab = "Number of Events (x)",
        ylab = "P(X = x)")

```



### Negative Binomial Graphing Code 

```r
# Parameters
r <- 10  # number of successes
p <- 0.5 # probability of success

# Define range of x (number of failures before 10 successes)
x_vals <- 0:30

# Compute the PMF
pmf_vals <- dnbinom(x_vals, size = r, prob = p)

barplot(pmf_vals,
        names.arg = x_vals,
        col = "red",
        main = paste("PMF of Negative Binomial (r =", r, ", p =", p, ")"),
        xlab = "Number of Failures (x)",
        ylab = "P(X = x)")
```


## Bias vs n values theta graph

```r
library(tidyverse)

bias_df <- combined_results %>%
  filter(term == "x") %>%  # 'x' is the name of beta1's term in tidy() output
  group_by(model, n, theta) %>%
  summarise(mean_estimate = mean(estimate),
            true_beta1 = unique(beta1),
            bias = mean_estimate - true_beta1,
            .groups = "drop") %>%
  mutate(theta = factor(theta))

# Plot
ggplot(bias_df, aes(x = n, y = bias, color = theta, group = theta)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~model, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = unique(bias_df$n)) +
  labs(title = "Bias of β₁ Estimate vs Sample Size",
       x = "Sample Size (n)",
       y = "Bias",
       color = expression(theta)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 14))


```

### MSE vs n values theta 

```r
library(tidyverse)

# Combine simulation results from all runs
# Make sure both poisson_out and negbin_out include a 'model' column
poisson_out$model <- "Poisson"
negbin_out$model <- "Negative Binomial"

combined_results <- bind_rows(poisson_out, negbin_out)

# Filter for beta1 term (usually "x" in tidy() output)
mse_df <- combined_results %>%
  filter(term == "x") %>%
  group_by(model, n, theta) %>%
  summarise(
    true_beta1 = unique(beta1),
    mse = mean((estimate - true_beta1)^2),
    .groups = "drop"
  ) %>%
  mutate(theta = factor(theta))

# Plot
ggplot(mse_df, aes(x = n, y = mse, color = theta, group = theta)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~model, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = unique(mse_df$n)) +
  labs(title = "MSE of β₁ Estimate vs Sample Size",
       x = "Sample Size (n)",
       y = "Mean Squared Error",
       color = expression(theta)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 14))


```


### Average Standard Error Graphs

```r
library(tidyverse)

# Tag model
poisson_out$model <- "Poisson"
negbin_out$model <- "Negative Binomial"

# Combine data
combined_results <- bind_rows(poisson_out, negbin_out)

# Filter: β₁ only ("x") and n = 25
se_df <- combined_results %>%
  filter(term == "x", n == 25) %>%
  group_by(model, theta) %>%
  summarise(avg_stderr = mean(std.err), .groups = "drop")

# Plot with x-axis breaks every 500
ggplot(se_df, aes(x = theta, y = avg_stderr, color = model)) +
  geom_line(size = 1.3) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  scale_x_continuous(breaks = seq(0, max(se_df$theta), by = 500)) +
  labs(title = "Average Standard Error of β₁ vs Dispersion (θ)",
       x = expression(theta),
       y = "Average Std. Error of β₁",
       color = "Model") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

```