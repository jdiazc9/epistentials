rm(list = ls())

require(testthat)
require(ggplot2)
library(tidyr)
library(gtools)

structureData <- function(df, function_cols = 'auto') {
  
  # ----------
  # 
  # README: This function takes an input data frame (df) consisting of a set of
  # compositional variables (the "elements", say mutations, species, etc.)
  # which should take values 0 or 1, and a "function" variable which can take
  # any scalar value. If there are multiple "function" variables (columns),
  # these are automatically understood to be replicate measurements.
  # Each combination of elements should appear exactly once, that is, there
  # should not be more than one row in the data frame with the same
  # combination of 0s and 1s for the compositional variables.
  # Function columns can be automatically detected if function_cols is set to
  # "auto". If instead function_cols is a number or array of numbers, those
  # number(s) will be taken as the input index(es) for the function column(s).
  # 
  # The function returns a data frame with all the delta_functions for each
  # element across all potential backgrounds. If replicates are
  # provided/detected, means and standard deviations across them are returned.
  # The case where one or more values are unknown (NA) is valid, in these cases
  # means and standard deviations are computed across the non-NA replicates
  # only.
  # 
  # ----------
  
  # auto-detect "function" columns:
  # if there are more than two unique values, we consider it a function column
  # otherwise we consider it a compositional column
  # (note this can fail if e.g. function is binary)
  if (function_cols == 'auto') {
    function_cols <- sapply(1:ncol(df),
                            FUN = function(i) length(unique(df[, i])))
    function_cols <- which(function_cols > 2)
  }
  
  # compositional columns: all but the function_cols
  comp_cols <- setdiff(1:ncol(df), function_cols)
  elements <- colnames(df)[comp_cols]
  
  # sanity checks:
  # are all element names unique?
  if(length(elements) != length(unique(elements))) warning('Duplicated element names in data frame columns')
  # are compositional variables 0s and 1s?
  if(!all(sort(unique(unlist(df[, elements]))) == c(0, 1))) warning('Compositional variables should take values 0 and 1')
  # are all rows in df unique?
  if(!all(unique(df[, elements]) == df[, elements])) warning('Row repeats in input data fram: multiple instances of a same element combination. Function will be averaged.')
  
  # split data frame by replicate
  # each replicate is stored in a different element of a list object
  # function column is systematically named "fun" regardless of original name
  # (facilitates downstream analysis)
  # within each replicate, if there are repeated compositions (that is, 2 or
  # more rows with identical compositions), their corresponding functions
  # are averaged (aggregated via mean())
  df <- lapply(function_cols,
               FUN = function(i) aggregate(fun ~ .,
                                           data = cbind(df[, comp_cols], fun = df[, i]),
                                           FUN = mean))
  
  # compare backgrounds and +1's
  delta_df <- do.call(rbind,
                      lapply(df,
                             FUN = function(df_i) {
                               
                               do.call(rbind,
                                       lapply(elements,
                                              FUN = function(element_i) {
                                                
                                                backgrounds <- df_i[df_i[, element_i] == 0, ]
                                                plusones <- df_i[df_i[, element_i] == 1, setdiff(colnames(df_i), element_i)]
                                                compare <- merge(backgrounds, plusones,
                                                                 by = setdiff(elements, element_i),
                                                                 suffixes = c('_background', '_plusone'))
                                                compare <- cbind(compare[, elements],
                                                                 focal_element = element_i,
                                                                 compare[, c('fun_background', 'fun_plusone')])
                                                compare$delta_fun <- compare$fun_plusone - compare$fun_background
                                                
                                                return(compare)
                                                
                                              }))
                               
                             }))
  
  # aggregate replicates, get mean and standard deviation
  delta_df <- do.call(data.frame,
                      aggregate(delta_df,
                                cbind(fun_background, fun_plusone, delta_fun) ~ .,
                                FUN = function(x) c(mean = mean(x, na.rm = T), sd = sd(x, na.rm = T))))
  
  # output
  return(delta_df)
  
}

getFEEs <- function(df, function_cols = 'auto', mode = 'delta', reg_type = 'ols') {
  
  # ----------
  # Takes an input data frame and returns FEE slopes and intercepts for all
  # elements, either in the deltaF-vs-F (mode = 'delta') or in the F'-vs-F
  # (mode = 'plusone') representation. FEEs are fit either via ordinary
  # least-squares (OLS) regression (reg_type = 'ols') or via total least-squares
  # (TLS) regression (reg_type = 'tls').
  # ----------
  
  # structure data
  df <- structureData(df, function_cols)
  
  # wrapper function for ordinary least-squeres (OLS) or total least-squares (TLS) regression
  makeRegression <- function(x, y, reg_type = 'ols') {
    
    if (reg_type == 'ols') {
      mylm <- lm(y ~ x, data.frame(x = x, y = y))
      slope <- as.numeric(mylm$coefficients[2])
      intercept <- as.numeric(mylm$coefficients[1])
    } else if (reg_type == 'tls') {
      mypc <- prcomp(data.frame(x = x, y = y))$rotation
      slope <- mypc[2, 1]/mypc[1, 1]
      intercept <- mean(y) - mypc[2, 1]/mypc[1, 1]*mean(x)
    } else {
      warning('\'reg_type\' must be one of \'ols\' or \'tls\'')
    }
    
    return(c(slope = slope,
             intercept = intercept))
    
  }
  
  # what goes as y-axis?
  if (mode == 'delta') {
    y_col <- which(colnames(df) == 'delta_fun.mean')
  } else if (mode == 'plusone') {
    y_col <- which(colnames(df) == 'fun_plusone.mean')
  } else {
    warning('\'mode\' must be one of \'delta\' or \'plusone\'')
  }
  
  regs <- do.call(rbind,
                  lapply(unique(df$focal_element),
                         FUN = function(element_i) {
                           
                           reg <- makeRegression(df[df$focal_element == element_i, 'fun_background.mean'],
                                                 df[df$focal_element == element_i, y_col],
                                                 reg_type = reg_type)
                           return(data.frame(focal_element = element_i,
                                             slope = reg['slope'],
                                             intercept = reg['intercept']))
                           
                         }))
  
  return(regs)
  
}

plotFEEs <- function(df, function_cols = 'auto', mode = 'delta', reg_type = 'ols') {
  
  # ----------
  # README: This function takes an input (df), in "raw" format (not generated
  # by the structureData() function), and produces the most basic
  # version of the global epistasis plots.
  # The "mode" and "reg_type" parameters are defined as specified in the
  # getFEEs() function.
  # ----------
  
  # make regressions
  regs <- getFEEs(df, function_cols = function_cols, mode = mode, reg_type = reg_type)
  
  # structure data frame to set up plots
  df <- structureData(df, function_cols)
  
  # regressions and plots
  if (mode == 'delta') {
    
    base_plot <- ggplot(df, aes(x = fun_background.mean, y = delta_fun.mean)) +
      geom_abline(slope = 0,
                  intercept = 0,
                  color = 'gray')
    if (!all(is.na(df$fun_background.sd))) {
      base_plot <- base_plot +
        geom_errorbarh(aes(xmin = fun_background.mean - fun_background.sd, xmax = fun_background.mean + fun_background.sd),
                       alpha = 0.25)
    }
    if (!all(is.na(df$delta_fun.sd))) {
      base_plot <- base_plot +
        geom_errorbar(aes(ymin = delta_fun.mean - delta_fun.sd, ymax = delta_fun.mean + delta_fun.sd),
                      alpha = 0.25)
    }
    base_plot <- base_plot +
      geom_point() +
      facet_wrap(~ focal_element) +
      geom_abline(data = regs,
                  aes(slope = slope, intercept = intercept),
                  color = 'firebrick',
                  linewidth = 1) +
      scale_x_continuous(name = expression(italic(F)[italic(B)])) +
      scale_y_continuous(name = expression(Delta*italic(F))) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0))
    
  } else if (mode == 'plusone') {
    
    base_plot <- ggplot(df, aes(x = fun_background.mean, y = fun_plusone.mean)) +
      geom_abline(slope = 1,
                  intercept = 0,
                  color = 'gray')
    if (!all(is.na(df$fun_background.sd))) {
      base_plot <- base_plot +
        geom_errorbarh(aes(xmin = fun_background.mean - fun_background.sd, xmax = fun_background.mean + fun_background.sd),
                       alpha = 0.25)
    }
    if (!all(is.na(df$fun_plusone.sd))) {
      base_plot <- base_plot +
        geom_errorbar(aes(ymin = fun_plusone.mean - fun_plusone.sd, ymax = fun_plusone.mean + fun_plusone.sd),
                      alpha = 0.25)
    }
    base_plot <- base_plot +
      geom_point() +
      geom_blank(aes(y = fun_background.mean, x = fun_plusone.mean)) + # small hack for clean, equal x and y axis limits & scales
      facet_wrap(~ focal_element) +
      geom_abline(data = regs,
                  aes(slope = slope, intercept = intercept),
                  color = 'firebrick',
                  linewidth = 1) +
      scale_x_continuous(name = expression(italic(F)[italic(B)])) +
      scale_y_continuous(name = expression(italic(F)[+1])) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0))
    
  }
  
  return(base_plot)
  
}

getInterCoefficients <- function(df, function_cols = 'auto', mode = 'Taylor') {
  
  # This function takes an input data frame (df) (important: df should NOT have
  # been generated through structureData(), it should be a "raw" data frame) and
  # returns the interaction coefficients at all orders, meaning, if the function
  # F is expressed as:
  #    F = f_0 + f_1 x_1 + f_2 x_2 + ... + f_12 x_1 x_2 + ... + f_123 x_1 x_2 x_3 + ...
  # with x_i = 0,1 (Taylor) or x_i = -1,+1 (Fourier), then the function returns
  # the (fitted) values of the coefficients f_0, f_i, f_ij, etc.
  # This, in principle, can only be used reliably with combinatorially complete
  # data frames (otherwise it will return a bunch of NAs for some coefficients,
  # and provide unreliable estimates of the non-NA ones).
  # Important: in the input df, compositional values should take values 0 and 1,
  # even if the Fourier coefficients are to be returned. The conversion of 0's
  # into -1's is done internally by the function if 'mode' is set to 'Fourier'
  # (this conversion is not done if mode = 'Taylor').
  
  # structure replicates, if any (this whole chunk is borrowed from the
  # structureData() function, see details there)
  if (function_cols == 'auto') {
    function_cols <- sapply(1:ncol(df),
                            FUN = function(i) length(unique(df[, i])))
    function_cols <- which(function_cols > 2)
  }
  comp_cols <- setdiff(1:ncol(df), function_cols)
  elements <- colnames(df)[comp_cols]
  df <- lapply(function_cols,
               FUN = function(i) aggregate(fun ~ .,
                                           data = cbind(df[, comp_cols], fun = df[, i]),
                                           FUN = mean))
  
  # Taylor or Fourier coefficients
  if (mode == 'Fourier') {
    
    df <- lapply(df,
                 FUN = function(df_i) {
                   
                   comp_vars <- df_i[, 1:(ncol(df) - 1)]
                   comp_vars[comp_vars == 0] <- -1
                   return(cbind(comp_vars, df_i$fun))
                   
                 })
    
  }
  
  # fit linear model with interaction terms up to the N-th (with N being the
  # number of compositional variables)
  # note that this may fail if any of the compositional variables include the
  # character ':' in their name
  output_coefs <- do.call(rbind,
                          lapply(df,
                                 FUN = function(df_i) {
                                   
                                   my_fit <- lm(as.formula(paste('fun ~ ', paste(rep('.', length(comp_cols)), collapse = ' * '))),
                                                data = df_i)
                                   coefs <- my_fit$coefficients
                                   names(coefs)[1] <- ''
                                   coefs <- data.frame(order = c(0, 1 + nchar(names(coefs[2:length(coefs)])) - nchar(gsub(':', '', names(coefs[2:length(coefs)])))),
                                                       index = names(coefs),
                                                       value = as.numeric(coefs))
                                   
                                 }))
  output_coefs <- do.call(data.frame,
                          aggregate(value ~ .,
                                    data = output_coefs,
                                    FUN = function(x) c(mean = mean(x), sd = sd(x))))
  return(output_coefs)
  
}

getEffectiveInteractions <- function(df, function_cols = 'auto') {
  
  # This function takes a "raw" data frame (df) as input, and returns a list of
  # two data frames:
  # 
  #    effInter: this data frame contains the average functional effects (deltaF)
  #              and the average second-order interaction (epsilon) for every
  #              pair of elements i and j (averages are computed across all
  #              potential backgrounds). It also contains the corresponding
  #              pairwise "effective interactions" and their weights.
  #    expected_slopes: expected slope for each element (weighted sum of
  #              effective interactions of that element).
  # 
  # In both data frames, if replicates are provided, means and SDs are returned.
  
  # structure replicates, if any (this whole chunk is borrowed from the
  # structureData() function, see details there)
  if (function_cols == 'auto') {
    function_cols <- sapply(1:ncol(df),
                            FUN = function(i) length(unique(df[, i])))
    function_cols <- which(function_cols > 2)
  }
  comp_cols <- setdiff(1:ncol(df), function_cols)
  elements <- colnames(df)[comp_cols]
  df <- lapply(function_cols,
               FUN = function(i) aggregate(fun ~ .,
                                           data = cbind(df[, comp_cols], fun = df[, i]),
                                           FUN = mean))
  
  # define all pairs of elements
  pairs <- permutations(length(elements), 2, elements)
  
  # for every replicate...
  epsilons <- lapply(df,
                     FUN = function(df_i) {
                       
                       # for every pair of elements...
                       do.call(rbind,
                               lapply(1:nrow(pairs),
                                      FUN = function(p) {
                                        
                                        element_i <- pairs[p, 1]
                                        element_j <- pairs[p, 2]
                                        
                                        # backgrounds, +1's (+element i or j), and +2's (both elements i and j)
                                        fun <- rbind(cbind(df_i[df_i[, element_i] == 0 & df_i[, element_j] == 0, ], pos = 'bg'),
                                                     cbind(df_i[df_i[, element_i] == 1 & df_i[, element_j] == 0, ], pos = 'i'),
                                                     cbind(df_i[df_i[, element_i] == 0 & df_i[, element_j] == 1, ], pos = 'j'),
                                                     cbind(df_i[df_i[, element_i] == 1 & df_i[, element_j] == 1, ], pos = 'ij'))
                                        fun <- fun[, c(setdiff(elements, c(element_i, element_j)), 'fun', 'pos')]
                                        fun <- spread(fun, key = pos, value = fun)
                                        
                                        # deltas and epsilons
                                        fun$deltaF_i <- fun$i - fun$bg
                                        fun$deltaF_j <- fun$j - fun$bg
                                        fun$epsilon_ij <- fun$ij - fun$i - fun$j + fun$bg
                                        
                                        # adjust formatting and output
                                        fun <- cbind(0, 0,
                                                     element_i = element_i, element_j = element_j,
                                                     fun)
                                        colnames(fun)[1:2] <- c(element_i, element_j)
                                        fun <- fun[!is.na(fun$epsilon_ij), ]
                                        fun <- fun[, c('element_i', 'element_j', 'deltaF_i', 'deltaF_j', 'epsilon_ij', elements)]
                                        
                                        return(fun)
                                        
                                      }))
                       
                     })
  
  # effective interactions
  effInter <- lapply(epsilons,
                     FUN = function(eps) {
                       
                       eps <- aggregate(cbind(deltaF_i, deltaF_j, epsilon_ij) ~ element_i + element_j,
                                        data = eps,
                                        FUN = mean)
                       eps <- merge(eps,
                                    aggregate(deltaF_j ~ element_i,
                                              data = eps,
                                              FUN = function(x) sum(x^2)),
                                    by = 'element_i',
                                    suffixes = c('', '.sum_of_squares'))
                       eps$weight_ij <- eps$deltaF_j^2 / eps$deltaF_j.sum_of_squares
                       eps$effInter_ij <- eps$epsilon_ij / eps$deltaF_j
                       
                       return(eps)
                       
                     })
  
  # expected slopes
  expected_slopes <- lapply(effInter,
                            FUN = function(eff) {
                              
                              eff$weighted_effInter <- eff$weight_ij * eff$effInter_ij
                              out <- aggregate(weighted_effInter ~ element_i,
                                               data = eff,
                                               FUN = sum)
                              colnames(out) <- c('element', 'expected_slope')
                              
                              return(out)
                              
                            })
  
  # means and SDs across replicates
  effInter <- do.call(data.frame,
                      aggregate(cbind(deltaF_i, deltaF_j, epsilon_ij, weight_ij, effInter_ij) ~ element_i + element_j,
                                data = do.call(rbind, effInter),
                                FUN = function(x) c(mean = mean(x), sd = sd(x))))
  expected_slopes <- do.call(data.frame,
                             aggregate(expected_slope ~ element,
                                       data = do.call(rbind, expected_slopes),
                                       FUN = function(x) c(mean = mean(x), sd = sd(x))))
  
  # output
  return(list(effInter = effInter,
              expected_slopes = expected_slopes))
  
}







                      
                                               
                      





### EXAMPLES

if (T) {
  
  function_cols <- 'auto'
  df <- read.csv('../data_sets/pyoverdine-training_Diaz-Colunga2024.csv')
  gedf <- structureData(df)
  plotFEEs(df, mode = 'delta')
  plotFEEs(df, mode = 'plusone', reg_type = 'tls')
  effInter <- getEffectiveInteractions(df)[['effInter']]
  expected_slopes <- getEffectiveInteractions(df)[['expected_slopes']]
  
  # df <- read.csv('../data_sets/amyl_Sanchez-Gorostiaga2019.csv')
  # plotFEEs(structureData(df), mode = 'delta')
  # plotFEEs(structureData(df), mode = 'plusone')
  
}












### TESTS


if (F) {
  
  # check that all element names are unique
  # otherwise return warning & halt execution
  elements <- colnames(df)[comp_cols]
  test_that("Are element names (mutations/species/etc.) unique in the data frame?", {
    expect_true(length(elements) == length(unique(elements)))
  })
  
  # check that compositional variables are 0s and 1s
  # otherwise return warning & halt execution
  test_that("Are compositional variables 0s and 1s?", {
    expect_true(all(sort(unique(unlist(df[, elements]))) == c(0, 1)))
  })
  
  # check that each unique combination appears exactly once
  # otherwise return warning & halt execution
  test_that("Are there duplicate element combinations in the input data frame?", {
    expect_true(all(unique(df[, elements]) == df[, elements]))
  })
  
}




