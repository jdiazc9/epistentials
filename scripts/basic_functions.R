require(testthat)
require(ggplot2)
require(tidyr)
require(gtools)
require(combinat)
require(MASS)

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
  
  rownames(regs) <- NULL
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
                   
                   comp_vars <- df_i[, 1:(ncol(df_i) - 1)]
                   comp_vars[comp_vars == 0] <- -1
                   return(cbind(comp_vars, fun = df_i$fun))
                   
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
  rownames(output_coefs) <- NULL
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
  return(effInter)
  
}















##################################################################
################### EXTREMELY HACKY CODE BELOW ###################
######################### NEEDS UPDATING #########################
##################################################################








plotLandscape <- function(df, function_cols = 'auto') {
  
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
  
  # mean and SD of function across replicates
  df <- do.call(data.frame,
                aggregate(fun ~ .,
                          data = do.call(rbind, df),
                          FUN = function(x) c(mean = mean(x), sd = sd(x))))
  
  # number of elements
  n_elements <- rowSums(df[, elements])
  
  # is an input assemblage a direct descendant of another?
  isDirectDescendant <- function(this_assemblage, of_this_assemblage) all(this_assemblage >= of_this_assemblage) & sum(this_assemblage) == (1 + sum(of_this_assemblage))
  
  # make edges of fitness graph
  edges <- do.call(rbind,
                   lapply(1:nrow(df),
                          FUN = function(i) {
                            
                            source_assemblage <- df[i, elements]
                            target_assemblages <- df[sapply(1:nrow(df),
                                                            FUN = function(j) isDirectDescendant(df[j, elements], source_assemblage)),
                                                     elements]
                            
                            if(nrow(target_assemblages)) {
                              edges <- data.frame(source = paste(source_assemblage, collapse = ''),
                                                  target = sapply(1:nrow(target_assemblages), FUN = function(j) paste(target_assemblages[j, ], collapse = '')),
                                                  source.nmut = n_elements[i],
                                                  target.nmut = 1 + n_elements[i])
                            } else {
                              edges <- data.frame(source = character(0),
                                                  target = character(0),
                                                  source.nmut = numeric(0),
                                                  target.nmut = numeric(0)) 
                            }
                            
                            return(edges)
                            
                          }))
  
  edges <- cbind(edge_id = paste('edge_', 1:nrow(edges), sep = ''),
                 edges)
  
  # plot graph
  ### FIXME: this NEEDS to be optimized
  mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')
  
  n_mut <- ncol(df) - 1
  
  landscape <- data.frame(genot = sapply(1:nrow(df),
                                         FUN = function(i) paste(df[i, elements], collapse = '')),
                          f = df$fun.mean,
                          f.sd = df$fun.sd)
  
  dfl <- cbind(edges,
               source.f = setNames(landscape$f, landscape$genot)[edges$source],
               target.f = setNames(landscape$f, landscape$genot)[edges$target],
               source.f.sd = setNames(landscape$f.sd, landscape$genot)[edges$source],
               target.f.sd = setNames(landscape$f.sd, landscape$genot)[edges$target])
  
  dfl$color <- 'A'
  
  dfx <- gather(dfl[, c(1, 4, 5)], position, nmut, source.nmut:target.nmut)
  dfx$position <- setNames(c('source', 'target'), c('source.nmut', 'target.nmut'))[dfx$position]
  
  dfy <- gather(dfl[, c(1, 6, 7)], position, f, source.f:target.f)
  dfy$position <- setNames(c('source', 'target'), c('source.f', 'target.f'))[dfy$position]
  
  dfxy <- merge(dfx, dfy, by = c('edge_id', 'position'))
  
  dfl <- merge(dfxy, dfl[, c('edge_id', 'color')], by = 'edge_id')
  
  dy <- min(c(max(landscape$f) - landscape$f[1], landscape$f[1] - min(landscape$f)))
  dy <- round(dy/0.1)*0.1
  ybreaks <- seq(landscape$f[1] - 10*dy, landscape$f[1] + 10*dy, by = dy)
  
  myplot <-
    ggplot(dfl, aes(x = nmut, y = f, group = edge_id, color = color)) +
    geom_abline(slope = 0,
                intercept = 0,
                color = '#d1d3d4') +
    geom_line() +
    scale_x_continuous(name = '# of elements',
                       breaks = 0:n_mut,
                       labels = as.character(0:n_mut)) +
    scale_y_continuous(name = 'Function',
                       expand = c(0.05, 0.05)) +
    scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = 'none',
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)
  
  return(myplot)

}

predictF <- function(df, target, infer_residuals = FALSE, function_cols = 'auto') {
  
  if (infer_residuals) {
    
    # apply closure condition to the WHOLE landscape to infer values of residuals
    # only advisable for <10 species (otherwise pseudoinverse computation explodes)
    
    ### --------------------------------------------------------------------------
    # aux. functions
    
    nSpecies <- function(community) {
      
      # returns the number of species in a community (given in 'sp1,sp2,sp3,...' format)
      # returns 0 if input is NA
      
      as.numeric(sapply(community,
                        FUN = function(comm) {
                          if (is.na(comm)) return(0)
                          else return(length(strsplit(comm, split = ',')[[1]]))
                        }))
      
    }
    
    orderName <- function(community) {
      
      # order community names with the form 'sp2,sp2,sp3,...' so that species consistently
      # appear in alphabetical order
      
      as.character(sapply(community,
                          FUN = function(s) paste(sort(strsplit(s, split = ',')[[1]]),
                                                  collapse = ',')))
    }
    
    string2matrix <- function(data) {
      
      # takes a data frame with 2 columns, the first one containing community names
      # (species present separated by commas with no spaces, e.g. 'sp1,sp2,sp3,...') and
      # the second one community functions, and converts it to binary matrix format
      
      species <- sort(unique(unlist(lapply(data[, 1],
                                           FUN = function(community) strsplit(community, split = ',')[[1]]))))
      
      data_out <- lapply(1:nrow(data),
                         FUN = function(i) matrix(as.numeric(species %in% strsplit(data[i, 1], split = ',')[[1]]), nrow = 1))
      data_out <- do.call(rbind, data_out)
      colnames(data_out) <- species
      
      data_out <- cbind(data_out, fun = data[, 2])
      
      return(data_out)
      
    }
    
    matrix2string <- function(data) {
      
      # inverse operation with respect to string2matrix
      
      species <- colnames(data)[-ncol(data)]
      
      communities <- sapply(1:nrow(data),
                            FUN = function(i) paste(species[data[i, -ncol(data)] == 1], collapse = ','))
      communities <- orderName(communities)
      
      return(data.frame(community = communities,
                        fun = data[, ncol(data)]))
      
    }
    
    containsSpecies <- function(species, community) {
      
      # returns TRUE if species is present in community, FALSE otherwise
      
      return(sapply(community,
                    FUN = function(comm) species %in% strsplit(comm, split = ',')[[1]]))
      
    }
    
    closestPaths <- function(target_comm, comm_list, species.removal = F) {
      
      # given a focal community (comm) and a list of other communities (comm_list), returns those in comm_list that are closest to comm (in terms of species addition)
      
      # fetch communities that are potential ancestors of the target
      isAncestor <- function(community_1, community_2) { # check if community_1 is an ancestor of community_2 (i.e. all species in community_1 are also in community_2)
        sapply(community_1,
               FUN = function(x) all(strsplit(x, split = ',')[[1]] %in% strsplit(community_2, split = ',')[[1]]))
      }
      if (!species.removal) comm_list <- comm_list[isAncestor(comm_list, target_comm)]
      
      if (length(comm_list) == 0) {
        
        warning(paste('For community\n', target_comm, '\nno ancestors in data. Returning <empty> community.', sep = ''))
        
        comm_list <- ''
        
      }
      
      all_comms <- data.frame(community = c(target_comm, comm_list),
                              fun = NA)
      all_comms <- string2matrix(all_comms)
      
      comm <- all_comms[1, 1:(ncol(all_comms) - 1), drop = F]
      comm_list <- all_comms[2:nrow(all_comms), 1:(ncol(all_comms) - 1), drop = F]
      
      dist <- sapply(1:nrow(comm_list),
                     FUN = function(i) sum(abs(comm_list[i, ] - comm)))
      
      which_min <- which(dist == min(dist))
      closest_comms <- comm_list[which_min, , drop = F]
      
      
      # paths to comm
      paths <- lapply(1:nrow(closest_comms),
                      FUN = function(i) comm - closest_comms[i, ])
      paths <- do.call(rbind, paths)
      
      sources <- matrix2string(cbind(as.data.frame(closest_comms), fun = NA))$community
      
      # return info
      return(cbind(source = sources, target = target_comm, dist = dist[which_min], as.data.frame(paths)))
      
    }
    
    pathSteps <- function(target, source, single.traj = F) {
      
      # get steps in trajectory from source to target communities through addition/removal of species
      # is the closure condition is satisfied, a single trajectory is sufficient (speeds up calculations)
      
      path <- closestPaths(target, source)
      
      traj <- path[1, 4:ncol(path), drop = F]
      traj <- colnames(traj)[abs(traj[1, ]) == 1]
      
      # should a single trajectory be returned? (this is enough if the closer condition is satisfied)
      if (single.traj) {
        
        traj <- matrix(traj, nrow = 1)
        
      } else { # otherwise, return all possible trajectories (orders of species addition)
        
        if (length(traj) <= 4) {
          
          traj <- do.call(rbind, permn(traj))
          
        } else {
          
          # if number of steps is too large, getting all possible trajectories is computationally out of reach
          # in that case, we take only 25 arbitrarily chosen paths
          traj <- do.call(rbind,
                          lapply(1:25,
                                 FUN = function(i) matrix(sample(traj, size = length(traj), replace = F), nrow = 1)))
          
        }
        
      }
      
      traj_steps <- lapply(1:nrow(traj),
                           FUN = function(i) {
                             
                             # trajectory
                             cumtraj <- orderName(sapply(1:length(traj[i, ]), FUN = function(j) paste(traj[i, 1:j], collapse = ',')))
                             cumtraj <- c(source,
                                          sapply(1:length(cumtraj), FUN = function(j) orderName(paste(c(strsplit(source, split = ',')[[1]],
                                                                                                        cumtraj[j]), collapse = ','))))
                             cumtraj <- sapply(cumtraj, # if there are repetitions of species, that means the species is being removed instead of added
                                               FUN = function(x) {
                                                 
                                                 xtab <- table(strsplit(x, split = ',')[[1]])
                                                 removed_sp <- names(xtab)[xtab > 1]
                                                 
                                                 xout <- strsplit(x, split = ',')[[1]]
                                                 xout <- paste(xout[!(xout %in% removed_sp)], collapse = ',')
                                                 
                                                 return(xout)
                                                 
                                               })
                             cumtraj <- as.character(cumtraj)
                             
                             # knock-ins at each step
                             knockins <- traj[i, ]
                             
                             # backgrounds
                             backgrounds <- cumtraj[-length(cumtraj)]
                             backgrounds <- sapply(1:length(backgrounds),
                                                   FUN = function(i) {
                                                     xout <- strsplit(backgrounds[i], split = ',')[[1]]
                                                     xout <- xout[xout != knockins[i]]
                                                     xout <- paste(xout, collapse = ',')
                                                     return(xout)
                                                   })
                             backgrounds <- as.character(backgrounds)
                             
                             # signs
                             signs <- sapply(knockins, FUN = function(x) path[1, x])
                             signs <- as.numeric(signs)
                             
                             return(list(trajectory = cumtraj,
                                         backgrounds = backgrounds,
                                         knock_ins = knockins,
                                         signs = signs))
                             
                           })
      
      if(single.traj) traj_steps <- traj_steps[[1]]
      
      return(traj_steps)
      
    }
    
    # end of aux. functions
    ### --------------------------------------------------------------------------
    
    
    
    ### --------------------------------------------------------------------------
    # preprocessing
    
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
    
    # mean and SD of function across replicates
    df <- do.call(data.frame,
                  aggregate(fun ~ .,
                            data = do.call(rbind, df),
                            FUN = function(x) c(mean = mean(x), sd = sd(x))))
    df <- df[, c(elements, 'fun.mean')]
    
    # get FEEs
    fits <- getFEEs(df)
    elements <- fits$focal_element
    fits <- data.frame(a = fits$intercept, b = fits$slope)
    rownames(fits) <- elements
    
    # add "empty" assemblage if necessary
    data <- matrix2string(df)
    if (!any(data$community == '')) data <- rbind(data,
                                                  data.frame(community = '',
                                                             fun = 0))
    
    # format targets
    if(!is.data.frame(target)) {
      if (is.matrix(target)) {
        target <- as.data.frame(target)
      } else {
        target <- as.data.frame(t(target))
        colnames(target) <- elements
      }
    }
    
    target_id <- matrix2string(cbind(target, fun = NA))$community
    
    # end of preprocessing
    ### --------------------------------------------------------------------------
    
    # make GE data
    ge_data <- structureData(df)
    ge_data <- data.frame(background = matrix2string(ge_data[, elements])$community,
                          knock_in = ge_data$focal_element,
                          background_f = ge_data$fun_background.mean,
                          d_f = ge_data$delta_fun.mean)
    
    # get linear fits
    # fits <- makeFEEs(ge_data)
    
    # fetch species names and build all combinatorial arrangements
    species <- sort(unique(ge_data$knock_in))
    communities <- orderName(unlist(sapply(1:length(species),
                                           FUN = function(i) sapply(combn(species,
                                                                          i,
                                                                          simplify = FALSE),
                                                                    FUN = paste, collapse = ','))))
    communities <- c('<empty>', communities) # attach the 'empty' community
    
    # list epsilons that are known from observations
    ge_data$background[ge_data$background == ''] <- '<empty>'
    known_eps <- setNames(ge_data$d_f- (fits[ge_data$knock_in, 'a'] + fits[ge_data$knock_in, 'b']*ge_data$background_f),
                          paste(ge_data$background, ge_data$knock_in, sep = '+'))
    
    # estimate standard deviations of epsilons
    sigma <- sapply(species,
                    FUN = function(sp) sd(known_eps[grepl(paste('\\+', sp, sep = ''), names(known_eps))]))
    
    # arrange epsilons in matrix form (rows: backgrounds, columns: knock-ins), this makes them easier to access later on
    eps_matrix <- matrix(NA, nrow = length(communities), ncol = length(species))
    colnames(eps_matrix) <- species
    rownames(eps_matrix) <- communities
    
    # some of the epsilons are known from the observations: add them to the matrix
    for (i in 1:length(known_eps)) {
      eps_matrix[strsplit(names(known_eps)[i], split = '\\+')[[1]][1],
                 strsplit(names(known_eps)[i], split = '\\+')[[1]][2]] <- known_eps[i]
    }
    
    # elements of the matrix where the knock-in is already present in the background are set to NaN, other elements are either numeric (known epsilons) or NA (incognitas)
    for (sp in species) {
      eps_matrix[containsSpecies(sp, rownames(eps_matrix)), sp] <- NaN
    }
    
    # epsilons in vector form (makes it easier to operate with them)
    eps <- as.data.frame(matrix(NA, nrow = length(eps_matrix), ncol = 2))
    for (i in 1:nrow(eps_matrix)) {
      for (j in 1:ncol(eps_matrix)) {
        eps[ncol(eps_matrix)*(i-1) + j, 1] <- paste(rownames(eps_matrix)[i], colnames(eps_matrix)[j], sep = '+')
        eps[ncol(eps_matrix)*(i-1) + j, 2] <- eps_matrix[i, j]
      }
    }
    eps <- eps[!is.nan(eps[, 2]), ]
    eps <- eps[order(eps[, 2], na.last = FALSE), ]
    rownames(eps) <- eps[, 1]
    eps <- matrix(eps[, 2], dimnames = list(eps[, 1], 'eps'))
    
    rownames(eps) <- gsub('<empty>', '', rownames(eps))
    communities[communities == '<empty>'] <- ''
    
    # build M matrix
    # i and j are two arbitrary species
    # B is an arbitrary background (not containing i or j), Bi and Bj are the same backgrounds containing also i and j respectively
    # a_i, b_i are the intercept and slope of the fit FEE for species i
    # e_i(B) is the epsilon for species i on background B
    # general form of the constraint:
    #   a_i*b_j + (1 + b_j)*e_i(B) + e_j(Bi) = a_j*b_i + (1 + b_i)*e_j(B) + e_i(Bj)
    constraints <- data.frame(B = character(0), i = character(0), j = character(0))
    
    # list all backgrounds with at least 2 species less than the most rich community (containing all species)
    bgs_valid <- communities[nSpecies(communities) <= length(species) - 2]
    
    # for each of those backgrounds, get all possible pairs of additional species (not in the background itself)
    sp_valid <- containsSpecies(species, bgs_valid)
    for (i in 1:length(bgs_valid)) {
      sp <- species[!sp_valid[, i]]
      pairs <- t(combn(sp, 2))
      constraints <- rbind(constraints,
                           data.frame(B = as.character(bgs_valid[i]),
                                      i = pairs[, 1],
                                      j = pairs[, 2]))
    }
    
    # some of the constraints are redundant; to identify them we need to check which edges of the 'community graph' are covered by each constraint
    isDescendant <- function(community_1, community_2) { # check if community_2 is a descendant of community_1 (i.e. all species in community_1 are also in community_2)
      sapply(community_2,
             FUN = function(x) all(strsplit(community_1, split = ',')[[1]] %in% strsplit(x, split = ',')[[1]]))
    }
    
    edges <- NULL
    nSpecies_comm <- nSpecies(communities)
    for (i in 1:length(communities)) {
      edges <- c(edges,
                 paste(communities[i],
                       communities[nSpecies_comm == (1 + nSpecies_comm[i]) & isDescendant(communities[i], communities)],
                       sep = ' / '))
    }
    edges <- data.frame(edge = edges,
                        covered = FALSE)
    
    redundant_constraints <- NULL
    for (i in 1:nrow(constraints)) {
      
      # each constraint covers 4 communities
      c1 <- constraints[i, 1]
      c2 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 2]), collapse = ','))
      c3 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 3]), collapse = ','))
      c4 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 2], constraints[i, 3]), collapse = ','))
      
      # each constraint covers 4 edges
      edges_i <- c(paste(c1, c2, sep = ' / '),
                   paste(c1, c3, sep = ' / '),
                   paste(c2, c4, sep = ' / '),
                   paste(c3, c4, sep = ' / '))
      
      # which of those are already covered?
      covered_i <- edges$covered[edges$edge %in% edges_i]
      
      if (all(covered_i)) { # if they are all covered, tag the i-th constraint as redundant
        redundant_constraints <- c(redundant_constraints, i)
      } else { # if they are not all covered, keep the constraint and set them all to covered
        edges$covered[edges$edge %in% edges_i] <- TRUE
      }
      
    }
    
    # remove redundant constraints
    constraints <- constraints[-redundant_constraints, ]
    
    # formulate each constraint in matrix form
    M <- matrix(0, nrow = nrow(constraints), ncol = nrow(eps))
    C <- matrix(0, nrow = nrow(constraints), ncol = 1)
    rownames(M) <- paste('C.', 1:nrow(constraints), sep = '')
    rownames(C) <- rownames(M)
    colnames(M) <- rownames(eps)
    
    for (k in 1:nrow(M)) {
      
      i <- constraints$i[k]
      j <- constraints$j[k]
      B <- constraints$B[k]
      
      Bi <- orderName(paste(c(strsplit(B, split = ',')[[1]], i), collapse = ','))
      Bj <- orderName(paste(c(strsplit(B, split = ',')[[1]], j), collapse = ','))
      
      a_i <- fits[i, 'a']
      b_i <- fits[i, 'b']
      a_j <- fits[j, 'a']
      b_j <- fits[j, 'b']
      
      e_i_B <- paste(B, i, sep = '+')
      e_j_B <- paste(B, j, sep = '+')
      e_j_Bi <- paste(Bi, j, sep = '+')
      e_i_Bj <- paste(Bj, i, sep = '+')
      
      M[k, e_i_B] <- 1 + b_j
      M[k, e_j_B] <- -(1 + b_i)
      M[k, e_j_Bi] <- 1
      M[k, e_i_Bj] <- -1
      
      C[k] <- a_j*b_i - a_i*b_j
      
    }
    
    # check if any constraint is automatically satisfied by the known epsilons
    check <- M %*% eps - C
    check <- sapply(1:nrow(check), FUN = function(i) check[i] < 1e-5) # tolerance: 1e-5 (elements won't be exactly 0 even when the constraint is satisfied due to floating-point precision)
    
    stopifnot(all(check[!is.na(check)] == TRUE)) # all checks should either be TRUE (constraint satisfied) or NA (constraint not evaluated due to missing epsilons), if this is not the case something has failed
    
    # the dimension of M can be further reduced by moving the known epsilons to the right-hand side of the equation M %*% eps = C
    eps_1 <- eps
    eps_1[!is.na(eps)] <- 0
    eps_2 <- eps
    eps_2[is.na(eps_2)] <- 0 # eps = eps_1 + eps_2 ---> M %*% eps_1 = C - M %*% eps_2
    
    C_prime <- C - (M %*% eps_2) # this makes it so the equation to solve is M %*% eps_1 = C_prime
    remove_these <- sapply(1:nrow(M), FUN = function(i) all(as.numeric(M[i, ]) == rep(0, ncol(M)))) # to remove constraints that are automatically satisfied by the known epsilons
    C_prime <- C_prime[!remove_these, ]
    M <- M[!remove_these, 1:sum(is.na(eps))] # also removes the columns corresponding to known epsilons
    
    # scale elements of M by the sigmas
    M_sigma <- M * matrix(rep(sigma[gsub('.*\\+' ,'', colnames(M))], nrow(M)), nrow = nrow(M), byrow = TRUE)
    
    # now we just need to solve the system of linear equations M_sigma %*% eps_sigma = C_prime where the elements of eps_sigma are the incognitas
    # use the pseudoinverse for the solution that minimizes sum of squares of the elements of eps_sigma
    M_sigma_pseudoinv <- ginv(M_sigma)
    eps_sigma_0 <- M_sigma_pseudoinv %*% C_prime
    rownames(eps_sigma_0) <- colnames(M)
    
    # undo the sigma transformation to recover the epsilons
    eps_1 <- eps_sigma_0 * sigma[gsub('.*\\+' ,'', rownames(eps_sigma_0))]
    
    # add these values to the matrix of epsilons
    for (i in 1:nrow(eps_1)) {
      B <- gsub('\\+.*', '', rownames(eps)[i])
      if (B == '') B <- '<empty>'
      knockin <- gsub('.*\\+', '', rownames(eps)[i])
      eps_matrix[B, knockin] <- eps_1[i]
    }
    
    
    
    
    
    
    
    
    ### MAKE PREDICTIONS
    
    # iterate over target (multiple target communities can be passed at once)
    predicted_f <- sapply(target_id,
                          FUN =  function(t) {
                            
                            # get closest path to target community (although any path should be equivalent if the closure condition is satisfied globally)
                            path <- closestPaths(t, data$community)
                            steps <- pathSteps(t, path$source[1], single.traj = T)
                            steps$backgrounds[steps$backgrounds == ''] <- '<empty>'
                            
                            # iterate
                            comm <- steps$trajectory[1]
                            fun <- mean(data$fun[data$community == comm])
                            for (s in 1:length(steps$knock_ins)) {
                              comm <- steps$trajectory[s]
                              fun <- fun + steps$signs[s]*(fits[steps$knock_ins[s], 'a'] + fits[steps$knock_ins[s], 'b']*fun + eps_matrix[steps$backgrounds[s], steps$knock_ins[s]])
                            }
                            
                            return(fun)
                            
                          })
    
    out <- cbind(target, predicted_fun = predicted_f)
    rownames(out) <- NULL
    
    return(out)
    
  } else {
    
    # predict the function of a community (target) starting from the function of an in-sample community
    # considers every unknown residual to be zero
    # first obtains the in-sample communities that are closer to the target community, then every possible path (i.e. order of addition/removal of species)
    # final prediction is the average across all paths
    
    # ----------------------------------------------------------------------------
    # internal, auxiliary functions
    
    orderName <- function(community) {
      
      # order community names with the form 'sp2,sp2,sp3,...' so that species consistently
      # appear in alphabetical order
      
      as.character(sapply(community,
                          FUN = function(s) paste(sort(strsplit(s, split = ',')[[1]]),
                                                  collapse = ',')))
    }
    
    string2matrix <- function(data) {
      
      # takes a data frame with 2 columns, the first one containing community names
      # (species present separated by commas with no spaces, e.g. 'sp1,sp2,sp3,...') and
      # the second one community functions, and converts it to binary matrix format
      
      species <- sort(unique(unlist(lapply(data[, 1],
                                           FUN = function(community) strsplit(community, split = ',')[[1]]))))
      
      data_out <- lapply(1:nrow(data),
                         FUN = function(i) matrix(as.numeric(species %in% strsplit(data[i, 1], split = ',')[[1]]), nrow = 1))
      data_out <- do.call(rbind, data_out)
      colnames(data_out) <- species
      
      data_out <- cbind(data_out, fun = data[, 2])
      
      return(data_out)
      
    }
    
    matrix2string <- function(data) {
      
      # inverse operation with respect to string2matrix
      
      species <- colnames(data)[-ncol(data)]
      
      communities <- sapply(1:nrow(data),
                            FUN = function(i) paste(species[data[i, -ncol(data)] == 1], collapse = ','))
      communities <- orderName(communities)
      
      return(data.frame(community = communities,
                        fun = data[, ncol(data)]))
      
    }
    
    containsSpecies <- function(species, community) {
      
      # returns TRUE if species is present in community, FALSE otherwise
      
      return(sapply(community,
                    FUN = function(comm) species %in% strsplit(comm, split = ',')[[1]]))
      
    }
    
    closestPaths <- function(target_comm, comm_list, species.removal = F) {
      
      # given a focal community (comm) and a list of other communities (comm_list), returns those in comm_list that are closest to comm (in terms of species addition)
      
      # fetch communities that are potential ancestors of the target
      isAncestor <- function(community_1, community_2) { # check if community_1 is an ancestor of community_2 (i.e. all species in community_1 are also in community_2)
        sapply(community_1,
               FUN = function(x) all(strsplit(x, split = ',')[[1]] %in% strsplit(community_2, split = ',')[[1]]))
      }
      if (!species.removal) comm_list <- comm_list[isAncestor(comm_list, target_comm)]
      
      if (length(comm_list) == 0) {
        
        warning(paste('For community\n', target_comm, '\nno ancestors in data. Returning <empty> community.', sep = ''))
        
        comm_list <- ''
        
      }
      
      all_comms <- data.frame(community = c(target_comm, comm_list),
                              fun = NA)
      all_comms <- string2matrix(all_comms)
      
      comm <- all_comms[1, 1:(ncol(all_comms) - 1), drop = F]
      comm_list <- all_comms[2:nrow(all_comms), 1:(ncol(all_comms) - 1), drop = F]
      
      dist <- sapply(1:nrow(comm_list),
                     FUN = function(i) sum(abs(comm_list[i, ] - comm)))
      
      which_min <- which(dist == min(dist))
      closest_comms <- comm_list[which_min, , drop = F]
      
      
      # paths to comm
      paths <- lapply(1:nrow(closest_comms),
                      FUN = function(i) comm - closest_comms[i, ])
      paths <- do.call(rbind, paths)
      
      sources <- matrix2string(cbind(as.data.frame(closest_comms), fun = NA))$community
      
      # return info
      return(cbind(source = sources, target = target_comm, dist = dist[which_min], as.data.frame(paths)))
      
    }
    
    pathSteps <- function(target, source, single.traj = F) {
      
      # get steps in trajectory from source to target communities through addition/removal of species
      # is the closure condition is satisfied, a single trajectory is sufficient (speeds up calculations)
      
      path <- closestPaths(target, source)
      
      traj <- path[1, 4:ncol(path), drop = F]
      traj <- colnames(traj)[abs(traj[1, ]) == 1]
      
      # should a single trajectory be returned? (this is enough if the closer condition is satisfied)
      if (single.traj) {
        
        traj <- matrix(traj, nrow = 1)
        
      } else { # otherwise, return all possible trajectories (orders of species addition)
        
        if (length(traj) <= 4) {
          
          traj <- do.call(rbind, permn(traj))
          
        } else {
          
          # if number of steps is too large, getting all possible trajectories is computationally out of reach
          # in that case, we take only 25 arbitrarily chosen paths
          traj <- do.call(rbind,
                          lapply(1:25,
                                 FUN = function(i) matrix(sample(traj, size = length(traj), replace = F), nrow = 1)))
          
        }
        
      }
      
      traj_steps <- lapply(1:nrow(traj),
                           FUN = function(i) {
                             
                             # trajectory
                             cumtraj <- orderName(sapply(1:length(traj[i, ]), FUN = function(j) paste(traj[i, 1:j], collapse = ',')))
                             cumtraj <- c(source,
                                          sapply(1:length(cumtraj), FUN = function(j) orderName(paste(c(strsplit(source, split = ',')[[1]],
                                                                                                        cumtraj[j]), collapse = ','))))
                             cumtraj <- sapply(cumtraj, # if there are repetitions of species, that means the species is being removed instead of added
                                               FUN = function(x) {
                                                 
                                                 xtab <- table(strsplit(x, split = ',')[[1]])
                                                 removed_sp <- names(xtab)[xtab > 1]
                                                 
                                                 xout <- strsplit(x, split = ',')[[1]]
                                                 xout <- paste(xout[!(xout %in% removed_sp)], collapse = ',')
                                                 
                                                 return(xout)
                                                 
                                               })
                             cumtraj <- as.character(cumtraj)
                             
                             # knock-ins at each step
                             knockins <- traj[i, ]
                             
                             # backgrounds
                             backgrounds <- cumtraj[-length(cumtraj)]
                             backgrounds <- sapply(1:length(backgrounds),
                                                   FUN = function(i) {
                                                     xout <- strsplit(backgrounds[i], split = ',')[[1]]
                                                     xout <- xout[xout != knockins[i]]
                                                     xout <- paste(xout, collapse = ',')
                                                     return(xout)
                                                   })
                             backgrounds <- as.character(backgrounds)
                             
                             # signs
                             signs <- sapply(knockins, FUN = function(x) path[1, x])
                             signs <- as.numeric(signs)
                             
                             return(list(trajectory = cumtraj,
                                         backgrounds = backgrounds,
                                         knock_ins = knockins,
                                         signs = signs))
                             
                           })
      
      if(single.traj) traj_steps <- traj_steps[[1]]
      
      return(traj_steps)
      
    }
    
    # end of internal, auxiliary functions
    # ----------------------------------------------------------------------------
    
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
    
    # mean and SD of function across replicates
    df <- do.call(data.frame,
                  aggregate(fun ~ .,
                            data = do.call(rbind, df),
                            FUN = function(x) c(mean = mean(x), sd = sd(x))))
    df <- df[, c(elements, 'fun.mean')]
    
    # get FEEs
    fits <- getFEEs(df)
    elements <- fits$focal_element
    fits <- data.frame(a = fits$intercept, b = fits$slope)
    rownames(fits) <- elements
    
    # add "empty" assemblage if necessary
    data <- matrix2string(df)
    if (!any(data$community == '')) data <- rbind(data,
                                                  data.frame(community = '',
                                                             fun = 0))
    
    #format targets
    if(!is.data.frame(target)) {
      if (is.matrix(target)) {
        target <- as.data.frame(target)
      } else {
        target <- as.data.frame(t(target))
        colnames(target) <- elements
      }
    }
    
    target_id <- matrix2string(cbind(target, fun = NA))$community
    
    fun <- sapply(target_id,
                  FUN = function(t) {
                    
                    paths <- closestPaths(t, data$community)
                    
                    # if there are too many possible paths from an in-sample community to the target community,
                    # we randomly choose 5 of them to avoid heavily increasing computation time
                    if (nrow(paths) > 5) paths <- paths[order(sample(1:nrow(paths), size = 5, replace = F)), ]
                    
                    fun <- lapply(1:nrow(paths),
                                  FUN = function(i) {
                                    
                                    path <- paths[i, ]
                                    steps <- pathSteps(t, path$source[1], single.traj = F)
                                    
                                    fun_i <- sapply(steps,
                                                    FUN = function(st) {
                                                      
                                                      st$backgrounds[st$backgrounds == ''] <- '<empty>'
                                                      
                                                      # iterate
                                                      comm <- st$trajectory[1]
                                                      fun <- mean(data$fun[data$community == comm])
                                                      for (s in 1:length(st$knock_ins)) {
                                                        comm <- st$trajectory[s]
                                                        fun <- fun + st$signs[s]*(fits[st$knock_ins[s], 'a'] + fits[st$knock_ins[s], 'b']*fun)
                                                      }
                                                      
                                                      return(fun)
                                                      
                                                    })
                                    
                                    return(fun_i)
                                    
                                  })
                    
                    return(mean(unlist(fun)))
                    
                  })
    
    out <- cbind(target, predicted_fun = fun)
    rownames(out) <- NULL
    
    return(out)
    
  }
  
}

