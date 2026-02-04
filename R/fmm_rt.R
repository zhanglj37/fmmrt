#' Bayesian Factor Mixture Model with Response Time
#'
#' Fit a Bayesian factor mixture model that jointly models item responses
#' and response times to classify respondents into attentive and careless groups.
#' The model is estimated using JAGS and returns parameter estimates,
#' convergence diagnostics, and posterior classification probabilities.
#'
#' @param responses A numeric matrix of item responses with rows representing respondents
#'   and columns representing items.
#' @param rt A numeric matrix of response times with the same dimensions as \code{responses}.
#'   Response times are internally log transformed.
#' @param iter Total number of MCMC iterations.
#' @param burnin Number of burn in iterations.
#' @param save.file Logical indicating whether to save a CSV file
#'   containing classification results and original data.
#' @param seed Random seed for reproducibility.
#'
#' @details
#' The model assumes two latent classes representing attentive and careless respondents.
#' Item responses are modeled using a factor analytic structure for the attentive group
#' and a reduced structure for the careless group.
#' Response times are modeled jointly to improve classification.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item epsr Maximum potential scale reduction factor across monitored parameters.
#'   \item estimates Posterior summaries including means, quantiles, and HPD intervals.
#'   \item probability A data frame containing posterior probabilities of being attentive
#'   and binary group classification for each respondent.
#' }
#'
#' @references
#' Zhang, L., Ulitzsch, E., & Domingue, B. W. (2025).
#' Bayesian factor mixture modeling with response time for detecting careless respondents.
#' \emph{Behavior Research Methods}, 57(10), 286.
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc HPDinterval coda.options
#' @importFrom rjags load.module
#' @export


fmm_rt <- function(responses, rt, iter = 10000, burnin = 5000, save.file = TRUE, seed = 123){


  ############################
  # prepare the data
  ############################
  J = ncol(responses)
  N = nrow(responses)
  logrt = log(rt)
  jags.data <- list("y"=responses, 'logrt'=logrt, "N"=N, "J"=J)
  entities.to.monitor <- c("loading_y1", "mu_y1", "var_response1", "mu_t1", "var_speed1", "var_rt1", "mu_y2", "var_response2", "mu_t2", "var_rt2", "prob")  # parameters of interest

  ############################
  # prepare the model
  ############################
  modelstring <- as.character("
model {
    ############################
    # priors for classification
    ############################
    for(i in 1:N) {
      ind[i] ~ dcat(prob[1:2,i])
      prob[1:2,i] ~ ddirch(alpha[1:2])
    }
    for(k in 1:2) {alpha[k] <- 1}

    ############################
    # parameters of the model for theattentive group
    ############################
    for(i in 1:N){
      # factor scores
      tau[i] ~ dnorm(0, invsig_tau)
      omega[i] ~ dnorm(0, 1)
    }
    invsig_tau ~ dgamma(0.01, 0.01)
    var_speed1 <- 1/invsig_tau

    for(j in 1:J) {
      # loading
      loading_y1[j] ~ dnorm(0, 0.01)T(0,) # avoid swiching
      # intercept
      mu_t1[j] <- mu_t2 + delta[j]
      delta[j] ~ dnorm(0, 0.01)T(0,)
      mu_y1[j] ~ dnorm(0, 0.01)
      # variance and inverse-variance
      invsig_y1[j] ~ dgamma(0.01, 0.01)
      var_response1[j] <- 1/invsig_y1[j]
      invsig_t1[j] ~ dgamma(0.01, 0.01)
      var_rt1[j] <- 1/invsig_t1[j]
      invsig_y[j, 1] <- invsig_y1[j]
      invsig_t[j, 1] <- invsig_t1[j]
    }

    ############################
    # parameters of the model for the careless group
    ############################
    mu_y2 ~ dnorm(0, 0.01)
    mu_t2 ~ dnorm(0, 0.01)
    invsig_t2 ~ dgamma(0.01, 0.01)
    var_rt2 <- 1/invsig_t2
    invsig_y2 ~ dgamma(0.01, 0.01)
    var_response2 <- 1/invsig_y2
    for(j in 1:J){
      invsig_y[j, 2] <- invsig_y2
      invsig_t[j, 2] <- invsig_t2
    }

    ############################
    # Mixture Model
    ############################
    for(i in 1:N) {
      for(j in 1:J) {
        # Mixture model for response
        mu_y_ij[i, j, 1] <- mu_y1[j] + loading_y1[j] * omega[i]
        mu_y_ij[i, j, 2] <- mu_y2
        y[i, j] ~ dnorm(mu_y_ij[i, j, ind[i]], invsig_y[j, ind[i]])

        # Mixture model or logrt
        mu_t_ij[i, j, 1] <- mu_t1[j] - tau[i]
        mu_t_ij[i, j, 2] <- mu_t2
        logrt[i, j] ~ dnorm(mu_t_ij[i, j, ind[i]], invsig_t[j, ind[i]])
      }
    }
}

") # closes the model as string
  model.file.name <- "fmm_rt.txt"
  write(x=modelstring, file=model.file.name, append=FALSE)

  ############################
  # Run the model
  ############################
  load.module("glm")
  set.seed(seed)
  modelfit <- jags(model.file=model.file.name,
                   data=jags.data,
                   n.chains=2,
                   parameters.to.save=entities.to.monitor,
                   n.iter=iter,
                   n.burn=burnin,
                   n.thin=1,
                   progress.bar="none")

  ############################
  # Obtain the results
  ############################
  # EPSR
  max.psr = max(modelfit$BUGSoutput$summary[,'Rhat'])

  # Parameter Estimation
  modelfit.mcmc <- as.mcmc(modelfit)
  modelfit.mcmc.onelist <- as.mcmc(do.call(rbind,modelfit.mcmc)) # combine chains
  coda.options(combine.stats=TRUE, combine.plots=TRUE)
  summary.statistics <- summary(modelfit.mcmc.onelist)
  statistic = as.data.frame(summary.statistics$statistic)
  quantiles = as.data.frame(summary.statistics$quantiles)
  HPD.interval <- HPDinterval(modelfit.mcmc.onelist, prob=.95)
  estimates = cbind(statistic, quantiles, as.data.frame(HPD.interval))
  attentive_prob <- estimates[grep("prob\\[1,(\\d+)\\]", rownames(estimates)), ]
  estimates <- estimates[!grepl("prob", rownames(estimates)), ]

  # Classification (based on the median estimates of the group probability)
  probability <- data.frame(
    rownum = as.numeric(gsub("prob\\[1,(\\d+)\\]", "\\1", rownames(attentive_prob))), # Extract row number
    prob = attentive_prob[, "50%"], # Extract the median estimation
    group = ifelse(attentive_prob[, "50%"] > 0.5, 1, 0)
  )
  probability = probability[order(probability$rownum), ]
  rownames(probability) <- NULL


  if(save.file) write.csv(cbind(probability, responses, rt), file="dat_group.csv", row.names=F)
  if(file.exists(model.file.name)) file.remove(model.file.name)

  results = list(epsr = max.psr, estimates = estimates, probability = probability)
  return(results)
}


