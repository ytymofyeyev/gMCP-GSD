# derive analysis time since start using input dataset 'D' and 'enrollment'
deriveCalendarTime <- function(D, enrollment) {
  Ts <- vector("list", nrow(D))
  if (!"enrollment" %in% names(D))
    D <-
      D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))
  for (i in 1:nrow(D)) {
    Ts[[i]] <- sapply(D$iaSpec[[i]], function(e) {
      n2Time(
        D$endpointParam[[e$H]],
        n = e$atIF * D$hypN[e$H],
        D$enrollment[[e$H]],
        ratio = D$allocRatio[e$H]
      )
    })
  }
  return(Ts)
}

# derivation sample size from calendar time at IA
derive_Ns <- function(D, enrollment, doRounding = TRUE) {
  Ns <- vector("list", nrow(D))
  if (!"enrollment" %in% names(D))
    D <-
      D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))
  for (i in 1:nrow(D)) {
    Ns[[i]] <- Time2n(
      x = D$endpointParam[[i]],
      T = D$iaTime[[i]],
      enrollment = D$enrollment[[i]],
      ratio = D$allocRatio[[i]]
    )
  }
  if (doRounding)
    Ns <- lapply(Ns, function(x)
      round(x))
  return(Ns)
}

# derivation of information fractions    # currently not used
deriveIF <- function(D, enrollment, digits = 2) {
  IFs <- list(nrow(D))
  for (i in 1:nrow(D)) {
    IFs[[i]] <- round(
      Time2n(
        x = D$endpointParam[[i]],
        T = D$iaTime[[i]],
        enrollment = enrollment,
        ratio = D$allocRatio[[i]]
      ) /
        D$hypN[[i]],
      digits = digits
    )
  }
  return(IFs)
}
###############################################################################


# extract delta (hypothesis parameter on the natural scale)
getDelta <- function (x, ...) {
  UseMethod("getDelta", x)
}

getDelta.default <- function(x) {
  x$p1 - x$p2
}

getDelta.tte_exp <- function(x) {
  x$p1 / x$p2     # return hazard ratio
}

###############################################################################
# extract effect size info (hypothesis parameter on the natural scale)
getEffectSizeDetails <- function (x, ...) {
  UseMethod("getEffectSizeDetails", x)
}

getEffectSizeDetails.normal <- function(x, digits = 2) {
  fmtStr <- sprintf("%%.%df", digits)
  sprintf(fmtStr, round(x$p1 - x$p2, digits = digits))
}

getEffectSizeDetails.binomial <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}

getEffectSizeDetails.binomial_pooled <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}

getEffectSizeDetails.binomial_unpooled <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}
getEffectSizeDetails.tte_exp <- function(x) {
  sprintf("HR = %.2f (mCntl = %.1f mo)", x$p1 / x$p2,-log(1 / 2) / x$p2)
}
###############################################################################

# compute standardization factor, psi,  that link natural parameter to the effect size, theta 
# that non-centrality parameter is in the form theta*sqrt(n)= psi*delta*sgrt(n)
# e.g, for log-rank test, psi = sqrt(c*(1-c)), where c=1/(1+ratio), ratio is the allocation ratio
getStandardizingCoef <- function (x, ratio=1, ...) {
    UseMethod("getStandardizingCoef", x )
}

getStandardizingCoef.default <- function(x, ratio = 1){
    sqrt(ratio)/(1+ratio) * 1
}

getStandardizingCoef.binomial <- function(x, ratio = 1){
    # x$alpha and x$beta use used to calibrate standardization coefficient 
    # when mapping binomial rates to normal effect size using n.I 
    if (is.null(x$alpha) | is.null(x$beta)){
        x$alpha <- 0.025
        x$beta  <- 0.85
    }
    nAux <- nBinomial( p1=x$p1, p2=x$p2, ratio=ratio,alpha=x$alpha, beta=x$beta)
    thetaAux <- (qnorm(1-x$alpha) + qnorm(1-x$beta))/sqrt(nAux)
    return(thetaAux/(x$p1-x$p2))
}

# binomial_pooled gives somewhat conservative power
getStandardizingCoef.binomial_pooled <- function(x, ratio = 1){
    if(!is.array(x$pPooled)) pPooled <- x$pPooled
    if (!is.null(x$p1) & !is.null(x$p2)) pPooled <- 1/(1+ratio)*(x$p1+ratio*x$p2)
    sqrt(ratio)/(1+ratio) / sqrt( pPooled *(1-pPooled)) 
}

getStandardizingCoef.binomial_unpooled <- function(x, ratio = 1){
    1/sqrt(x$p1*(1-x$p1)*(1+ratio) + x$p2*(1-x$p2)*(1+ratio)/ratio) 
}

################################################################################
# combine events from control and treatment arms
eEvents_total <-
  function(hr = 1,
           ratio = 1,
           lambdaC = 1,
           eta = 0,
           gamma = 1,
           R = 1,
           S = NULL,
           T,
           Tfinal = NULL,
           minfup = 0,
           digits = 4,
           target = 0)
# assuming a common dropout
  {
    if (T == 0)
      return(0)
    Qe <- ratio / (1 + ratio)
    eDC <-
      eEvents(
        lambda = lambdaC,
        eta = eta,
        gamma = gamma * (1 - Qe),
        R = R,
        S = S,
        T = T,
        Tfinal = Tfinal,
        minfup = minfup
      )
    eDE <-
      eEvents(
        lambda = lambdaC * hr,
        eta = eta,
        gamma = gamma * Qe,
        R = R,
        S = S,
        T = T,
        Tfinal = Tfinal,
        minfup = minfup
      )
    return(sum(eDC$d + eDE$d) - target)
  }
eEvents_totalVec <- Vectorize(eEvents_total, c("hr", "T"))

# get time when a given number of events will be reached
tEvents <-
  function (n,
            hr = 1,
            ratio = 1,
            lambdaC = 1,
            eta = 0,
            gamma = 1,
            R = 1,
            S = NULL,
            T,
            Tfinal = NULL,
            minfup = 0,
            tol = .Machine$double.eps ^ 0.25)
# n is the target number of events
  {
    # solve eEvents_total for parameter T
    z <-
      stats::uniroot(
        f = eEvents_total,
        interval = c(1e-04, 1e5),
        hr = hr,
        ratio = ratio,
        lambdaC = lambdaC,
        eta = eta,
        gamma = gamma,
        R = R,
        S = S,
        Tfinal = Tfinal,
        minfup = minfup,
        target = n,
        tol = tol
      )
    z$root
  }
tEventsVec <- Vectorize(tEvents, c("hr", "n"))

# Given sample size (or number events) obtain calendar time
n2Time <-
  function (x, n, enrollment, ratio) {
    UseMethod("n2Time", x)
  }

n2Time.default <- function(x, n, enrollment, ratio = 1) {
  fun <- function(t) {
    Time2n(x, t, enrollment = enrollment, ratio = ratio) - n
  }
  if (is.null(x$maturityTime))
    x$maturityTime <- 0
  # accrualDuration <- group_by(enrollment, stratum) %>%
  #     dplyr::mutate( totDuration = sum(duration))  %>% ungroup() %>%
  #     dplyr::select(totDuration) %>% max()
  accrualDuration <- sum(enrollment$duration)
  uniroot(fun, c(0, accrualDuration + x$maturityTime))$root
}

n2Time.tte_exp <- function(x, n, enrollment, ratio = 1) {
  tEvents(
    n = n,
    hr = x$p1 / x$p2,
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration
  )
}

# Given calendar time obtain sample size (or number events) available
Time2n <-
  function (x, T, enrollment, ratio) {
    UseMethod("Time2n", x)
  }
# x is an endpoint object

Time2n.default <- function(x, T, enrollment, ratio = 1) {
  if (is.null(x$dropoutHazard))
    eta <- 0
  else
    eta <- x$dropoutHazard
  if (is.null(x$maturityTime))
    maturatyTime <- 0
  else
    maturatyTime <- x$maturityTime
  
  timeAux <- c(0, cumsum(enrollment$duration))
  N_rand <- cumsum(c(0, enrollment$rate * c(enrollment$duration)))
  
  approx(
    x = timeAux + maturatyTime,
    y = N_rand * exp(as.numeric(-eta) * (maturatyTime)),
    rule = 2,
    xout = T
  )$y
  
}

Time2n.tte_exp <- function(x, T, enrollment, ratio = 1) {
  eEvents_totalVec(
    T = T,
    hr = x$p1 / x$p2,
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration
  )
}

################################################################################

plot_iaTiming <-
  function(D,
           enrollment,
           Tmax = 80,
           plotInPercent = TRUE) {
    if (!"enrollment" %in% names(D))
      D <-
        D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))
    
    
    timeGrid <- seq(0, Tmax, 1)
    Nmax <- sum(enrollment$rate * enrollment$duration)
    N_Rand <- approx(
      x = c(0, cumsum(enrollment$duration)),
      y = c(0, cumsum(enrollment$rate * enrollment$duration)),
      xout = timeGrid,
      rule = 2
    )$y
    Randomized <- if (plotInPercent)
      N_Rand / Nmax * 100
    else
      N_Rand
    
    Y <- list(nrow(D))
    for (i in 1:nrow(D)) {
      N <-
        Time2n(
          x = D$endpointParam[[i]],
          T = timeGrid,
          enrollment = D$enrollment[[i]],
          ratio = D$allocRatio[[i]]
        )
      Y[[i]] <- if (plotInPercent)
        N / D$hypN[i] * 100
      else
        N
    }
    
    
    maturity <- sapply(D$endpointParam,
                       function(x)
                         ifelse(
                           !is.null(x$maturityTime),
                           paste0(", ", x$maturityTime, " months data"),
                           ""
                         ))
    # names(Y) <- paste0( D$ep, " (",D$id,")", maturity) # to be a plot legend
    names(Y) <-
      paste0(D$id, ": ", D$tag, " ", maturity) # to be a plot legend
    dat <- tibble(
      Time = timeGrid,
      Randomized,
      as.data.frame(Y, check.names = FALSE),
      .name_repair = "minimal"
    ) %>%
      pivot_longer(!Time, names_to = "Type", values_to = "Y") %>%
      mutate(Type = factor(Type, levels = c(
        "Randomized",
        setdiff(unique(Type), "Randomized")
      )))
    
    if (plotInPercent)
      dat <- dat %>% dplyr::filter(Y <= 105)
    
    ylabStr <-
      ifelse(plotInPercent,
             "Statistical Information, in %",
             "Statistical Information")
    iaTimes <- unique(unlist(D$iaTime))
    ans <-
      ggplot(data = dat, aes(x = Time, y = Y, col = Type)) + geom_line() +
      geom_vline(xintercept = iaTimes, linetype = "dashed") +
      labs(x = "Time, in months",
           y = ylabStr) +
      scale_x_continuous(breaks = seq(0, Tmax, 6), limits = c(0, Tmax)) +
      annotate(
        "text",
        x = iaTimes,
        y = 0,
        label = round(iaTimes, 0),
        angle = 45,
        size = 3.5,
        hjust = 0
      )
    if (plotInPercent)
      ans <- ans + scale_y_continuous(breaks = seq(0, 100, 10))
    
    ans
  }



print_sfInfo <- function(x, digits = 4)
{
  if (!is.null(x)) {
    if (is.character(x$sfu)) {
      sfName <- case_match(
        x$sfu,
        "OF" ~ "O'Brien - Fleming",
        "WT" ~ " Wang - Tsiatis",
        "Pocock" ~ "Pocock",
        .default  = "mis-specified"
      )
      if (sfName == "mis-specified")
        stop("Mis-specified boundary")
    } else
      sfName <-
        do.call(x$sfu, list(0.025, 1, x$sfupar))$name # dummy call
    
    if (!is.null(x$sfupar) & length(x$sfupar) <= 4) {
      sfParam <- paste0(", parameter = ",
                        paste(lapply(x$sfupar, function(y) {
                          if (is.numeric(y))
                            round(y, digits)
                        }), collapse = " "))
    } else
      sfParam <- ""
    paste0(sfName, sfParam)
  } else {
    "No group sequential testing"
  }
}

# given graph G, for each hypothesis generate possible weights and list what
# rejections lead to that weight
getPossibleWeightsInfo <- function(G, # G is gMCP graph object
                                   numDigitsToRound = 5) {
  # aux function to form a character string of reject Hj
  formRejStr <- function (x) {
    aux <- sapply(x, function(y) {
      if (all(y == 1))
        "Initial allocation"
      else
        paste0("Successful ",
               paste("H",
                     which(y == 0),
                     collapse = " & ",
                     sep = ""))
    })
    paste(aux, collapse = " or ", sep = "_")
  }  # formRejStr
  
  K <- length(getNodes(G)) # number of hypotheses in G
  W <- generateWeights(G)
  # For each Hj generate possible weights and info of rejection scenarios that lead to that weight
  res <-
    tibble::tibble(Hint = apply(W[, 1:K], 1, c, simplify = FALSE),
                   as.data.frame(W[,-c(1:K)])) %>%
    mutate(# number of elementary hyp in 'Hint'
      mj = map_dbl(Hint, sum)) %>%
    pivot_longer(-c(Hint, mj), names_to = "Hj", values_to = "possibleWeight") %>%
    mutate(possibleWeight = round(possibleWeight, numDigitsToRound)) %>%
    # dplyr::filter(possibleWeight>0) %>%  ## ?? YT Sep. 2023
    group_by(Hj, possibleWeight) %>%
    reframe(max_mj = max(mj),
            mj = mj,
            Hint = Hint) %>%
    group_by(Hj, possibleWeight) %>%
    dplyr::filter(mj == max_mj) %>%
    reframe(listHj = list(Hint)) %>%
    group_by(Hj) %>%
    mutate(rejectedHypInfo = map_chr(listHj, formRejStr)) %>%
    dplyr::select(Hj, possibleWeight, rejectedHypInfo) %>%
    mutate(Hj = as.numeric(gsub(".*?([0-9]+).*", "\\1", Hj)))
  return(res)
}
# aux function that put spending function info from list to char string
report_MT_grSeq <- function( 
    G,        # graph for H1, ..., Hm
    D,
    # infoFr,   # list of length m setting info fraction for the Hi
    # spendFun,  # list of spending functions for Hi
    sigLevel  = 0.025, # overall significance level usually (1-sided)
    pdigits   = 5,     # p-value digits
    ddigits   = 3,     # delta digits
    idigits   = 3,     # info fraction digits
    powdigits = 2      # power digits
)
{
  nodes <- getNodes(G)
  m <- length(nodes)
  
  # possibleWeight <- apply(generateWeights(G)[,-c(1:m)],2,unique,simplify = FALSE)
  wInfo <- getPossibleWeightsInfo(G)
  possibleWeight <- split(wInfo$possibleWeight, wInfo$Hj)
  scenarioInfo   <-
    wInfo %>%
    group_by(Hj, possibleWeight) %>%
    #    summarise(wInfo = rejectedHypInfo) %>% ?? why I had this line
    mutate(hypNames = nodes[Hj])
  
  paramData <- tibble(
    hypNames = nodes,
    test.type = 1,
    k = sapply(D$infoFr, length),
    timing = D$infoFr,
    sfu    = sapply(D$grSeqTesting, function(x)
      x$sfu),
    sfupar = lapply(D$grSeqTesting, function(x)
      x$sfupar),
    sfInfo = sapply(D$grSeqTesting, print_sfInfo),
    possibleWeight
  ) %>% unnest(possibleWeight) %>%
    mutate(
      alpha = possibleWeight * sigLevel,
      weight_alpha = paste0(possibleWeight, " (", alpha, ")")
    ) %>%
    filter(alpha > 0)
  
  gsDesign_m <- function(...) {
    args <- list(...)
    if (args$k == 1) {
      # if a single analysis spend all alpha at once
      return(list(upper = list(
        bound = qnorm(args$alpha, lower.tail = FALSE)
      )))
    } else{
      return(do.call(gsDesign, args))
    }
  }
  gsDesignArgs <-
    paramData %>% dplyr::select(alpha, test.type, k, timing, sfu, sfupar)
  gsDesignList <- pmap(gsDesignArgs, gsDesign_m)
  # extract nominal p-values
  nominalPval <-
    lapply(gsDesignList, function(x)
      pnorm(x$upper$bound, lower.tail = FALSE))
  
  res <-
    paramData %>% add_column(nominalPval) %>% left_join(scenarioInfo)
  res$Analysis <-
    sapply(res$k, function(x)
      paste(1:x, collapse = "<br>"))
  res$timing <-
    sapply(res$timing, function(x)
      paste(round(x, digits = idigits), collapse = "<br>"))
  # res$nominalPvalx2 <- sapply(res$nominalPval, function(x) paste(round( 2*x, digits=5),collapse = "<br>") )
  res$nominalPval   <-
    sapply(res$nominalPval, function(x)
      paste(round(x, digits = pdigits), collapse = "<br>"))
  res$scenarioInfo <- filter(scenarioInfo, possibleWeight > 0)
  # run power evaluation if given effect size and sample size
  if (!is.null(D$theta) & !is.null(D$hypN)) {
    # extract nominal p-values
    aux <- tibble::tibble(
      hypNames = paramData$hypNames,
      k = sapply(paramData$timing, length),
      timing =    paramData$timing,
      a = lapply(gsDesignList, function(x) {
        if (!is.null(x$lower$bound))
          x$lower$bound
        else
          rep(-20, length(x$upper$bound))
      }),
      b = lapply(gsDesignList, function(x)
        x$upper$bound)
    ) %>%
      left_join(dplyr::select(D, hypNames, theta, hypN, standFactor, logDelta))
    aux$n.I <-
      mapply(FUN = "*", aux$timing, aux$hypN, SIMPLIFY = FALSE)
    
    aux$thetaHat <-
      mapply(
        FUN = function(z, n)
          z / sqrt(n),
        aux$b,
        aux$n.I,
        SIMPLIFY = FALSE
      )
    deltaHat <- mapply(function(th, c, e) {
      ans <- th / c
      if (e)
        exp(-ans)
      else
        ans
    },
    aux$thetaHat,
    aux$standFactor,
    aux$logDelta,
    SIMPLIFY = FALSE)
    
    gsProbabilityArgs <-
      aux %>%  dplyr::select(k, theta, n.I, a, b)
    gsProbabilityList <- pmap(gsProbabilityArgs, gsProbability)
    pow <-
      lapply(gsProbabilityList, function(x)
        cumsum(x$upper$prob))
    
    res <- res %>% add_column(deltaHat)
    res$deltaHat <-
      sapply(deltaHat, function(x)
        paste(round(x, digits = ddigits), collapse = "<br>"))
    
    res <- res %>% add_column(pow)
    res$pow <-
      sapply(res$pow, function(x)
        paste(round(x, digits = powdigits), collapse = "<br>"))
    
  }
  res %>% dplyr::arrange(as.numeric(gsub(".*?([0-9]+).*", "\\1", hypNames)), alpha)
}

knit_MT_table <- function(hyp_testing_dataset, digits = 5) {
  df <- hyp_testing_dataset %>%
    dplyr::select(hypNames, alpha, possibleWeight, rejectedHypInfo) %>%
    dplyr::rename(
      'Local alpha level' = alpha,
      'Weight' = possibleWeight,
      'Testing Scenario' = rejectedHypInfo,
    )
  
  if (is_html_output()) {
    df[, -1] %>%
      kable("html",
            escape = F,
            #  align=rep("c",5),
            digits = digits,
            caption = "List of possible local alpha levels following the graphical testing procedure") %>%
      kable_styling() %>%
      pack_rows(index = table(fct_inorder(df$hypNames))) %>%
      column_spec(1, latex_valign = "m")
    # collapse_rows(columns = 1, valign = "top")
  } else if (is_latex_output()) {
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    df[, -1] %>%
      mutate_all(linebreak) %>%
      kable(
        "latex",
        booktabs = T,
        escape = F,
        longtable = TRUE,
        caption = "Efficacy p-value Boundaries"
      ) %>%
      kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
      pack_rows(index = table(fct_inorder(df$hypNames)))
  } else if (knitr::pandoc_to("docx")) {
    require(flextable)
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    flextable(df)
  }
}

# prepare table that report scenarios as for cut-off nominal p-vals for all Hi at all analyses
knit_MT_grSeq_table <- function(hyp_testing_dataset, digits = 5) {
  df <- hyp_testing_dataset %>%
    dplyr::select(hypNames, alpha,
                  Analysis, timing, nominalPval,  deltaHat, pow) %>%
    dplyr::rename(
      'Local alpha level'       = alpha,
      'Info fraction'           = timing,
      'Nominal p-val (1-sided)' = nominalPval,
      # '2 x Nominal p-val'      = nominalPvalx2
      'Hurdle delta'            = deltaHat,
      'Power'                   = pow
    )
  if (is_html_output()) {
    df[,-1] %>%
      kable("html",
            escape = F,
            #  align=rep("c",5),
            digits = digits,
            caption = "Efficacy p-value Boundaries") %>%
      kable_styling() %>%
      pack_rows(index = table(fct_inorder(df$hypNames)))  %>%
      column_spec(
        column = 1,
        underline = TRUE,
        width = '3cm',
        latex_valign = "m"
      )
    # collapse_rows(columns = 1, valign = "middle")
  } else if (is_latex_output()) {
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    df[,-1] %>%
      mutate_all(linebreak) %>%
      kable(
        "latex",
        booktabs = T,
        escape = F,
        longtable = TRUE,
        caption = "Efficacy p-value Boundaries"
      ) %>%
      kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
      pack_rows(index = table(fct_inorder(df$hypNames)))
  } else if (knitr::pandoc_to("docx")) {
    require(flextable)
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    flextable(df)
  }
}

# Need to set a special spending function to account for the delayed recycling
# code Xi and Tamhane (2013) delayed alpha recycling for MT in grSeq
sfRecycle <- function (alpha, t, param)
{
  x <- list(
    name = "User-defined 2-stage",
    param = param,
    parname =  c(
      "initSigLevel",
      "tStar",
      "sfInit",
      "sfInitParam",
      # spend fun before tStar
      "sfStar",
      "sfStarParam"
    ),
    # spend fun after  tStar
    sf = sfRecycle,
    spend = NULL,
    bound = NULL,
    prob = NULL
  )
  class(x) <- "spendfn"
  
  checkScalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  checkVector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  if (length(param) != 6)
    stop("sfRecycle : parameter must be of length 6")
  checkScalar(param$initSigLevel, "numeric", c(0, Inf), c(FALSE, FALSE))
  checkVector(param$tStar, "numeric", c(0, Inf), c(TRUE, FALSE))
  
  t[t > 1] <- 1
  newSigLevel  <- alpha
  initSigLevel <- param$initSigLevel
  tStar        <- param$tStar
  sfInit       <- param$sfInit
  sfInitParam  <- param$sfInitParam
  sfStar       <- param$sfStar
  sfStarParam  <- param$sfStarParam
  
  
  if (sfStar(alpha, t = tStar, sfStarParam)$spend == sfStar(alpha, 1, sfStarParam)$spend) {
    stop("No increase of sfStar from tStar to 1\n")
    # x$spend <- ifelse(t<=tStar, sfInit(alpha,t=t,param=sfInitParam)$spend,alpha)
    # return(x)
  }
  
  # to implement (8), (9) equations in Xi (2015)
  getGammaStar <- function(y) {
    # Note that step type spending function might not be solved for gammaStar
    y - do.call(sfStar, list(
      alpha = y,
      t = tStar,
      param = sfStarParam
    ))$spend - newSigLevel + alpha1
  }
  
  if (alpha <= param$initSigLevel) {
    # no delayed recycling if no increase in available alpha
    x$spend <-
      do.call(sfInit, list(
        t     = t,
        alpha = alpha,
        param = sfInitParam
      ))$spend
    return(x)
  }
  
  alpha1 <-
    do.call(sfInit,
            list(alpha = initSigLevel, t = tStar, param = sfInitParam))$spend
  gammaStar <-
    uniroot(getGammaStar, c(.Machine$double.eps, 1e16))$root
  # gammaStar <- (newSigLevel - alpha1) / (1-log(1+(exp(1)-1)*tStar))
  x$spend <- ifelse(
    t >= tStar,
    do.call(
      sfInit,
      list(t = tStar, alpha = initSigLevel, param = sfInitParam)
    )$spend +
      do.call(sfStar, list(
        t = t,     alpha = gammaStar, param = sfStarParam
      ))$spend -
      do.call(sfStar, list(
        t = tStar, alpha = gammaStar, param = sfStarParam
      ))$spend,
    do.call(sfInit, list(
      t = t, alpha = initSigLevel, param = sfInitParam
    ))$spend
  )
  x
}

sfSuperPower <- function (alpha, t, param)
{
  checkScalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  par <-  param[1] - param[2] * alpha
  checkScalar(par, "numeric", c(0, 15), c(FALSE, TRUE))
  checkVector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  x <-
    list(
      name    = "Augmented Power",
      param   = param,
      parname = c("A", "B"),
      sf      = sfPower,
      spend   = alpha * t ^ par,
      bound   = NULL,
      prob    = NULL
    )
  class(x) <- "spendfn"
  x
}

if(FALSE){
    t <- 0:100 / 100
    sfParam <-list(
        initSigLevel = 0.05, # to spend all local alpha at IA
        tStar = 0.6,
        sfInit=sfPower,  sfInitParam = 2,
        sfStar=sfLinear, sfStarParam = c(1,1)
    )
    
    plot(t, sfPower(0.1, t, 2)$spend,
         xlab = "Proportion of sample size",
         ylab = "Cumulative Type I error spending",
         main = "Spending Function Examples",
         type = "l", cex.main = .9
    )
    
    lines(t, sfRecycle(0.1, t, param =sfParam )$spend, lty = 2)
    lines(t, sfPower( 0.05, t, param =2)$spend, lty = 3)
}  

# Example of delayed alpha recycling specs

# H1   = list(
#     sfu = sfRecycle,               # ORR
#     sfupar = list(
#         initSigLevel = w[1]*alphaTotal, # to spend all local alpha at IA
#         tStar = 0.8, 
#         sfInit=sfLinear,  sfInitParam = c(0.001,0.99),
#         sfStar=sfLinear, sfStarParam = c(1,1)
#     )
# )

if(FALSE){ # some archive thing
    unique(D$iaSpec) %>%  # to remove multi-arm hypothesis
        lapply(unlist, recursive=FALSE) %>% do.call(what=rbind) %>% data.frame() %>% tibble() %>%
        mutate(id=paste0("H",H)) %>% 
        left_join(y =D %>% dplyr::select(id,ep,hypN,iaTime),by="id") %>%
        group_by(ep) %>% arrange(id, atIF) %>%
        mutate(
            ia = row_number(atIF),
            n.I = atIF*hypN,
        )

    iaDat <- unique(unlist(D$iaSpec, recursive = FALSE)) %>% 
        lapply(unlist) %>% do.call(what=rbind) %>% data.frame() %>% tibble() %>%
        mutate(id=paste0("H",H)) %>% 
        left_join(y =D %>% dplyr::select(id,ep,hypN,iaTime),by="id") %>%
        group_by(ep) %>% arrange(id, atIF) %>%
        mutate(
            ia = row_number(atIF),
            n.I = atIF*hypN,
        )
}    

# generate gtable of timelines of IA by hypotheses
timeline_gtable <- function(D, startDate = "2022-10-12", lpi = NULL) {
  D <-
    dplyr::filter(D, regiment == unique(D$regiment)[1]) # limit plotting to just a single regiment
  # Set hypothesis names in single vector
  hypothesis_names       <- D$hypNames
  J                      <- length(hypothesis_names)
  times <- sort(unique(unlist(D$iaTime)))
  # Set number of analyses for each hypothesis
  K_j                    <- sapply(D$infoFr, length)
  # Set maximum number of analyses
  K                      <- max(K_j)
  # Set colours for each stage to use in final plot
  if (K == 2) {
    colours              <- c("#9ECAE1FF", "#3182BDFF")
  } else {
    colours              <-
      grDevices::adjustcolor(RColorBrewer::brewer.pal(n = K, name = "Blues"))
  }
  times_with_zero        <- c(0, times)
  # Set vector containing stage names
  event                  <- paste0("Stage ", 1:K)
  # Start at year 0 as will only plot time information in months from FPI (could
  # modify to use actual dates if desired)
  if (is.null(startDate))
    startDate <- "0000-01-01"
  
  start <- end <- as.POSIXct(rep(startDate, K))
  
  # Set start and end times of each stage
  for (k in 1:K) {
    start[k]             <- start[k] %m+% months(times_with_zero[k])
    end[k]               <-
      end[k] %m+% months(times_with_zero[k + 1])
  }
  # Group for the above information is the Stage and colours are those specified
  # earlier
  group                  <- rep("Stages", K)
  color                  <- colours
  # Add information on LPI
  if (!is.null(lpi)) {
    event                  <- c(event, paste0(lpi, " mo"))
    lpi_time               <-
      as.POSIXct(startDate) %m+% months(lpi)
    start                  <- c(start, lpi_time)
    end                    <- c(end, lpi_time)
    group                  <- c(group, "LPI")
    color                  <- c(color, "black")
  }
  # Loop over the hypotheses and add the information for each of them
  for (j in 1:J) {
    # Extract sample size / event information up to data maturity
    ss_j                 <-
      D$hypN[[j]] * D$infoFr[[j]]              # hypotheses[[j]]$ss[1:K_j[j]]
    # Convert to IF (rounding so it plots better)
    if_j                 <-
      round(D$infoFr[[j]] * 100) # round(100*ss_j/ss_j[K_j[j]])
    # The group is just the hypothesis name
    group                <-
      c(group, rep(D$hypNames[j], K_j[j])) # c(group, rep(hypothesis_names[j], K_j[j]))
    color                <- c(color, colours[1:K_j[j]])
    # Event information is based on combining various strings together
    if (K_j[j] > 1) {
      event_j          <- c(paste0("IA", 1:(K_j[j] - 1), ": "), "FA: ")
    } else {
      event_j          <- "FA: "
    }
    if (D$endpointType[j] == "Binomial") {
      text_j           <- "pts"
    } else if (D$endpointType[j] == "TTE") {
      text_j           <- "ev"
    }
    if (K_j[j] > 1) {
      event_j          <-
        paste0(
          event_j,
          times_with_zero[2:(K_j[j] + 1)],
          " mo,\n",
          ss_j,
          " ",
          text_j,
          c(rep(", ", K_j[j] - 1), ""),
          c(if_j[-K_j[j]], ""),
          c(rep("%IF", K_j[j] - 1), "")
        )
    } else {
      event_j          <-
        paste0(
          event_j,
          times_with_zero[2:(K_j[j] + 1)],
          " mo,\n",
          ss_j,
          " ",
          text_j,
          c(rep(", ", K_j[j] - 1), ""),
          c(if_j[-K_j[j]], ""),
          c(rep("%IF", K_j[j] - 1), "")
        )
    }
    event                <- c(event, event_j)
    start                <- c(start, end[1:K_j[j]])
    end                  <- c(end, end[1:K_j[j]])
  }
  # Build data frame for plotting
  timeline_data  <- data.frame(
    event = event,
    start = start,
    group = group,
    end   = end,
    color = color
  )
  
  # Set x limit for the plot based on FA time
  x_limit <- lubridate::interval(startDate,
                                 max(as.POSIXct(timeline_data$end)) %m+% months(1)) %/%  months(12)
  # Build initial timeline plot
  p_timeline             <-
    vistime::gg_vistime(timeline_data, linewidth = 5) +
    ggplot2::scale_x_datetime(
      breaks =
        seq(as.POSIXct(min(startDate)),
            max(as.POSIXct(
              timeline_data$end
            )) %m+%
              months(1), "years"),
      labels = 12 * seq(0, x_limit, 1)
    ) +
    ggplot2::xlab("Months") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  # Not always needed, but ggplot_build() can help avoid overlapping labels
  p_build                <- ggplot2::ggplot_build(p_timeline)
  p_build$data[[4]]$size <- 3
  p_build$data[[5]]$size <- 3
  p_timeline             <- ggplot2::ggplot_gtable(p_build)
  return(p_timeline)
}

checkInput <- function(inputD, G, enrollment = NULL) {
  N_rand <- cumsum(c(0, enrollment$rate * c(enrollment$duration)))
  numHyp <- nrow(inputD)
  Hs  <-
    inputD$iaSpec |> unlist(recursive = FALSE) |> map_dbl(pluck("H"))
  IFs <-
    inputD$iaSpec |> unlist(recursive = FALSE) |> map_dbl(pluck("atIF"))
  if (!all(Hs %in% 1:numHyp))
    stop("H index outside of range in inputD$iaSpec")
  if (!all(IFs > 0 &
           IFs <= 1))
    stop("Info fraction must be in (0,1]")
  # TODO
  # 1) check that enrollment and sample sizes for binary EPs are in agreement
  
  if (is.null(inputD$tag))
    stop("Missing tag field in inputD")
}

# main function that calculates data for report
exec_calc <- function(inputD) {
  checkInput(inputD, G)    # check the inputs TODO
  
  D <- inputD %>% dplyr::mutate(
    # effect on a natural scale, e.g, mean group difference or HR
    delta       = map_dbl(endpointParam, getDelta),
    deltaStr    = map_chr(endpointParam, getEffectSizeDetails),
    # get standardization factor to calculate effect sizes 'theta'
    standFactor = map2_dbl(endpointParam,  allocRatio, getStandardizingCoef),
    logDelta  = map_lgl(endpointParam, function(x)
      class(x) %in% "tte_exp"),
    theta =  if_else(logDelta,-log(delta) * standFactor, delta * standFactor),
    hypNames = paste(id, tag, sep = ": ")  # combine ID and tag names into full name
  )
  
  # derive calendar time of IAs usisng 'iaSpec', enrollment, and endpoint maturity time in 'endpointParam'
  D$iaTime <- deriveCalendarTime(D)
  # derive sample sizes (number of events)
  D$n.I   <- derive_Ns(D)
  # derive information fractions
  D$infoFr  <-
    lapply(
      D$n.I,
      FUN = function(x)
        unique(round(x / x[length(x)], idigits))
    )
  D$hypN    <- sapply(
    D$n.I,
    FUN = function(x)
      x[length(x)]
  )
  
  # add description columns
  descr_fields <-
    match(c("regiment", "ep", "suffix"), names(D), nomatch = 0)
  if (!any(descr_fields))
    descr_fields <- match(c("tag"), names(D), nomatch = 0)
  D$descr <- select(D, descr_fields) %>% pmap_chr(., paste)
  D <-
    add_column(D, grSeqTestingCh = sapply(D$grSeqTesting, print_sfInfo))
  
  ia_details <- D %>%
    dplyr::distinct(id, iaSpec, .keep_all = TRUE)     %>%
    dplyr::mutate(id_tag = paste0(id, " (", tag, ")"))   %>%
    dplyr::select(id_tag, iaSpec, iaTime, n.I, infoFr)        %>%
    dplyr::mutate(iaSpec = map(iaSpec, function(x) {
      lapply(x, as_tibble)
    })) %>%
    unnest(c(iaSpec, iaTime, n.I, infoFr)) %>%
    unnest(iaSpec) %>%
    dplyr::group_by(id_tag) %>%
    dplyr::mutate(ia = row_number(),
                  criterion = paste0("H", H, " at information fraction ", round(atIF, idigits))) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(iaTime) %>% 
    mutate(ia_ind = cur_group_id()) %>%
    dplyr::ungroup()
  
  # the main call to calculate details of MT in group sequential testing
  hyp_testing_dataset <-
    report_MT_grSeq(G, D, pdigits = pdigits, idigits = idigits)
  return(list(
    D                   = D ,
    ia_details          = ia_details,
    hyp_testing_dataset = hyp_testing_dataset
  ))
}


