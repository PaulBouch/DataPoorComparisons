# 23/03/2015 tsa.fit - tidy up argument calling of empty matrices
# 23/03/2015 tsa.fit - allow both residuals and prediction errors
# 24/03/2015 tsa.fit - allows external input of variance for surveys
# 30/04/2015 make survey model into two components: mean and variance
# 03/05/2015 tsa.input - allow natural mortality parameters and input survey variance
# 18/02/2016 tidy up error prediction and residual output code
# 18/02/2016 give parameters a scale and ndeps attribute
# 22/02/2016 replace tsa.fit with optim version
# 26/03/2016 add function to update parameters
# 01/04/2016 add fsd.y.cvmult
# 13/04/2016 allow e.g. bycatch to be read in using auxiliary construct
# 16/04/2016 tsa.update.param - fix bug when params include discard.rate.step
# 09/12/2016 tsa.input, tsa.input.data - make path platform independent
# 09/12/2016 tsa.setup - correct bug in ylim setup when only one survey
# 09/12/2016 tsa.setup - correct bug when alim in survey is greater or equal to age of plus group
# 12/12/2016 tsa.setup - rename subset as years for consistency with other structures
# 14/04/2017 tsa.fit - correct bug in stock summary when discards present - didn't add in older ages when provided
#            as individual ages rather than as a plus group

tsa.input <- function(path = ".", catch, discards, landings, survey, auxiliary) {
  
  if (missing(catch) & (missing(discards) | missing(landings))) 
    stop('either catch or landings and discards must be given')

  if (missing(auxiliary)) stop('auxiliary data must be given')
  

  # construct output structure 

  out <- vector(5, mode = "list")
  names(out) <- c("catch", "discards", "landings", "auxiliary", "survey")
  
  
  # wrapper function for input
  
  process.data <- function(infiles, required, optional = NULL, survey = FALSE) {
    ok <- required %in% names(infiles)
    if (!all(ok)) stop('data missing: ', paste(required[!ok], collapse = ", "))
    infiles <- infiles[intersect(c(required, optional), names(infiles))]
    sapply(infiles, tsa.input.data, simplify = FALSE, path = path, survey = survey)
  }
  
  
  # read in catch, discards and landings data
  
  id <- c("catch", "discards", "landings")
  id <- id[c(!missing(catch), !missing(discards), !missing(landings))]
  
  out[id] <- lapply(id, function(x) {
    cat('\n', x, ' data\n', sep = "")
    process.data(get(x), required = c("numbers", "weights"))
  })

  if (missing(catch)) {
    cat('\ncatch data\n\tinferred from discards and landings data\n')
    
    out$catch <- vector(2, mode = "list")
    names(out$catch) <- c("numbers", "weights")
    
    out$catch$numbers <- with(out, discards$numbers + landings$numbers)
    out$catch$weights <- with(out$discards, numbers * weights) + 
      with(out$landings, numbers * weights)
    
    out$catch <- within(
      out$catch, weights <- ifelse(numbers == 0, 0, weights / numbers))
  }
  
  
  # read in auxiliary data
  # natural mortality and maturity must be there, but natural mortality can be either the
  # name of a file or parameters of a power relationship
  
  cat('\nauxiliary data\n')

  stopifnot(c("maturity", "natural.mortality") %in% names(auxiliary))
  
  id <- c("maturity", "natural.mortality")
  id <- id[sapply(id, function(x) is.character(auxiliary[[x]]))]
  
  optionalID <- setdiff(names(auxiliary), c("maturity", "natural.mortality"))
  
  out$auxiliary <- process.data(auxiliary, required = id, optional = optionalID)
  
  if (!("stock.weights" %in% names(out$auxiliary))) {
    cat('\tstock weights inferred from catch weights\n')
    out$auxiliary$stock.weights <- out$catch$weights
  }  
  
  if (!("natural.mortality" %in% names(out$auxiliary))) {
    cat ("\tcomputing natural mortality from stock weights\n")
    wk.weights <- colMeans(out$auxiliary$stock.weights)
    wk.M <- with(auxiliary$natural.mortality, a * (wk.weights ** b))
    
    out$auxiliary <- within(out$auxiliary, {
      natural.mortality <- outer(rep(1, nrow(stock.weights)), wk.M)
      dimnames(natural.mortality) <- dimnames(stock.weights)
    })
  }
  
  
  # check dimensions all consistent
  
  wk <- c(out$catch, out$discards, out$landings, out$auxiliary)
  wk <- lapply(wk, dimnames)
  if (!all(sapply(wk[-1], identical, y = wk[[1]])))
      stop('dimensions of data sets incompatible')
  
  if (missing(survey)) return(out)
  

  # read in survey data - messy code at the moment
  
  out$survey <- sapply(names(survey), simplify = FALSE, FUN = function(x) {
    cat('\nsurvey', x, 'data\n')
    result <- process.data(survey[[x]], required = "indices", optional = "variances", 
                           survey = TRUE)
    if ("variances" %in% names(result)) 
      c(result$indices, variances = list(result$variances$indices))
    else
      result$indices
  })
  
  out
}


tsa.input.data <- function(infile, path = ".", survey = FALSE) {

  infile <- file.path(path, infile)
  cat ("\treading from '", infile, "'\n", sep = "")
  
  # read in year limits of data

  wk <- unlist(read.table(infile, header = FALSE, skip = 1, nrows = 1))
  if (length(wk) != 2 | wk[1] >= wk[2]) stop('second input line does not give minimum and maximum years')
  years <- as.character(seq(round(wk[1]), round(wk[2])))
  
  
  # read in age limits of data
  # also check for plus group - get rid of any plus group info for surveys, but can keep for response and 
  # auxiliary data are these are assumed to work with plus group data anyway

  wk <- unlist(read.table(infile, header = FALSE, skip = 2, nrows = 1, stringsAsFactors = FALSE))
  plus.group <- FALSE
  if (is.character(wk))
  {
    if (substring(wk[2], nchar(wk[2])) != "+") stop('error in age information')
    plus.group <- TRUE
    wk <- as.numeric(substring(wk, 1, nchar(wk) - c(0, 1)))
  }
  if (length(wk) != 2 | wk[1] >= wk[2]) stop('third input line does not give minimum and maximum ages')
  ages <- as.character(seq(round(wk[1]), round(wk[2])))
  if (plus.group) ages[length(ages)] <- paste(ages[length(ages)], "+", sep = "")

  if (survey)
  {
    wk <- unlist(read.table(infile, header = FALSE, skip = 3, nrows = 1))
    if (length(wk) != 2 | wk[1] >= wk[2] | wk[1] < 0 | wk[2] > 1) 
      stop('third input line does not give valid survey dates')
    season <- mean(wk)
  }

  data <- as.matrix(read.table(infile, header = FALSE, skip = if (survey) 4 else 3))
  if (nrow(data) != length(years) | ncol(data) != length(ages)) 
    stop ('dimensions of data inconsistent with specified years and ages')
  
  dimnames(data) <- list(years, ages)
#  if (survey & plus.group) 
#  {
#    cat('\tinfo: dropping plus group data\n')
#    data <- data[, -ncol(data)]
#  }
  
  if (survey) list(indices = data, season = season) else data
}


tsa.setup <- function(
  data, model, ylim, alim, 
  surveys.include = if (is.null(data$survey)) "none" else names(data$survey), 
  surveys.control = list(), 
  f.est, f.range,   
  recruitment = c("average", "ricker", "hockey", "random walk"), 
  recruitment.distribution = c("log normal", "normal"), 
  large.yearclass, gudmundssonH1, aux.ave = 3) {

  # constructs tsa object for assessment
  
  # argument checking  
  
  is.survey <- !identical(surveys.include, "none")
  
  if (is.survey) {

    ok <- surveys.include %in% names(data$survey)
    if (!all(ok)) stop("unknown surveys in 'surveys.include': ", 
                       paste(surveys.include[!ok], collapse = ", "))
  
    ok <- names(surveys.control) %in% names(data$survey)
    if (!all(ok)) stop("unknown surveys in 'surveys.control': ", 
                       paste(names(surveys.control[!ok]), collapse = ", "))

  }

  recruitment <- match.arg(recruitment)
  recruitment.distribution <- match.arg(recruitment.distribution)
  
  
  # set up output list
  
  out <- list(original.data = data)

  
  # set up year and age limits, with error checking

  wk <- switch(model$response, catch = data$catch$numbers, 
               discards = data$discards$numbers)
  wk.years <- as.integer(rownames(wk))
  wk.ages <- sub("+", "", colnames(wk), fixed = TRUE)
  plus.group <- !identical(colnames(wk), wk.ages)
  wk.ages <- as.integer(wk.ages)
  
  
  # slight inconsistency here - use survey data to set up ylim, but then use
  # ylim to constrain years of surveys - need to resolve, but unlikely to be
  # serious in practice
  

  if (missing(ylim)) {
    ylim <- range(wk.years)
    if (is.survey) 
      ylim[2] <- max(ylim[2], sapply(data$survey[surveys.include], function(i) 
        max(as.integer(rownames(i$indices)))))
  }

  if (ylim[1] < min(wk.years)) {
    warning('lower year limit is less than first catch (discards and landings) year ',
            'and will be adjusted')
    ylim[1] <- min(wk.years)
  }

  out$ylim <- ylim


  if (missing(alim)) alim <- range(wk.ages)

  if (alim[1] > 1 & !(recruitment %in% c("average", "random walk"))) 
    stop('lower age limit > 1: only average or random walk recruitment allowed')
  
  if (alim[2] == max(wk.ages) & !plus.group) {
    warning('upper age limit = maximum age and will be decreased by one to form a plus group', 
            .call = FALSE)
    alim[2] <- alim[2] - 1
  }

  out$alim <- alim


  # sort out response variables - need catch data anyway 

  if (model$response == "catch") {
    catch <- tsa.setup.data(
      data$catch$numbers, ylim, alim, type = "response", x.weights = data$catch$weights,
      subset.years = model$years, subset.ages = NULL, scaling = model$scaling, 
      aux.ave = aux.ave)
    out$catch <- catch
    out$model <- model
  }
  else {
    # ensure scaling is the same for discards and landings - need catch data
    # anyway for simplifying discard data

    subset.landings <- subset.discards <- subset.catch <- NULL

    if ("years" %in% names(model)) {
      if (is.numeric(model$years)) 
        subset.landings <- subset.discards <- subset.catch <- model$years
      else {
        subset.landings <- model$years$landings
        subset.discards <- model$years$discards
        subset.catch <- if (is.null(subset.discards)) subset.landings 
          else if (is.null(subset.landings)) subset.discards
          else intersect(subset.landings, subset.discards)
      }
    }

    catch <- tsa.setup.data(
      data$catch$numbers, ylim, alim, type = "response", x.weights = data$catch$weights, 
      subset.years = subset.catch, subset.ages = NULL, scaling = model$scaling, 
      aux.ave = aux.ave)

    out$landings <- tsa.setup.data(
      data$landings$numbers, ylim, alim, type = "response", 
      x.weights = data$landings$weights, subset.years = subset.landings, 
      subset.ages = NULL, scaling = catch$scaling, aux.ave = aux.ave) 

    out$discards <- tsa.setup.data(
      data$discards$numbers, ylim, alim, type = "response", 
      x.weights = data$discards$weights, subset.years = subset.discards, 
      subset.ages = NULL, scaling = catch$scaling, aux.ave = aux.ave) 

    if (!("model" %in% names(model))) stop('model must be specified in model argument')

    if (model$model == "ogive" & !("log.age" %in% names(model))) 
      stop('log.age must be specified in model argument')
    

    # deal with age ranges for both landings and discards
    
    wk <- seq(alim[1], alim[2])
  
    if (is.null(model$ages)) 
      model$ages <- list(landings = wk, discards = wk)
    else
      model <- within(model, {
        if (is.null(ages$discards)) ages$discards <- wk
        if (is.null(ages$landings)) ages$landings <- wk
        ages$discards <- sort(ages$discards)
        ages$landings <- sort(ages$landings)
        if (!setequal(wk, c(ages$discards, ages$landings)) 
            | min(ages$discards) > min(ages$landings)
            | max(ages$discards) > max(ages$landings))
          stop (
            'landings or discard ages misspecified in ages component of model argument')
      })
    
    out$model <- model
    
    
    if (length(model$ages$landings) < length(wk)) {
      id <- as.character(setdiff(wk, model$ages$landings))
      out$landings$data[, id] <- 0
      out$landings$weights[, id] <- 0
      out$landings$include[, id] <- 0
      out$landings$cvmult[, id] <- 0
      out$discards$data[, id] <- catch$data[, id]
      out$discards$weights[, id] <- catch$weights[, id]
    }
    
    if (length(model$ages$discards) < length(wk)) {
      id <- as.character(setdiff(wk, model$ages$discards))
      out$discards$data[, id] <- 0
      out$discards$weights[, id] <- 0
      out$discards$include[, id] <- 0
      out$discards$cvmult[, id] <- 0
      out$landings$data[, id] <- catch$data[, id]
      out$landings$weights[, id] <- catch$weights[, id]
    }

    if (model$model == "step") {
      if (is.null(model$step)) 
        stop ("step component missing from model argument (list giving initial ', 
              'discard ages and step year)")

      if (is.null(model$step$year) | is.null(model$step$initial.ages)) 
        stop ("year or initial.ages missing from model$step")
    
      wk.y <- seq(ylim[1], ylim[2])

      ida <- as.character(setdiff(wk, model$step$initial.ages))
      idy <- as.character(wk.y[wk.y < model$step$year])
      out$discards$data[idy, ida] <- 0
      out$discards$weights[idy, ida] <- 0
      out$discards$include[idy, ida] <- 0
      out$discards$cvmult[idy, ida] <- 0
      out$landings$data[idy, ida] <- catch$data[idy, ida]
      out$landings$weights[idy, ida] <- catch$weights[idy, ida]
    }
  }


  # set up survey data

  if (is.survey) {
    
    out$survey <- sapply(surveys.include, simplify = FALSE, FUN = function(i) {

      info <- data$survey[[i]]
      control <- if (i %in% names(surveys.control)) surveys.control[[i]] else list()

        
      # get non-plus group ages for survey that are within alim bounds
      
      ages <- colnames(info$indices)
      plusGroup <- grepl("+", ages, fixed = TRUE)
      
      if (sum(plusGroup) >= 2) stop ("multiple plus groups in survey", i)
      if (any(plusGroup) & !plusGroup[length(ages)]) 
        stop ("plus group not oldest age in survey", i)
      
      ages <- as.integer(ages[!plusGroup])
      
      if (max(ages) >= alim[2])
        warning("upper age limit of survey ", i, " reduced so that lower than plus group", call. = FALSE)
      ages <- intersect(seq(alim[1], alim[2] - 1), ages)
      
      
      # get survey years that are within ylim bounds, ensuring survey indices 
      # begin after first catch year
      
      years <- as.integer(rownames(info$indices))
      years <- intersect(seq(ylim[1] + 1, ylim[2]), years)
      
      
      # set up default model structure, alter any values specified in 
      # model.control, and check values are sensible
      
      model <- list(mean = "full", variance = "estimate", scaling = 0, ages = ages, 
                    years = years, logScale = TRUE, amat = NA)
      
      ok <- names(control) %in% names(model)
      if (any(!ok)) stop("unknown names in survey.model: ", 
                         paste(names(control)[!ok], collapse = ", "))
      
      model[names(control)] <- control
      
      match.arg(model$mean, c("full", "loglinear"))
      match.arg(model$variance, c("estimate", "known"))
      
      if (!all(model$ages %in% ages)) {
        model$ages <- intersect(model$ages, ages)
        warning("ages misspecified for survey ", i, ": following ages used: ", 
                paste(model$ages, collapse = ", "), call. = FALSE)
      }
      
      if (!all(model$years %in% years)) {
        model$years < intersect(model$years, years)
        warning("years misspecified for survey ", i, ": following years used: ", 
                paste(model$years, collapse = ", "), call. = FALSE)
      }
      
      if (model$mean == "loglinear") {
        if (is.na(model$amat)) 
          stop("amat must be specified for loglinear model for survey ", i)
        if (model$amat < min(model$ages) + 2 | model$amat > max(model$ages))
          stop("amat should be between ", min(model$ages) + 2, " and ", 
               max(model$ages), " for survey ", i)
      }
      
      result <- list(model = model)
      
      tsa.setup.args <- list(
        x = info$indices, ylim = ylim, alim = alim, type = "survey", 
        scaling = model$scaling, subset.years = model$years, subset.ages = model$ages, 
        season = info$season)
      if (model$variance %in% "known") tsa.setup.args$x.variances <- info$variances
      
      c(result, do.call("tsa.setup.data", tsa.setup.args))
    }
    )
  }

  
  # set up auxiliary data

  out$auxiliary <- sapply(names(data$auxiliary), simplify = FALSE, FUN = function(i) {
      switch(i, 
             "stock.weights" = switch(model$response, 
                                      catch = out$catch$weights, 
                                      discards = catch$weights),
             tsa.setup.data(data$auxiliary[[i]], ylim, alim, type = "auxiliary", 
                            aux.ave = aux.ave))
  })



  # set up other structures 
  
  amat <- tail(f.est, 1)
  if (amat < alim[1] + 2 | amat >= max(alim)) 
    stop('amat should be between ', alim[1] + 2, ' and ', max(alim) - 1)


  out$gudmundssonH1 <- 
    if (missing(gudmundssonH1)) rep(1, diff(alim)+1) else gudmundssonH1
  out$f.range <- f.range

  out$recruitment <- recruitment <- match.arg(recruitment)
  if (!(out$recruitment %in% c("average", "random walk")) & alim[1] > 1) 
      stop('recruitment model not implemented for minimum age > 1')
  
  out$recruitment.distribution <- recruitment.distribution 
  
  if (!missing(large.yearclass)) out$large.yearclass <- large.yearclass

  out <- within(out, {
    fsd.v.cvmult <- fsd.y.cvmult <- rep(1, diff(ylim)+1)
    names(fsd.v.cvmult) <- names(fsd.y.cvmult) <- as.character(seq(ylim[1], ylim[2]))
  })

  param.setup <- function(id) 
    matrix(NA, nrow = length(id), ncol = 6, 
           dimnames = list(id, c("estimate", "lower", "upper", "scale", "ndeps", "active")))
  
  out$f.est <- f.est

  out$param <- list(
    fishing.selection = param.setup(paste("age", f.est)),
    fishing.sd = param.setup(c("F", "U", "V", "Y")),
    cv = switch(model$response, 
                catch = param.setup("catch"), 
                discards = param.setup(c("landings", "discards"))), 
    recruitment = switch(recruitment, average = param.setup(c("average", "cv")), 
      ricker = param.setup(c("alpha (slope at origin)", "beta (density dep)", "cv")),
      hockey = param.setup(c("max recruits", "ssb (change point)", "cv")),
      "random walk" = param.setup(c("log mean recruit", "sd persistent", "cv")))
  )    

  # only relevant if random walk

  out$recruit.cvmult <- rep(1, diff(ylim)+1)
  names(out$recruit.cvmult) <- as.character(seq(ylim[1], ylim[2]))
  

  if (model$response == "discards") {
    
    if (model$model %in% c("full", "step")) 
      out$param$discard.rate.sd <- param.setup(c("transitory", "persistent"))

    if (model$model == "step") {
      wk.ages <- with(model$ages, intersect(discards, landings))
      out$param$discard.rate.step <- param.setup(paste("age", wk.ages))
    }
        
    if (model$model == "ogive")
      out$param$discard.ogive.sd <- param.setup(
        c("transitory residual", "transitory intercept", "persistent intercept", 
          "transitory slope", "persistent slope"))
  }  

  
  if (is.survey) {
    out$param$survey <- sapply(
      out$survey, USE.NAMES = TRUE, simplify = FALSE, function(i) {

        ages <- i$model$ages
        ages <- switch(i$model$mean, 
                       full = ages, 
                       loglinear = c(min(ages), min(ages) + 1, i$model$amat))

        result <- list(selection = param.setup(paste("age", ages)))
        
        cv <- switch(i$model$variance, 
                     estimate = c("sigma", "eta", "omega", "beta"),
                     known = c("omega", "beta"))
        
        result$cv <- param.setup(cv)

        result
      }
    )
  }
  

  if ("misrep" %in% names(model)) 
    out$param$misrep <- param.setup(c("sd.transitory", "sd.persistent"))
  
  out
}  


tsa.setup.data <- function(
  x, ylim, alim, type = c("response", "survey", "auxiliary"), 
  x.weights, subset.years, subset.ages, scaling, season, x.variances, aux.ave) {

  years <- as.character(seq(ylim[1], ylim[2]))
  ages <- as.character(seq(alim[1], alim[2]))

  data <- include <- cvmult <- weights <- 
    matrix(0, nrow = length(years), ncol = length(ages), dimnames = list(years, ages))

  row.id <- intersect(rownames(x), years)
  
  colnames(x) <- sub("+", "", colnames(x), fixed = TRUE)
  col.id <- intersect(colnames(x), ages)
  data[row.id, col.id] <- x[row.id, col.id]


  # add in weight data, and adjust for plus group for response data
  
  if (type == "response")
  {
    colnames(x.weights) <- colnames(x)

    weights[row.id, col.id] <- x.weights[row.id, col.id]
  
    max.age <- max(as.numeric(colnames(x)))
    ia <- as.character(alim[2])
    ja <- as.character(alim[2]:max.age)

    wk <- rowSums(x[row.id, ja, drop = FALSE])
    data[row.id, ia] <- wk
    weights[row.id, ia] <- 
      rowSums(x.weights[row.id, ja, drop = FALSE] * x[row.id, ja, drop = FALSE]) / wk
    weights[row.id, ia][wk == 0] <- 0
  }


  # calculate scaling and set up include and cvmult structures for response and survey data
  # ensure survey data start after first catch year

  if (type %in% c("response", "survey"))
  {
    include[row.id, col.id] <- 1

    if (!is.null(subset.years))
    {
      include[setdiff(rownames(include), as.character(subset.years)), ] <- 0
      if (type == "response" & any(include[1,] == 0)) stop('all data should be included in first year')
    }     

    if (type == "survey")
      include[, setdiff(colnames(include), as.character(subset.ages))] <- 0

    if (all(include == 0)) stop('no data included')

    cvmult <- include
    
    if (is.null(scaling)) 
      scaling <- - ceiling(median(log10(data[include == 1 & data > 0]), na.rm = TRUE))
    data <- data * (10^scaling)

    if (type %in% "survey" & !missing(x.variances)) {
      cvmult[row.id, col.id] <- sqrt(x.variances)[row.id, col.id] * (10^scaling) 
      cvmult <- ifelse(include == 1, cvmult, 0)
    }
  }
  

  # deal with auxiliary data (and weight data) in non-response years

  if (type %in% c("response", "auxiliary"))
  {
    ymax <- max(as.numeric(row.id))
    if (ymax < ylim[2])
    {
      wk <- switch(type, response = weights, auxiliary = data)
      for (iy in seq(ymax + 1, ylim[2]))
        wk[as.character(iy), ] = colMeans(wk[as.character(seq(ymax - (aux.ave - 1), ymax)), ])
      if (type == "response") weights <- wk else data <- wk
    }
  }

  switch(
    type, 
    response = list(data = data, include = include, cvmult = cvmult, weights = weights, scaling = scaling), 
    survey = list(data = data, season = season, include = include, cvmult = cvmult, scaling = scaling),
    auxiliary = data)
}


tsa.fit <- function(tsa.ob, dll = "TSA.dll", filterType = c("estimate", "forward", "backward"), 
                    hessian = FALSE, restart) {
  
  filterType <- match.arg(filterType)
  
  nyear <- diff(tsa.ob$ylim) + 1
  nage <- diff(tsa.ob$alim) + 1
  
  null.mat <- matrix(0, nrow = nyear, ncol = nage)
  
  # cols 1 and 2 are for misreporting, 3 and 4 for random walk recruitment
  
  null.stateOutput <- matrix(0, nrow = nyear, ncol = 4)  
  
  
  # sort out discards stuff - check if to be included
  
  is.discards <- tsa.ob$model$response == "discards"
  
  discard_model <- 1
  dstep <- 0
  dstepparameter <- matrix(0, nrow = nage, ncol = 6, 
                           dimnames = list(seq(tsa.ob$alim[1], tsa.ob$alim[2]), NULL))
  dparameter <- matrix(0, nrow = 5, ncol = 6)
  
  
  if (is.discards) {
    discards <- tsa.ob$discards
    
    wk.model <- tsa.ob$model
    # assume for now that log.age = FALSE
    discard_model <- switch(wk.model$model, full = 1, ogive = 2, step = 3)  
    if (wk.model$model %in% c("full", "step"))
      dparameter[1:2, ] <- tsa.ob$param$discard.rate.sd
    else
      dparameter <- tsa.ob$param$discard.ogive.sd
    
    
    if (wk.model$model == "step") {
      dstep <- wk.model$step$year
      wk.ages <- with(tsa.ob$model$ages, intersect(landings, discards))
      dstepparameter[as.character(wk.ages), ] <- tsa.ob$param$discard.rate.step
    }
    
    landings <- tsa.ob$landings
  }
  else {
    discards <- list(data = null.mat, include = null.mat, cvmult = null.mat, weights = null.mat)
    landings <- tsa.ob$catch
  }
  
  
  # sort out survey stuff - first check if survey data are to be included and
  # provide dummy values if not
  
  is.survey <- "survey" %in% names(tsa.ob)
  
  nsurvey <- if (!is.survey) 1 else length(tsa.ob$survey)
  
  null.smat <- array(0, dim = c(nsurvey, nyear, nage))
  survey <- survey_include <- survey_cvmult <- null.smat
  survey_season <- rep(0, nsurvey)
  cparameter <- array(0, dim = c(nsurvey, nage, 6))
  sparameter <- array(0, dim = c(nsurvey, 4, 6))
  survey_model <- survey_amat <- survey_logScale <- rep(0, nsurvey)
  
  if (is.survey) {
    for (i in 1:nsurvey) {
      
      wk <- tsa.ob$survey[[i]]
      survey[i,,] <- wk$data
      survey_include[i,,] <- wk$include
      survey_cvmult[i,,] <- wk$cvmult
      survey_season[i] <- wk$season
      
      survey_model[i] <- switch(wk$model$mean, full = 1, loglinear = 2) +
        switch(wk$model$variance, estimate = 0, known = 2)
      
      wk.param <- tsa.ob$param$survey[[i]]
      
      if (wk$model$mean %in% "full") {
        wk.age <- wk$model$ages - tsa.ob$alim[1] + 1
        cparameter[i, wk.age, ] <- wk.param$selection
      } 
      else if (wk$model$mean %in% "loglinear") {
        cparameter[i, 1:3, ] <- wk.param$selection
        survey_amat[i] <- wk$model$amat
      }
      
      if (wk$model$variance %in% "estimate") {
        sparameter[i, , ] <- wk.param$cv
      }
      else if (wk$model$variance %in% "known") {
        sparameter[i, 3:4, ] <- wk.param$cv
        sparameter[i,1,1] <- 1
      }
      
      survey_logScale[i] <- as.integer(wk$model$logScale)
    }
  }
  
  
  is.seals <- "seals" %in% names(tsa.ob)
  if (is.seals) {
    seals <- tsa.ob$seals
    seals$parameter <- tsa.ob$param$seals
    seals$model <- 'constant rate' 
  }
  else {
    null.mat <- matrix(0, nrow = nyear, ncol = nage)
    seals <- list(
      data = null.mat, include = null.mat, cvmult = null.mat, abundance = rep(0, nyear),
      parameter = matrix(0, nrow = 6, ncol = 6), model = 'constant rate')
  }
  
  
  large.recruitment <- rep(0, nyear)
  if ("large.yearclass" %in% names(tsa.ob)) {
    names(large.recruitment) <- seq(tsa.ob$ylim[1], tsa.ob$ylim[2])
    large.recruitment[as.character(tsa.ob$large.yearclass + tsa.ob$alim[1])] <- 1
  }  
  
  
  if ("misrep" %in% names(tsa.ob$model))
    misrep <- list(ylim = tsa.ob$model$misrep, parameter = tsa.ob$param$misrep)
  else
    misrep <- list(ylim = c(1, 0), parameter = matrix(0, nrow = 2, ncol = 6))

  
  # combines these together to be passed through to the fortran code
    
  other.cvmult <- with(tsa.ob, cbind(fsd.v.cvmult, fsd.y.cvmult, recruit.cvmult))


  fselectivity <- tsa.ob$param$fishing.selection
  fsd <- tsa.ob$param$fishing.sd
  cvparameter <- tsa.ob$param$cv
  rparameter <- tsa.ob$param$recruitment
  
  nrparameter <- nrow(rparameter)
  
  setupParam <- function(id) {
    
    out <- c(
      fselectivity[, id], 
      fsd[, id], 
      cvparameter[, id], 
      rparameter[, id]
    )

    if (is.discards) {
      out <- c(out, switch(
        tsa.ob$model$model, 
        full = dparameter[1:2, id],
        ogive = dparameter[, id],
        step = dparameter[1:2, id]
      ))
      
      if (tsa.ob$model$model %in% "step") {
        wk.ages <- with(tsa.ob$model$ages, intersect(landings, discards))
        out <- c(out, dstepparameter[as.character(wk.ages), id])
      }
    }

    if (is.survey) {
      for (iSurvey in 1:nsurvey) {
        
        model <- tsa.ob$survey[[iSurvey]]$model
        
        out <- c(out, switch(
          model$mean, 
          full = {
            ages <- model$ages - tsa.ob$alim[1] + 1
            cparameter[iSurvey, ages, id]
          },  
          loglinear = cparameter[iSurvey, 1:3, id]
        ))
        
        out <- c(out, switch(
          model$variance, 
          estimate = sparameter[iSurvey, , id],
          known = sparameter[iSurvey, 3:4, id] 
        ))
      }
    }
    
    if ("misrep" %in% names(tsa.ob$model))
      out <- c(out, misrep$parameter[, id])

    out
  }
  
  parAll <- setupParam(1)
  lower <- setupParam(2)
  upper <- setupParam(3)
  scale <- setupParam(4)
  ndeps <- setupParam(5)
  estimate <- setupParam(6)
  
  estimate <- estimate == 1
  start <- parAll[estimate]
  lower <- lower[estimate]
  upper <- upper[estimate]
  scale <- scale[estimate]
  ndeps <- ndeps[estimate]
  
  scale[is.na(scale)] <- 0.01
  ndeps[is.na(ndeps)] <- 0.001
  
  if (!missing(restart)) start <- restart
  
  likCalc <- function(par, allOutput = FALSE) {
    
#    save(par, file = "par monitoring.RData")
    
    parAll[estimate] <- par
    
    ict <- 0
    
    n <- nrow(fselectivity)
    fselectivity[, 1] <- parAll[ict + 1:n]
    ict <- ict + n
    
    n <- nrow(fsd)
    fsd[, 1] <- parAll[ict + 1:n]
    ict <- ict + n
    
    n <- nrow(cvparameter)
    cvparameter[, 1] <- parAll[ict + 1:n]
    ict <- ict + n
    
    n <- nrow(rparameter)    
    rparameter[, 1] <- parAll[ict + 1:n]
    ict <- ict + n
 
    if (is.discards) {
      if (tsa.ob$model$model %in% c("full", "step")) {
        n <- 2
        dparameter[1:2, 1] <- parAll[ict + 1:n]
        ict <- ict + n
      }     
      else if (tsa.ob$model$model %in% "ogive") {
        n <- nrow(dparameter)
        dparameter[, 1] <- parAll[ict + 1:n]
        ict <- ict + n
      }

      if (tsa.ob$model$model %in% "step") {
        ages <- with(tsa.ob$model$ages, intersect(landings, discards))
        n <- length(ages)
        dstepparameter[as.character(ages), 1] <- parAll[ict + 1:n]
        ict <- ict + n
      }
    }
    
    if (is.survey) {
      for (iSurvey in 1:nsurvey) {
        
        model <- tsa.ob$survey[[iSurvey]]$model
        
        if (model$mean == "full") {
          ages <- model$ages - tsa.ob$alim[1] + 1
          n <- length(ages)
          cparameter[iSurvey, ages, 1] <- parAll[ict + 1:n]
          ict <- ict + n
        }
        else if (model$mean == "loglinear") {
          n <- 3
          cparameter[iSurvey, 1:n, 1] <- parAll[ict + 1:n]
          ict <- ict + n
        }
        
        if (model$variance == "estimate") {
          n <- 4
          sparameter[iSurvey, , 1] <- parAll[ict + 1:n]
          ict <- ict + n
        }
        else if (model$variance == "known") {
          n <- 2
          sparameter[iSurvey, 3:4, 1] <- parAll[ict + 1:n]
          ict <- ict + n
        }
      }
    }
      
    if ("misrep" %in% names(tsa.ob$model)) {
      n <- 2
      misrep$parameter[, 1] <- parAll[ict + 1:n]
      ict <- ict + n
    }

    id <- c(1, 2, 3, 6)
    
    onlyFilter <- if (allOutput) 1 else 2
    
    out <- .Fortran(
      "tsa", 
      ylim = as.integer(tsa.ob$ylim), 
      alim = as.integer(tsa.ob$alim), 
      nsurvey = as.integer(nsurvey), 
      landings = as.double(landings$data), 
      linclude = as.integer(landings$include), 
      lcvmult = as.double(landings$cvmult), 
      lweights = as.double(landings$weights),
      discards = as.double(discards$data), 
      dinclude = as.integer(discards$include), 
      dcvmult = as.double(discards$cvmult), 
      dweights = as.double(discards$weights),
      survey = as.double(survey), 
      sinclude = as.integer(survey_include), 
      scvmult = as.double(survey_cvmult), 
      sseason = as.double(survey_season),
      naturalmortality = as.double(tsa.ob$auxiliary$natural.mortality), 
      maturity = as.double(tsa.ob$auxiliary$maturity), 
      stockweights = as.double(tsa.ob$auxiliary$stock.weights), 
      fmodel = as.integer(c(length(tsa.ob$f.est), tail(tsa.ob$f.est, 1))), 
      fselectivity = as.double(fselectivity[, id]), 
      fsd = as.double(fsd[, id]),      
      other_cvmult = as.double(other.cvmult),  
      recruitmentmodel = as.integer(c(
        switch(tsa.ob$recruitment, average = 1, ricker = 2, hockey = 5, 
               "random walk" = 6), 
        switch(tsa.ob$recruitment.distribution, "log normal" = 1, normal = 2))), 
      nrparameter = nrparameter, 
      rparameter = as.double(rparameter[, id]), 
      largerecruitment = as.integer(large.recruitment), 
      ncvparameter = as.integer(nrow(tsa.ob$param$cv)), 
      cvparameter = as.double(cvparameter[, id]),
      surveymodel = as.integer(survey_model), 
      cparameter = as.double(cparameter[, , id]), 
      sparameter = as.double(sparameter[, , id]),
      survey_amat = as.integer(survey_amat), 
      survey_logScale = as.integer(survey_logScale),
      seals = as.double(seals$data),
      sealsInclude = as.integer(seals$include),
      sealsCVmult = as.double(seals$cvmult),
      sealsAbundance = as.double(seals$abundance),
      sealsModel = as.integer(
        switch(seals$model, 'constant rate' = 1, 'constant weight' = 2)),
      nSealsParameter = as.integer(
        switch(seals$model, 'constant rate' = 6, 'constant weight' = 6)),
      sealsParameter = as.double(seals$parameter[, id]),
      gudmundssonH1 = as.double(tsa.ob$gudmundssonH1), 
      dmodel = as.integer(discard_model),
      dparameter = as.double(dparameter[, id]),
      dstep = as.integer(dstep),
      dstepparameter = as.double(dstepparameter[, id]),
      misreportingYlim = as.integer(misrep$ylim),
      misreportingParameter = as.double(misrep$param[, id]),
      mean_f_range = as.integer(tsa.ob$f.range),
      control = as.double(unlist(control)), 
      loglik = double(1), 
      lresiduals = double(2 * nyear * nage), 
      dresiduals = double(2 * nyear * nage), 
      sresiduals = double(2 * nsurvey * nyear * nage), 
      sealsResiduals = double(2 * nyear * nage),
      n_est = double(nyear * nage), 
      n_se = double(nyear * nage),
      f_est = double(nyear * nage), 
      f_se = double(nyear * nage), 
      scatchability = as.double(array(0, dim = c(nsurvey, nyear, 4))), 
      logit_p_est = double(nyear * nage),
      logit_p_se = double(nyear * nage),
      stateOutput = as.double(null.stateOutput),
      stockSummaryOutput = as.double(array(0, dim = c(nyear, 17))),
      onlyFilter = as.integer(onlyFilter)
    )
    
    if (allOutput) return(out) else return(out$loglik)
  }  
  
  
  options(warn = -1)
  if (!is.loaded("tsa", PACKAGE = "TSA")) dyn.load(dll)
  options(warn = 0)
  
  control <- list(major_iter_lim = 0, step_limit = 0.01, fun_prec = 300, 
                  big_lig_lik = 10000)      
  
  
  
  out <- list(data = tsa.ob)
  
  
  # estimate parameters if required
  
  if (filterType == "estimate") {
    out$convergence <- optim(
      start, likCalc, method = "L-BFGS-B", lower = lower, upper = upper, 
      control = list(trace = 1, REPORT = 1, parscale = scale, ndeps = ndeps)
    )
    with(out$convergence, cat("\nconvergence:", convergence, "\nmessage:", message, "\n"))
  }  
  
  
  # run through (a final time) to get all the other output structures
  
  finalPar <- if (filterType == "estimate") out$convergence$par else start
  
  if (filterType == "forward") {
    out <- likCalc(finalPar)
    dyn.unload(dll)      
    return(out)
  }
  
  wk <- likCalc(finalPar, allOutput = TRUE)
  out$whole <- wk
  
  
  # calculate hessian if required
  
  if (hessian) {
    out$hessian <- optimHess(finalPar, likCalc, control = list(parscale = scale, ndeps = ndeps))
    out <- within(out, {
      vcov <- 2 * solve(hessian)
      se <- sqrt(diag(vcov))
    })
  }
  
  dyn.unload(dll)       # needed because there is a bug if you go from survey to no survey data
  
  
  wk.errors <- array(wk$lresiduals, dim = c(2, nyear, nage))
  
  prediction.errors <- list(catch = tsa.errorMat(wk$lresiduals, wk$ylim, wk$alim, wk$linclude, "predErrors"))
  
  residuals <- list(catch = tsa.errorMat(wk$lresiduals, wk$ylim, wk$alim, wk$linclude, "residuals"))
  
  if (is.discards)
  {
    names(prediction.errors)[1] <- "landings"
    names(residuals)[1] <- "landings"
    
    prediction.errors$discards <- tsa.errorMat(wk$dresiduals, wk$ylim, wk$alim, wk$dinclude, "predErrors")
    residuals$discards <- tsa.errorMat(wk$dresiduals, wk$ylim, wk$alim, wk$dinclude, "residuals")
  }
  
  if (is.survey)
  {
    wk.errors <- array(wk$sresiduals, dim = c(2, nsurvey, nyear, nage))
    prediction.errors$survey <- lapply(1:nsurvey, function(i) 
    {
      out <- tsa.outmat(wk.errors[1,i,,], wk$ylim, wk$alim)
      out[c(survey_include[i,,]) == 0] <- NA
      out
    })
    names(prediction.errors$survey) <- names(tsa.ob$survey)
    
    residuals$survey <- lapply(1:nsurvey, function(i) 
    {
      out <- tsa.outmat(wk.errors[2,i,,], wk$ylim, wk$alim)
      out[c(survey_include[i,,]) == 0] <- NA
      out
    })
    names(residuals$survey) <- names(tsa.ob$survey)
  }
  
  if (is.seals)
  {
    wk.errors <- array(wk$sealsResiduals, dim = c(2, nyear, nage))
    prediction.errors$seals <- tsa.outmat(wk.errors[1,,], wk$ylim, wk$alim)
    prediction.errors <- within(prediction.errors, seals[c(wk$sealsInclude) == 0] <- NA)
    
    residuals$seals <- tsa.outmat(wk.errors[2,,], wk$ylim, wk$alim)
    residuals <- within(residuals, seals[c(wk$sealsInclude) == 0] <- NA)
  }
  
  parameter.summary <- tsa.summary(wk, tsa.ob)
  
  if (hessian) {
    parameter.summary$se <- NA
    parameter.summary$se[estimate] <- out$se
  }
  
  n <- list(estimate = tsa.outmat(wk$n_est, wk$ylim, wk$alim), s.error = tsa.outmat(wk$n_se, wk$ylim, wk$alim))
  f <- list(estimate = tsa.outmat(wk$f_est, wk$ylim, wk$alim), s.error = tsa.outmat(wk$f_se, wk$ylim, wk$alim))
  
  out <- c(out, list(twice.loglik = - wk$loglik, summary = parameter.summary, 
                     prediction.errors = prediction.errors, residuals = residuals, n = n, f = f))
  
  
  stock.summary <- matrix(
    wk$stockSummaryOutput, nrow = diff(wk$ylim) + 1, ncol = 17, 
    dimnames = list(as.character(seq(wk$ylim[1], wk$ylim[2])), 
                    c(paste(rep(c("catch", "landings", "discards"), each = 3), c("", ".est", ".se"), sep = ""), 
                      paste(rep(c("meanF", "SSB", "TSB", "recruit"), each = 2), c(".est", ".se"), sep = ""))))
  
  stock.summary <- data.frame(stock.summary)
  
  ymax_hc <- as.numeric(tail(rownames(tsa.ob$original.data$catch$numbers), 1))
  ymax <- tsa.ob$ylim[2]
  
  if (ymax > ymax_hc) 
    stock.summary[as.character((ymax_hc+1):ymax), c("catch", "landings", "discards")] <- NA
  
  if (!is.discards) 
    out$stock.summary <- subset(
      stock.summary, select = -c(landings, landings.est, landings.se, discards, discards.est, discards.se))
  else {
    out$stock.summary.modelled <- stock.summary
    
    ywk <- intersect(rownames(stock.summary), rownames(tsa.ob$original.data$landings$numbers))
    # awk <- as.character(wk$alim[1]:wk$alim[2])
    
    tmp <- with(tsa.ob$original.data$landings, numbers * weights) * 10^tsa.ob$landings$scaling
    colnames(tmp) <- sub("+", "", colnames(tmp), fixed = TRUE)
    stock.summary[ywk, "landings"] <- rowSums(tmp[ywk, ])
    
    tmp <- with(tsa.ob$original.data$discards, numbers * weights) * 10^tsa.ob$landings$scaling
    colnames(tmp) <- sub("+", "", colnames(tmp), fixed = TRUE)
    stock.summary[ywk, "discards"] <- rowSums(tmp[ywk, ])
    
    stock.summary <- within(stock.summary, catch <- landings + discards)
    
    out$stock.summary <- stock.summary
  }
  
  wk.output <- matrix(wk$stateOutput, nrow = diff(wk$ylim) + 1, 4, 
                      dimnames = list(as.character(seq(wk$ylim[1], wk$ylim[2])), rep(c("est", "se"), 2)))
  
  out$misrep <- wk.output[, 1:2]
  out$meanRecruitment <- wk.output[, 3:4]
  
  if (is.discards)
  {
    out$logit_p <- list(estimate = tsa.outmat(wk$logit_p_est, wk$ylim, wk$alim), 
                        s.error = tsa.outmat(wk$logit_p_se, wk$ylim, wk$alim))
    #    out$trend_p <- list(estimate = tsa.outmat(wk$trend_p_est, wk$ylim, wk$alim), 
    #                        s.error = tsa.outmat(wk$trend_p_se, wk$ylim, wk$alim))
  }
  
  if (is.survey)
  {
    survey_catchability <- array(wk$scatchability, dim = c(nsurvey, nyear, 4))
    survey_catchability <- lapply(1:nsurvey, function(i) 
    {
      data <- data.frame(survey_catchability[i,,], row.names = as.character(seq(wk$ylim[1], wk$ylim[2])))
      names(data) <- c("omega.est", "omega.se", "beta.est", "beta.se")
      data
    })
    names(survey_catchability) <- names(tsa.ob$survey)
    out$survey_catchability <- survey_catchability
  }    
  
  out
}


tsa.outmat <- function(x, ylim, alim) 
  matrix(x, nrow = diff(ylim) + 1, ncol = diff(alim) + 1, 
         dimnames = list(as.character(seq(ylim[1], ylim[2])), as.character(seq(alim[1], alim[2]))))


tsa.errorMat <- function(errors, ylim, alim, include, type = c("predErrors", "residuals")) {
  nyear <- diff(ylim) + 1
  nage <- diff(alim) + 1
  errors <- array(errors, dim = c(2, nyear, nage))
  errors <- tsa.outmat(errors[switch(type, predErrors = 1, residuals = 2), , ], ylim, alim)
  errors[c(include) == 0] <- NA
  errors[as.character(ylim[1]),] <- NA
  errors
}


tsa.outpar <- function(x, rownames) {
  out <- data.frame(matrix(x, ncol = 4))
  names(out) <- c("estimate", "lower.bound", "upper.bound", "active")
  out$active <- round(out$active) == 1
  out$on.bound <- with(out, abs(estimate - lower.bound) < 1.0e-6 | abs(estimate - upper.bound) < 1.0e-4)
  out$estimate <- round(out$estimate, 4)
  if (!missing(rownames)) rownames(out) <- rownames
  out
}

tsa.update.engine <- function(x, type) {

  out <- data.frame(matrix(x, ncol = 4))
  names(out) <- c("estimate", "lower", "upper", "active")
  
  out$type <- type

  with(out, {
    active <- round(active) == 1
    
    on.lower <- abs(estimate - lower) < 1.0e-6
    on.upper <- abs(estimate - upper) < 1.0e-4

    update.upper <- active & on.upper
    update.lower <- active & on.lower & !(type %in% c("nonnegative", "positive") & abs(lower) < 1.0e-6)

    if (!any(update.upper | update.lower)) 
      return(list(boundUpdate = FALSE, update = as.matrix(data.frame(estimate, lower, upper))))
    
    new.upper <- ifelse(type %in% c("positive", "nonnegative"), 
                        upper * 1.2, 
                        upper + 0.2 * (upper - lower))
    new.lower <- ifelse(type %in% c("positive", "nonnegative"), 
                        lower * 0.5, 
                        lower - 0.2 * (upper - lower))
    new.lower <- ifelse(type %in% c("nonnegative") & new.lower <= 0.01, 0, new.lower)
    
    upper <- ifelse(update.upper, new.upper, upper)
    lower <- ifelse(update.lower, new.lower, lower)
    
    list(boundUpdate = TRUE, update = as.matrix(data.frame(estimate, lower, upper)))
  })
}



tsa.summary <- function(tsa.out, tsa.ob) {

  out <- with(tsa.out, {
    rbind(
      tsa.outpar(fselectivity, paste("F age", tsa.ob$f.est)),
      tsa.outpar(fsd, paste("sd", c("F", "U", "V", "Y"))),
      switch(
        tsa.ob$model$response, 
        catch = tsa.outpar(cvparameter, "cv catch"),
        discards = tsa.outpar(cvparameter, c("cv landings", "cv discards"))), 
      switch(
        tsa.ob$recruitment, 
        average = tsa.outpar(rparameter, c("recruit mean", "recruit cv")),   
        ricker = tsa.outpar(rparameter, c("recruitment alpha (slope at origin)", 
                                          "recruitment beta (density dependence)", "recruitment cv")),
        hockey = tsa.outpar(rparameter, c("recruitment value at change point", 
                                          "recruitment ssb at change point", 
                                          "recruitment cv")),
        "random walk" = tsa.outpar(rparameter, c("log mean recruitment at start", "sd of random walk", 
                                                 "recruitment cv")))
    )
  })

  if (tsa.ob$model$response == "discards")
  {
    if (tsa.ob$model$model %in% c("full", "step"))
    {
      dparameter <- matrix(tsa.out$dparameter, ncol = 4)[1:2, ]
      out <- rbind(out, tsa.outpar(dparameter, paste("discards sd", rownames(tsa.ob$param$discard.rate.sd))))
    }
    else
      out <- rbind(out, tsa.outpar(tsa.out$dparameter, 
                                   paste("discards sd", rownames(tsa.ob$param$discard.ogive.sd))))

    if (tsa.ob$model$model == "step")
    {
      dstepparameter <- matrix(tsa.out$dstepparameter, ncol = 4, 
                               dimnames = list(with(tsa.ob, seq(alim[1], alim[2])), NULL))
      wk.ages <- with(tsa.ob$model$ages, intersect(discards, landings))
      out <- rbind(out, tsa.outpar(dstepparameter[as.character(wk.ages), ], 
                                   paste("step", rownames(tsa.ob$param$discard.rate.step))))
    }
  }
    
  if ("survey" %in% names(tsa.ob))
  {
    nsurvey <- length(tsa.ob$survey)
    nage <- diff(tsa.out$alim) + 1
    
    cparameter <- array(tsa.out$cparameter, dim = c(nsurvey, nage, 4))
    sparameter <- array(tsa.out$sparameter, dim = c(nsurvey, 4, 4))
   
    wk <- lapply(1:nsurvey, function(i)
    {
      param <- tsa.ob$param$survey[[i]]
      model <- tsa.ob$survey[[i]]$model
      id <- names(tsa.ob$survey)[i]

      cv <- tsa.outpar(
        switch(model$variance, estimate = sparameter[i,,], known = sparameter[i,3:4,]), 
        paste(id, rownames(param$cv)))
      
      n <- nrow(param$selection)
      wk.ages <- switch(model$mean, 
                        full = model$ages - tsa.ob$alim[1] + 1, 
                        loglinear = 1:3)
      selection <- tsa.outpar(cparameter[i, wk.ages, ], 
                              paste(id, "selection", rownames(param$selection)))

      rbind(selection, cv)
    })
    wk <- do.call("rbind", wk)
    out <- rbind(out, wk)
  }
  
  if ("seals" %in% names(tsa.ob))
  {
    wk <- c("cv", paste("age", 1:4), "CodDep")
    out <- rbind(out, tsa.outpar(tsa.out$sealsParameter, wk))
  }

  if ("misrep" %in% names(tsa.ob$model))
    out <- rbind(out, tsa.outpar(tsa.out$misreportingParameter, c("misrep transitory", "misrep persistent")))
  
  out
}


tsa.update.estimates <- function(tsa.ob, tsa.fit) {
  
  getEstimates <- function(x) matrix(x, ncol = 4)[, 1]
  
  param <- within(tsa.ob$param, {
    fishing.selection[, 1] <- getEstimates(tsa.fit$fselectivity)
    fishing.sd[, 1] <- getEstimates(tsa.fit$fsd)
    cv[, 1] <- getEstimates(tsa.fit$cvparameter)
    recruitment[, 1] <- getEstimates(tsa.fit$rparameter)
  })

  if (tsa.ob$model$response == "discards") {
    if (tsa.ob$model$model %in% c("full", "step")) 
      param$discard.rate.sd[, 1] <- getEstimates(tsa.fit$dparameter)[1:2]
    else
      param$discard.ogive.sd[, 1] <- getEstimates(tsa.fit$dparameter)

    if (tsa.ob$model$model == "step")
      param$discard.rate.step[, 1] <- getEstimates(tsa.fit$dstepparameter)
  }
  
  if ("survey" %in% names(tsa.ob)) {

    nsurvey <- length(tsa.ob$survey)
    nage <- diff(tsa.ob$alim) + 1
    
    cparameter <- array(tsa.fit$cparameter, dim = c(nsurvey, nage, 4))
    sparameter <- array(tsa.fit$sparameter, dim = c(nsurvey, 4, 4))
    
    param$survey <- lapply(1:nsurvey, function(i) {
      out <- tsa.ob$param$survey[[i]]
      model <- tsa.ob$survey[[i]]$model

      out$cv[, 1] <- switch(model$variance, estimate = sparameter[i, , 1], known = sparameter[i, 3:4, 1])

      ages <- switch(model$mean, full = model$ages - tsa.ob$alim[1] + 1, loglinear = 1:3)
      out$selection[, 1] <- cparameter[i, ages, 1]
      
      out
    })
    
    names(param$survey) <- names(tsa.ob$param$survey)
  }
  
  if ("seals" %in% names(tsa.ob)) 
    param$seals[, 1] <- getEstimates(tsa.fit$sealsParameter)

  if ("misrep" %in% names(tsa.ob$model))
    param$misrep[, 1] <- getEstimates(tsa.fit$misreportingParameter)
  
  param
}


tsa.update.param <- function(tsa.ob, tsa.fit, update = c("both", "estimates", "bounds")) {

  update <- match.arg(update)
  
  id <- switch(update, both = 1:3, estimates = 1, bounds = 2:3)
  
  boundUpdates <- FALSE
  
  param <- tsa.ob$param
  
  # usually need only the fortran output, but if surveys have disappeared from the setup, then also need 
  # survey names
  
  if ("survey" %in% names(tsa.ob$param)) 
    tsa.fit <- c(tsa.fit$whole, list(surveyNames = names(tsa.fit$data$survey)))
  else
    tsa.fit <- tsa.fit$whole
  
  wk <- tsa.update.engine(tsa.fit$fselectivity, type = "positive")
  boundUpdates <- any(boundUpdates, wk$boundUpdate)
  param$fishing.selection[, id] <- wk$update[, id]

  wk <- tsa.update.engine(tsa.fit$fsd, type = "nonnegative")
  boundUpdates <- any(boundUpdates, wk$boundUpdate)
  param$fishing.sd[, id] <- wk$update[, id]

  wk <- tsa.update.engine(tsa.fit$cvparameter, type = "nonnegative")
  boundUpdates <- any(boundUpdates, wk$boundUpdate)
  param$cv[, id] <- wk$update[, id]

  wk <- tsa.update.engine(
    tsa.fit$rparameter, 
    type = switch(tsa.ob$recruitment, 
                  average = c("positive", "nonnegative"),
                  ricker = c("positive", "nonnegative", "nonnegative"),
                  hockey = c("positive", "positive", "nonnegative"), 
                  "random walk" = c("any", "nonnegative", "nonnegative")
    )
  )
  boundUpdates <- any(boundUpdates, wk$boundUpdate)
  param$recruitment[, id] <- wk$update[, id]

  if (tsa.ob$model$response == "discards") {

    if (tsa.ob$model$model %in% c("full", "step")) {
      wk <- tsa.update.engine(matrix(tsa.fit$dparameter, ncol = 4)[1:2, ], type = "nonnegative")
      boundUpdates <- any(boundUpdates, wk$boundUpdate)
      param$discard.rate.sd[, id] <- wk$update[, id]
    }
    else {
      wk <- tsa.update.engine(tsa.fit$dparameter, type = "nonnegative")
      boundUpdates <- any(boundUpdates, wk$boundUpdate)
      param$discard.ogive.sd[, id] <- wk$update[, id]
    }

    if (tsa.ob$model$model == "step") {
      ages <- with(tsa.ob$model$ages, intersect(landings, discards))
      ages <- ages - tsa.ob$alim[1] + 1
      wk <- tsa.update.engine(matrix(tsa.fit$dstepparameter, ncol = 4)[ages, ], type = "any")
      boundUpdates <- any(boundUpdates, wk$boundUpdate)
      param$discard.rate.step[, id] <- wk$update[, id]
    }

  }
      
  if ("survey" %in% names(tsa.ob)) {

    nsurvey <- length(tsa.ob$survey)
    nage <- diff(tsa.ob$alim) + 1
    
    nOldSurvey <- length(tsa.fit$surveyNames)
    cparameter <- array(tsa.fit$cparameter, dim = c(nOldSurvey, nage, 4))
    sparameter <- array(tsa.fit$sparameter, dim = c(nOldSurvey, 4, 4))
  
    surveyOK <- tsa.fit$surveyNames %in% names(tsa.ob$survey)
    cparameter <- cparameter[surveyOK, , , drop = FALSE]
    sparameter <- sparameter[surveyOK, , , drop = FALSE]
    
    for (i in 1:nsurvey) {
      
      model <- tsa.ob$survey[[i]]$model
      
      ip <- switch(model$variance, estimate = 1:4, known = 3:4)
      wk <- tsa.update.engine(sparameter[i, ip, ], type = "nonnegative")
      boundUpdates <- any(boundUpdates, wk$boundUpdate)
      param$survey[[i]]$cv[, id] <- wk$update[, id]

      ages <- switch(model$mean, full = model$ages - tsa.ob$alim[1] + 1, loglinear = 1:3)
      wk <- tsa.update.engine(cparameter[i, ages, ], type = "positive")
      boundUpdates <- any(boundUpdates, wk$boundUpdate)
      param$survey[[i]]$selection[, id] <- wk$update[, id]
    }

  }
  
  if ("seals" %in% names(tsa.ob)) {
    wk <- tsa.update.engine(tsa.fit$sealsParameter, type = "positive")
    boundUpdates <- any(boundUpdates, wk$boundUpdate)
    param$seals[, id] <- wk$update[, id]
  }

  if ("misrep" %in% names(tsa.ob$model)) {
    wk <- tsa.update.engine(tsa.fit$misreportingParameter, type = "nonnegative")
    boundUpdates <- any(boundUpdates, wk$boundUpdate)
    param$misrep[, id] <- wk$update[, id]
  }

  list(boundUpdates = boundUpdates, param = param)
}



tsa.plot.errors <- function(data, ylab = c("prediction errors", "residuals"), var.chk = FALSE, 
                            outfile, wk.title = NULL) {

  require(lattice)

  ylab = match.arg(ylab)
  
  data <- data.frame(
    prediction.errors = c(data), 
    year = rep(as.numeric(rownames(data)), ncol(data)), 
    age = rep(as.numeric(colnames(data)), each = nrow(data)))
  data <- within(data, age <- ordered(age))
  if (!var.chk)
    wk.plot <- xyplot(prediction.errors ~ year | age, data = data, scales = list(alternating = FALSE),
                      main = wk.title, ylab = ylab,
                      panel = function(x, y) {lpoints(x, y, pch = 16); panel.abline(h = 0)})
  else
    wk.plot <- xyplot(
      sqrt(abs(prediction.errors)) ~ year | age, data = data, scales = list(alternating = FALSE),
      main = wk.title, ylab = ylab,
      panel = function(x, y) {lpoints(x, y, pch = 16); panel.loess(x, y)})

  if (!missing(outfile)) {
    win.metafile(outfile, width = 11, height = 8)
    print(wk.plot)
    invisible()
  }
  else 
    print(wk.plot)
}


tsa.update.cvmult <- function(cvmult, points, id) {
  for (wk in points) {
    iy <- as.character(wk[1])
    ia <- as.character(wk[2])
    if (!(iy %in% rownames(cvmult))) {
      warning('year ', iy, ' out of range in cvmult ', id)
      next()
    }
    if (!(ia %in% colnames(cvmult))) {
      warning('age ', ia, ' out of range in cvmult ', id)
      next()
    }
    cvmult[iy, ia] <- cvmult[iy, ia] * wk[3]
  }
  cvmult
}


tsa.update <- function(tsa.ob, large.yearclass, gudmundssonH1, cvmult, survey.cvmult, seals.cvmult, 
                       fsd.v.cvmult, fsd.y.cvmult, recruit.cvmult) {

  if (!missing(large.yearclass)) tsa.ob$large.yearclass <- large.yearclass
  if (!missing(gudmundssonH1)) tsa.ob$gudmundssonH1 <- gudmundssonH1
  
  if (!missing(cvmult) && tsa.ob$model$response == "catch")
  {
    if ("ages" %in% names(cvmult))
    {
      if (length(cvmult$ages) != ncol(tsa.ob$catch$cvmult)) stop('ages has wrong length in cvmult')
      for (i in 1:nrow(tsa.ob$catch$cvmult)) 
        tsa.ob$catch$cvmult[i,] <- tsa.ob$catch$cvmult[i,] * cvmult$ages
    }  
    if ("points" %in% names(cvmult)) 
      tsa.ob$catch$cvmult <- tsa.update.cvmult(tsa.ob$catch$cvmult, cvmult$points, "catch")
  }   

  if (!missing(cvmult) && tsa.ob$model$response == "discards")
  {
    for (id in intersect(c("discards", "landings"), names(cvmult)))
    {
      wk.cvmult <- cvmult[[id]]
      if ("ages" %in% names(wk.cvmult))
      {
        wk.ages <- tsa.ob$model$ages[[id]]
        if (length(wk.cvmult$ages) != length(wk.ages)) stop('ages has wrong length in cvmult[[', id, ']]')
        ia <- as.character(wk.ages)
        for (i in 1:nrow(tsa.ob[[id]]$cvmult)) 
          tsa.ob[[id]]$cvmult[i, ia] <- tsa.ob[[id]]$cvmult[i, ia] * wk.cvmult$ages
      }
      if ("points" %in% names(wk.cvmult))
        tsa.ob[[id]]$cvmult <- tsa.update.cvmult(tsa.ob[[id]]$cvmult, cvmult[[id]]$points, id)
    }  
  }   

  if (!missing(survey.cvmult))
  {
    for (id in intersect(names(tsa.ob$survey), names(survey.cvmult)))
    {
      wk.cvmult <- survey.cvmult[[id]]
      if ("ages" %in% names(wk.cvmult))
      {
        wk.ages <- tsa.ob$survey[[id]]$model$ages
        if (length(wk.cvmult$ages) != length(wk.ages)) 
          stop('ages has wrong length in survey.cvmult[[', id, ']]')
        ia <- as.character(wk.ages)
        for (i in 1:nrow(tsa.ob$survey[[id]]$cvmult)) 
          tsa.ob$survey[[id]]$cvmult[i, ia] <- tsa.ob$survey[[id]]$cvmult[i, ia] * wk.cvmult$ages
      }
      if ("points" %in% names(wk.cvmult))
        tsa.ob$survey[[id]]$cvmult <- tsa.update.cvmult(tsa.ob$survey[[id]]$cvmult, wk.cvmult$points, id)
    }  
  }   

  if (!missing(seals.cvmult))
  {
    if ("ages" %in% names(seals.cvmult))
    {
      wk.ages <- tsa.ob$seals$model$ages
      if (length(seals.cvmult$ages) != length(wk.ages)) stop('ages has wrong length in seals.cvmult')
      ia <- as.character(wk.ages)
      for (i in 1:nrow(tsa.ob$seals$cvmult)) 
        tsa.ob$seals$cvmult[i, ia] <- tsa.ob$seals$cvmult[i, ia] * seals.cvmult$ages
    }
    if ("points" %in% names(seals.cvmult))
      tsa.ob$seals$cvmult <- tsa.update.cvmult(tsa.ob$seals$cvmult, seals.cvmult$points, "seals")
  }   


  if (!missing(fsd.v.cvmult)) {
    if ("years" %in% names(fsd.v.cvmult)) {
      for (wk in fsd.v.cvmult$years) {
        iy <- as.character(wk[1])
        tsa.ob$fsd.v.cvmult[iy] <- tsa.ob$fsd.v.cvmult[iy] * wk[2]
      }
    }  
  }   
  
  if (!missing(fsd.y.cvmult)) {
    if ("years" %in% names(fsd.y.cvmult)) {
      for (wk in fsd.y.cvmult$years) {
        iy <- as.character(wk[1])
        tsa.ob$fsd.y.cvmult[iy] <- tsa.ob$fsd.y.cvmult[iy] * wk[2]
      }
    }  
  }   
  
  if (!missing(recruit.cvmult)) {
    if ("years" %in% names(recruit.cvmult)) {
      for (wk in recruit.cvmult$years) {
        iy <- as.character(wk[1])
        tsa.ob$recruit.cvmult[iy] <- tsa.ob$recruit.cvmult[iy] * wk[2]
      }
    }  
  }   
  
  tsa.ob
}


tsa.plot.stock.summary <- function(tsa.fit.ob, year.at, recruitScale = c("log", "raw"), outfile) {

  require(lattice)
  require(grid)

  recruitScale <- match.arg(recruitScale)
  
  wk.summary <- tsa.fit.ob$stock.summary
  
  tsa.stack <- function(data, type, label = type) {

    out <- if (type %in% c("catch", "landings", "discards")) data[paste(type, c(".est", ".se", ""), sep = "")]
      else data.frame(data[paste(type, c(".est", ".se"), sep = "")], rep(NA, nrow(data)))
    
    if (type %in% c("catch", "landings", "discards", "recruit")) out <- out   
    
    out <- data.frame(as.numeric(rownames(data)), out, type = rep(label, nrow(data)))
    names(out) = c("year", "response", "se", "extra", "type")
    out
  }

  wk.Flabel <- with(tsa.fit.ob$data, paste("mean F (ages ", f.range[1], "-", f.range[2], ")", sep = ""))
  data <- rbind(tsa.stack(wk.summary, "catch"), tsa.stack(wk.summary, "meanF", wk.Flabel), 
                tsa.stack(wk.summary, "SSB"), tsa.stack(wk.summary, "recruit", "recruitment"))
    
  is.discards <- tsa.fit.ob$data$model$response == "discards"
  if (is.discards) {
    data <- rbind(data, tsa.stack(wk.summary, "landings"), tsa.stack(wk.summary, "discards"))
    # needed to get things in the correct order
    data <- within(data, type <- factor(as.character(type)))          
  }

  idx <- data$type != wk.Flabel
  idy <- c("response", "se", "extra") 
  scaling <- with(tsa.fit.ob$data, if (is.discards) discards$scaling else catch$scaling)
  data[idx, idy] <- data[idx, idy] * 10^(-scaling-3)

  if (!is.discards)
    data$include <- c(rowSums(tsa.fit.ob$data$catch$include) > 0, rep(NA, 3 * nrow(data) / 4))
  else {
    wk.d <- rowSums(tsa.fit.ob$data$discards$include) > 0
    wk.l <- rowSums(tsa.fit.ob$data$landings$include) > 0
    data$include <- c(wk.d | wk.l, wk.d, wk.l, rep(NA, nrow(data) / 2))
  }
  

  data <- within(data, {
    lower <- response - 2 * se
    upper <- response + 2 * se
    extra[extra == 0] <- NA
  })

  if (recruitScale == "log") {
    id <- data$type == "recruitment"
    data[id, ] <- within(data[id, ], {
      lower <- response * exp(- 2 * se / response)
      upper <- response * exp(2 * se / response)
    })
  }
  
  ylim <- lapply(levels(data$type), function(i) {

    data <- subset(data, type == i)
    out <- with(data, range(lower, upper, se, extra, 0, na.rm = TRUE))
    out[1] <- max(0, out[1])
    out
  })

  if (missing(year.at))
    year.at <- list(major = pretty(data$year, n = 4), minor = pretty(data$year, n = 4))
  else if (!is.list(year.at)) 
    year.at <- list(major = year.at, minor = year.at)
  else if (!all(c("major", "minor") %in% names(year.at))) 
    stop('year.at incorrectly specified')

  wk.plot <- with(data, 
    xyplot(response ~ year | type, 
      ylim = ylim, as.table = TRUE,
      xlim = extendrange(year),
      scales = list(alternating = FALSE, y = list(relation = "free")),
      par.settings = list(axis.line = list(col = "transparent")),
      ylab = "", xlab = "",
      strip = FALSE, 
      between = list(x = 2, y = 4),
      axis = function(side, ...)
      {
        switch(side,
          bottom = 
          {
            grid.lines(x = c(0, 1), y = c(0, 0), default.units = "npc", gp = gpar(col = "black", lwd = 1)) 
            panel.axis(side = side, outside = TRUE, at = year.at$minor, 
                       labels = FALSE, line.col = "black", text.cex = 1, rot = 0, tck = 0.5)
            panel.axis(side = side, outside = TRUE, at = year.at$major, 
                       labels = current.row() == nrow(trellis.currentLayout()), line.col = "black", 
                       text.cex = 1, tck = 1, rot = 0)
          },
        left = 
        {
          grid.lines(x = c(0, 0), y = c(0, 1), default.units = "npc", gp = gpar(col = "black", lwd = 1))
          panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
        })
      },
      panel = function(x, y, subscripts)
      {
        id <- as.character(type[subscripts])[1]
        if (id %in% "recruitment")
        {
          lsegments(x, lower[subscripts], x, upper[subscripts], col = grey(0.6), lwd = 2)
          lpoints(x, y, lwd = 2, col = "red", pch = 16)
        }
        else
        {
          lpolygon(c(x, rev(x)), c(lower[subscripts], rev(upper[subscripts])), col = grey(0.8), border = FALSE)
          llines(x, y, lwd = 2, col = "red")
          lpoints(x, extra[subscripts], pch = ifelse(include, 16, 1), col = "black")
        }
      
        pushViewport(viewport(clip = "off"))
        grid.text(id, 0, 1.05, just = c("left", "bottom"), gp = gpar(cex = 1.0))
        upViewport()
      }))

    if (!missing(outfile))
    {
      win.metafile(outfile, width = 11, height = 8)
      print(wk.plot)
      #dev.off()
      invisible()
    }
    else 
      print(wk.plot)
}
  

tsa.plot.stock.recruit <- function(
  tsa.fit.ob, outfile, type = c("standard", "timeseries"), recruitScale = c("log", "raw")) {

  require(lattice)
  require(grid)

  type <- match.arg(type)
  recruitScale = match.arg(recruitScale)
  
  data <- tsa.fit.ob$stock.summary
  arec <- tsa.fit.ob$data$alim[1]
  
  is.discards <- tsa.fit.ob$data$model$response == "discards"
  scaling <- with(tsa.fit.ob$data, if (is.discards) discards$scaling else catch$scaling)

  data <- with(data, data.frame(
    SSB = c(rep(NA, arec), SSB.est * 10^(-scaling-3)), 
    recruitment = c(recruit.est * 10^(-scaling-3), rep(NA, arec)), 
    se = c(recruit.se * 10^(-scaling-3), rep(NA, arec)), 
    yearClass = c(rep(NA, arec), as.numeric(rownames(data)))))

  if (arec >= 1) data$yearClass[1:arec] <- data$yearClass[arec+1]-(arec:1)
  
  wk.plot <- with(data, switch(
    type, 
    standard = {
      wk.xlim <- extendrange(c(0, SSB))
      wk.ylim <- extendrange(c(0, recruitment))
  
      xyplot(
        recruitment ~ SSB, xlim = wk.xlim, ylim = wk.ylim,
        par.settings = list(axis.line = list(col = "transparent")),
        axis = function(side, ...) {
          switch(
            side,
            bottom = {
              grid.lines(x = c(0, 1), y = c(0, 0), default.units = "npc", gp = gpar(col = "black", lwd = 1)) 
              panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
            },
            left = {
              grid.lines(x = c(0, 0), y = c(0, 1), default.units = "npc", gp = gpar(col = "black", lwd = 1))
              panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
            })
        },
        panel = function(x, y) {

          ltext(x, y, substring(yearClass, 3), col = "black")
          
          SSB <- seq(0, wk.xlim[2], length = 100)
          wk <- SSB * 10^(scaling + 3)
          
          fit <- switch(
            tsa.fit.ob$data$recruitment, 
            average = rep(tsa.fit.ob$summary["recruit mean", "estimate"], length(wk)), 
            ricker = tsa.fit.ob$summary["recruitment alpha (slope at origin)", "estimate"] * wk * 
              exp(- tsa.fit.ob$summary["recruitment beta (density dependence)", "estimate"] * wk),
            hockey = tsa.fit.ob$summary["recruitment value at change point", "estimate"] * 
              pmin(1, wk / tsa.fit.ob$summary["recruitment ssb at change point", "estimate"]))
          fit <- fit * 10^(-scaling-3)
          
          llines(SSB, fit, col = "black", lwd = 2)
        }
      )
    }, 
    timeseries = {
      wk.xlim <- extendrange(yearClass)
      wk.ylim <- extendrange(
        c(0, switch(recruitScale, log = recruitment * exp(2 * se / recruitment), raw = recruitment + 2 * se)))
      xyplot(
        recruitment ~ yearClass, xlim = wk.xlim, ylim = wk.ylim,
        par.settings = list(axis.line = list(col = "transparent")),
        axis = function(side, ...) {
          switch(
            side,
            bottom = {
              grid.lines(x = c(0, 1), y = c(0, 0), default.units = "npc", gp = gpar(col = "black", lwd = 1)) 
              panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
            },
            left = {
              grid.lines(x = c(0, 0), y = c(0, 1), default.units = "npc", gp = gpar(col = "black", lwd = 1))
              panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
            })
        },
        panel = function(x, y) {
          if (tsa.fit.ob$data$recruitment == "random walk") {
            fit <- tsa.fit.ob$meanRecruit[, "est"]
            fit.se <- tsa.fit.ob$meanRecruit[, "se"]
            lower <- exp(fit - 2 * fit.se) * 10^(-scaling-3)
            upper <- exp(fit + 2 * fit.se) * 10^(-scaling-3)
            fit <- exp(fit) * 10^(-scaling-3)
            id <- !is.na(y)
            lpolygon(c(x[id], rev(x[id])), c(lower, rev(upper)), border = FALSE, col = grey(0.6))
            llines(x[id], fit, lwd = 2, col = "black")
          }
        
          if (recruitScale == "log")
            lsegments(x, y * exp(- 2 * se / y), x, y * exp(2 * se / y), col = grey(0.3), lwd = 2)
          else 
            lsegments(x, y - 2 * se, x, y + 2 * se, col = grey(0.3), lwd = 2)
          lpoints(x, y, lwd = 2, col = "red", pch = 16)
          
          wk <- SSB * 10^(scaling + 3)
          fit <- switch(
            tsa.fit.ob$data$recruitment, 
            average = rep(tsa.fit.ob$summary["recruit mean", "estimate"], length(wk)), 
            ricker = tsa.fit.ob$summary["recruitment alpha (slope at origin)", "estimate"] * wk * 
              exp(- tsa.fit.ob$summary["recruitment beta (density dependence)", "estimate"] * wk),
            hockey = tsa.fit.ob$summary["recruitment value at change point", "estimate"] * 
              pmin(1, wk / tsa.fit.ob$summary["recruitment ssb at change point", "estimate"]))
          fit <- fit * 10^(-scaling-3)
                  
          llines(yearClass, fit, col = "black", lwd = 2)
        }
      )
    }))     

  if (!missing(outfile)) {
    win.metafile(outfile, width = 11, height = 8)
    print(wk.plot)
    #dev.off()
    invisible()
  }
  else 
    print(wk.plot)
}
  

tsa.plot.discards <- function(tsa.fit.ob, outfile) {

  require(lattice)

  if (tsa.fit.ob$data$model$response != "discards") stop('not a discards model')
  model <- tsa.fit.ob$data$model

  data <- with(tsa.fit.ob$logit_p, 
               data.frame(logit.prop = c(estimate), se = c(s.error), 
                          year = rep(as.numeric(rownames(estimate)), ncol(estimate)), 
                          age = rep(as.numeric(colnames(estimate)), each = nrow(estimate))))
  
  data <- within(data, 
  {
    discards = c(tsa.fit.ob$data$discards$data)
    landings = c(tsa.fit.ob$data$landings$data)
    include = c(tsa.fit.ob$data$discards$include)
    obs.prop <- 100 * discards / (discards + landings)
    age <- ordered(age)
    prop <- 100 * exp(logit.prop) / (1 + exp(logit.prop))
    upper <- 100 * exp(logit.prop + 2 * se) / (1 + exp(logit.prop + 2 * se))
    lower <- 100 * exp(logit.prop - 2 * se) / (1 + exp(logit.prop - 2 * se))
  })

  
  data <- subset(data, age %in% intersect(model$ages$landings, model$ages$discards))
  if (model$model == "step") 
    data <- subset(data, !(age %in% setdiff(model$ages$discards, model$step$initial.ages) & year < model$step$year))
 
  wk.plot <- with(data, xyplot(switch(model$model, ogive = prop ~ age | as.factor(year), prop ~ year | age), 
    ylim = extendrange(c(0, 100), f = 0.04), 
    scales = list(alternating = FALSE, 
                  y = list(at = switch(model$model, ogive = c(0, 50, 100), c(0, 25, 50, 75, 100)))), 
    ylab = "proportion discarded", xlab = "",
    panel = function(x, y, subscripts) 
    {
      lpolygon(c(x, rev(x)), c(lower[subscripts], rev(upper[subscripts])), col = grey(0.8), border = FALSE)
      llines(x, y, lwd = 2, col = "red")
      lpoints(x, obs.prop[subscripts], pch = ifelse(include[subscripts] == 1, 16, 1), col = "black")
    }))
    
  if (!missing(outfile))
  {
    win.metafile(outfile, width = 11, height = 8)
    print(wk.plot)
    #dev.off()
    invisible()
  }
  else 
    print(wk.plot)
}


tsa.plot.catchability <- function(tsa.fit.ob, outfile, ...) {

  require(lattice)

  if (!("survey" %in% names(tsa.fit.ob$data))) stop('survey data not included in model fit')

  data <- tsa.fit.ob$survey_catchability
  survey.id <- names(data)
  
  data <- do.call("rbind", 
                  lapply(names(data), function(i) 
                    data.frame(data[[i]], year = as.numeric(rownames(data[[i]])), survey = i)))
  data <- within(data, 
  {
    omega.se[is.na(omega.se)] <- 0   # deals with some numerically v small numbers which come out as NaNs
    beta.se[is.na(beta.se)] <- 0
    lower <- beta.est - 2 * beta.se
    upper <- beta.est + 2 * beta.se
  })

  for (i in survey.id)
  {
    include <- tsa.fit.ob$data$survey[[i]]$include
    year.min <- min(as.numeric(rownames(include)[rowSums(include) > 0]))
    data <- subset(data, !(survey == i & year < year.min))
  }
 
  wk.ylim <- extendrange(with(data, range(lower, upper, omega.est)), f = 0.04)
  wk.lab <- pretty(100 * (exp(wk.ylim) - 1))
  wk.lab <- wk.lab[wk.lab > -100]
  wk.at <- log(wk.lab / 100 + 1)
  
  wk.plot <- with(data, xyplot(beta.est ~ year | survey, ylim = wk.ylim,

    ##### Added to change the survey names (AJ)
    strip = strip.custom(factor.levels = c("ScoGFS-WIBTS-Q1", "ScoGFS-WIBTS-Q4", "IGFS-WIBTS-Q4", "UK-SCOWCGFS-Q1", "UK-SCOWCGFS-Q4")), 
 
                               ylab = "percentage change in catchability", xlab = "", ..., 
    scales = list(alternating = FALSE, y = list(at = wk.at, labels = wk.lab)),
    panel = function(x, y, subscripts) 
    {
      lpolygon(c(x, rev(x)), c(lower[subscripts], rev(upper[subscripts])), col = grey(0.8), border = FALSE)
      panel.abline(h = 0, lty = "dashed")
      llines(x, y, lwd = 2, col = "red")
      
      lpoints(x, omega.est[subscripts], pch = 16, col = "black")
    }))
    
  if (!missing(outfile))
  {
    win.metafile(outfile, width = 11, height = 8)
    print(wk.plot)
    #dev.off()
    invisible()
  }
  else 
    print(wk.plot)
}


tsa.setup.seals <- function(tsa.ob, infile, model) {

  seals <- read.csv(infile, row.names = "year")
  colnames(seals) <- gsub("age", "", colnames(seals))
  
  if (as.numeric(rownames(seals)[1]) > tsa.ob$ylim[1]) stop('no seal abundance data at start of time series')

  max.y <- as.numeric(tail(rownames(seals), 1))
  if (max.y < tsa.ob$ylim[2]) 
  {
    'info: estimating seal abundance in later years by latest estimate'
    seals[as.character((max.y+1):tsa.ob$ylim[2]),] <- NA
    seals[as.character((max.y+1):tsa.ob$ylim[2]), "abundance"] <- seals[as.character(max.y), "abundance"]
  }

  seals <- seals[as.character(tsa.ob$ylim[1]:tsa.ob$ylim[2]), c("abundance", tsa.ob$alim[1]:tsa.ob$alim[2])]


  seals <- list(abundance = seals$abundance, data = subset(seals, select = - abundance))
  
  tsa.ob$seals <- within(seals, 
  {
    names(abundance) <- rownames(data)

    include <- !is.na(data)
    include[, setdiff(colnames(include), as.character(model$ages))] <- 0

    cvmult <- include


    # data are consumption per seal, seal numbers in thousands, so to convert to same scale as landings and 
    # discards data (which are thousands * 10^scaling; here millions), need to multiply individual seal 
    # consumption by number of seals and by 10^scaling

    data <- apply(data, 2, function(i) 
    {
      out <- i
      out[is.na(i)] <- 0
      out * abundance * (10^tsa.ob$landings$scaling)
    })

    # needed to get catchabilities of a sensible magnitude

    abundance <- abundance * (10^model$abundance.scaling)      
  
    model <- model
  })
  
  tsa.ob$param$seals <- matrix(NA, nrow = 2 + length(model$ages), ncol = 4, 
                               dimnames = list(c("cv", paste("age", model$ages), "codDep"), 
                                               c("estimate", "lower", "upper", "active")))

  tsa.ob
}


tsa.retro.fit <- function(tsa.ob, retro.year, dll = "TSA.dll", 
                          updateType = c("both", "estimates", "bounds", "neither"), ...) {
  
  updateType <- match.arg(updateType)
  
  setup <- tsa.ob$data
  
  current.year <- setup$ylim[2]
  first.year <- setup$ylim[1]
  retro.years <- current.year:retro.year
  retros <- vector("list", length(retro.years))
  names(retros) <- retro.years
  retros[[1]] <- retro <- tsa.ob
  
  for ( last.year in retro.years[-1] )
  {
    cat("\nfitting last year ", last.year, ":\n", sep = ""); flush.console()
    setup <- tsa.subset(setup, c(first.year, last.year), ...)
    
    if (updateType %in% c("both", "estimates")) 
      setup$param <- tsa.update.param(setup, retro, update = "estimates")$param
  
    newFit <-tsa.fit(setup, dll = dll)
    print(newFit$summary)
    
    if (updateType %in% c("both", "bounds")) {
      wk <- tsa.update.param(setup, newFit, update = "bounds")
      while (wk$boundUpdates) {
        cat("\nrefitting last year ", last.year, " with updated bounds:\n", sep = ""); flush.console()
        setup$param <- wk$param
        newFit <- tsa.fit(setup, dll = dll)
        print(newFit$summary)
        wk <- tsa.update.param(setup, newFit, update = "bounds")
      }
    }
    
    retros[[paste(last.year)]] <- retro <- newFit
  }
  retros
}

## plot retro summary
tsa.retro.plot.stock.summary <- function (retros, year.at, outfile) {

  require(lattice)
  require(grid)
  
  wk.summary <- do.call(rbind, lapply(retros, function(x) {
    out <- x $ stock.summary
    out $ last.year <- x $ data $ ylim[2]
    out $ year <- as.numeric(rownames(x $ stock.summary))
    out
  }))
  rownames(wk.summary) <- NULL
    
  tsa.stack <- function(data, type, label = type) {
    out <- if (type %in% c("catch", "landings", "discards")) 
      data[paste(type, c(".est", ".se", ""), sep = "")]
    else 
      data.frame(data[paste(type, c(".est", ".se"), sep = "")], rep(NA, nrow(data)))
      
    if (type %in% c("catch", "landings", "discards", "recruit")) out <- out   
      
    out <- data.frame(data $ year, data $ last.year, out, type = rep(label, nrow(data)))
    names(out) = c("year", "last.year", "response", "se", "extra", "type")
    out
  }
    
  wk.Flabel <- with(retros[[1]] $data, paste("mean F (ages ", f.range[1], "-", f.range[2], ")", sep = ""))
  data <- rbind(tsa.stack(wk.summary, "catch"), tsa.stack(wk.summary, "meanF", wk.Flabel), tsa.stack(wk.summary, "SSB"), 
                tsa.stack(wk.summary, "recruit", "recruitment"))
  
  is.discards <- retros[[1]] $data $ model $ response == "discards"
  if (is.discards) {
    data <- rbind(data, tsa.stack(wk.summary, "landings"), tsa.stack(wk.summary, "discards"))
    data <- within(data, type <- factor(as.character(type)))          # needed to get things in the correct order
  }
  
  idx <- data$type != wk.Flabel
  idy <- c("response", "se", "extra") 
  scaling <- with(retros[[1]] $ data, if (is.discards) discards $ scaling else catch $ scaling)
  data[idx, idy] <- data[idx, idy] * 10^(-scaling-3)
  
  data <- within(data, {
    lower <- response - 2 * se
    upper <- response + 2 * se
    extra[extra == 0] <- NA
  })
  
  ylim <- lapply(levels(data$type), function(i) {
    data <- subset(data, type == i)
    out <- with(data, range(lower, upper, se, extra, 0, na.rm = TRUE))
    out[1] <- max(0, out[1])
    out
  })
    
  if (missing(year.at))
    year.at <- list(major = pretty(data$year, n = 4), minor = pretty(data$year, n = 4))
  else if (!is.list(year.at)) 
    year.at <- list(major = year.at, minor = year.at)
  else if (!all(c("major", "minor") %in% names(year.at))) 
    stop('year.at incorrectly specified')
  
  wk.plot <- with(data, { 
    xyplot(
      response ~ year | type,
      ylim = ylim, as.table = TRUE,
      xlim = extendrange(year),
      scales = list(alternating = FALSE, y = list(relation = "free")),
      par.settings = list(axis.line = list(col = "transparent")),
      ylab = "", xlab = "",
      strip = FALSE, 
      between = list(x = 2, y = 4),
      axis = function(side, ...) {
        switch(
          side,
          bottom = {
            grid.lines(x = c(0, 1), y = c(0, 0), default.units = "npc", gp = gpar(col = "black", lwd = 1)) 
            panel.axis(side = side, outside = TRUE, at = year.at$minor, labels = FALSE, line.col = "black", text.cex = 1, rot = 0, tck = 0.5)
            panel.axis(side = side, outside = TRUE, at = year.at$major, labels = current.row() == 2, line.col = "black", text.cex = 1, tck = 1, rot = 0)
          },
          left = {
            grid.lines(x = c(0, 0), y = c(0, 1), default.units = "npc", gp = gpar(col = "black", lwd = 1))
            panel.axis(side = side, outside = TRUE, line.col = "black", text.cex = 1, rot = 0)
          })
      },
      panel = function(x, y, subscripts) {
        yrs <- sort(unique(last.year), decreasing = TRUE)
        col <- colorRampPalette(c("darkblue", "red"))(length(yrs))
        col1 <- paste(col, "CC", sep = "")        
        col5 <- paste(col, "55", sep = "")
        
        for (i in rev(seq(along = yrs))) {
          filt <- last.year[subscripts] == yrs[i]
          llines(x[filt], lower[subscripts][filt], col = col5[i])
          llines(x[filt], upper[subscripts][filt], col = col5[i])          
          llines(x[filt], y[filt], lwd = 2, col = col1[i])
          wm <- which.max(x[filt])
          llines(rep(max(x[filt]), 2), c(lower[subscripts][filt][wm], upper[subscripts][filt][wm]), col = col5[i])
          lpoints(max(x[filt]), y[filt][wm], pch = 1, col = col1[i])
        }
        lpoints(x[filt], extra[subscripts][filt], pch = 16, col = "black")
                             
        id <- type[subscripts][1]
        pushViewport(viewport(clip = "off"))
        grid.text(id, 0, 1.05, just = c("left", "bottom"), gp = gpar(cex = 1.0))
        upViewport()
      })
    }
  )
    
  if (!missing(outfile)) {
    win.metafile(outfile, width = 11, height = 8)
    print(wk.plot)
    #dev.off()
    invisible()
  }
  else 
    print(wk.plot)
}



tsa.subset <- function(setup, ylim, minSurveyYears = 3, earlySurvey = NULL) {
  
  setup$ylim <- ylim
  
  yrs <- as.character(ylim[1]:ylim[2])
  
  id <- intersect(c("catch", "landings", "discards"), names(setup))
  
  setup[id] <- lapply(setup[id], function(x) 
    within(x, {
      data <- data[yrs, ]
      include <- include[yrs, ] 
      cvmult <-  cvmult[yrs, ]
      weights <- weights[yrs, ]
      
      if (!is.null(earlySurvey)) include[as.character(ylim[2]), ] <- 0
    })
  )    
  
  if ("survey" %in% names(setup)) {
    
    if (is.null(earlySurvey)) 
      excludeLastYear <- rep(FALSE, length(setup$survey))
    else 
      excludeLastYear <- !(names(setup$survey) %in% earlySurvey)
    
    setup$survey[] <- mapply(function(survey, excludeLY) 
      within(survey, {
        data <- data[yrs, ]
        include <- include[yrs, ]
        cvmult <-  cvmult[yrs, ]
        model <- within(model, years <- intersect(years, as.numeric(yrs)))
        
        if (excludeLY) {
          include[as.character(ylim[2]), ] <- 0
          model <- within(model, years <- setdiff(years, ylim[2]))
        }
      }), 
      setup$survey, excludeLastYear, 
      SIMPLIFY = FALSE
    )
    
    # check sufficient years of data
    
    notOk <- sapply(setup$survey, function(x) sum(apply(x$include, 1, max))) < minSurveyYears
    
    if (all(notOk)) {
      setup$survey <- NULL
      setup$param$survey <- NULL
    } 
    else if (any(notOk)) {
      setup$survey[notOk] <- NULL
      setup$param$survey[notOk] <- NULL
    }
  }
  
  setup$auxiliary <- within(setup$auxiliary, {
    natural.mortality <- natural.mortality[yrs, ]
    maturity <- maturity[yrs, ]
    stock.weights <- stock.weights [yrs, ]
  })
  
  setup <- within(setup, {
    fsd.v.cvmult <- fsd.v.cvmult[yrs]
    fsd.y.cvmult <- fsd.y.cvmult[yrs]
    recruit.cvmult <- recruit.cvmult[yrs]
  })
  
  setup  
}


tsa.cvmult.byHaul <- function(listID, hauls, refHaul = hauls[1], years) {
  stopifnot("cvmult" %in% names(listID))
  years <- as.character(years)
  stopifnot(years %in% row.names(listID$cvmult))
  wt <- ifelse(hauls > 0, sqrt(refHaul / hauls), 0)
  within(listID, cvmult[years, ] <- apply(cvmult[years, ], 2, "*", wt))
}
