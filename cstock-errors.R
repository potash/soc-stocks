
## probably need to do this via error prorogation / multivariate simulation

# https://cran.r-project.org/web/packages/errors/vignettes/rjournal.pdf

library(aqp)
library(soilDB)
library(lattice)
library(errors)

library(fitdistrplus)



.getData <- function(comp, top, bottom) {
  
  # get component data for a given series / component name
  ## TODO: normalize etc.
  
  # this is slow, we only really need it for the fragment totals
  # TODO: re-write using SDA_query
  # 2x speed of SDA query
  y <- get_chorizon_from_SDA(WHERE = sprintf("compname = '%s' AND areasymbol != 'US'", comp), childs = FALSE)
  
  # init SPC
  depths(y) <- cokey ~ hzdept_r + hzdepb_r
  
  ## TODO: use fragment volume low / high
  
  # re-name fragvol_r / replace NA with 0
  y$total_frags_pct <- ifelse(is.na(y$fragvol_r), 0, y$fragvol_r)
  
  # unique instances of components
  y <- unique(y, vars = c(horizonDepths(y)))
  
  # truncate depth range
  x <- trunc(y, z1 = top, z2 = bottom)
  
  # keep required columns
  vars <- c('cokey', 'hzdept_r', 'hzdepb_r', 'om_l', 'om_r', 'om_h', 'total_frags_pct', 'dbthirdbar_l', 'dbthirdbar_r', 'dbthirdbar_h')
  
  # working with horizon data from here on out
  h <- horizons(x)[, vars]
  
  # may not be required
  # NULL -> 0
  h$total_frags_pct <- ifelse(is.na(h$total_frags_pct), 0, h$total_frags_pct)
  
  # remove NA
  idx <- which(complete.cases(h))
  h <- h[idx, ]
  
  # local copies of key variables
  id <- h$cokey
  thick <- h$hzdepb_r - h$hzdept_r
  
  # convert to approximate soil organic carbon
  # encode as fraction
  soc_l <- (h$om_l / 100) * 0.58
  soc_h <- (h$om_h / 100) * 0.58
  soc <- (h$om_r / 100) * 0.58
  
  # soil fraction based on total coarse fragment volume
  fraction <- (1 - (h$total_frags_pct / 100))
  
  # Db
  db_l <- h$dbthirdbar_l
  db_h <- h$dbthirdbar_h
  db <- h$dbthirdbar_r
  
  # compile into a DF and done
  d <- data.frame(
    id,
    thick,
    soc,
    fraction,
    db,
    soc_l,
    soc_h,
    db_l,
    db_h
  )
  
  return(d)
}



## bootstrapping approach
.Cstock.boot <- function(d, n.sim) {
  
  # sample with replacement of entire profiles
  ids <- sample(unique(d$id), size = n.sim, replace = TRUE)
  d.boot <- lapply(ids, function(i) {
    
    # compute SOC stock for this profile
    d.i <- d[d$id == i, ]
    res <- sum(
      with(d.i, thick * soc * fraction * db * 10), 
      na.rm = TRUE
    )
    
    return(res)
  })
  
  s <- do.call('c', d.boot)
  
  # 5th-50th-95th percentiles
  s.q <- quantile(s, probs = c(0.05, 0.5, 0.95))
  return(s.q)
}

## error propagation approach (symmetric errors)
.Cstock.errors <- function(d) {
  
  # constant errors applied to thickness and fine earth fraction
  errors(d$thick) <- 5
  errors(d$fraction) <- 0.05
  
  # use 1/2 l-H range for now
  errors(d$soc) <- (d$soc_h - d$soc_l) / 2
  errors(d$db) <- (d$db_h - d$db_l) / 2
  
  # double-check
  # format(d$soc[1], digits = 2, notation = 'plus-minus')
  
  # graphical check on errors
  plot(d$soc, d$db, las = 1)
  
  # compute stock
  d$stock <- d$thick * d$soc * d$fraction * d$db * 10
  
  # sum by profile
  s <- tapply(d$stock, d$id, sum, simplify = FALSE)
  s <- do.call('c', s)
  
  # keep median
  # kg / m^2
  return(median(s))
}



SOC.stock <- function(comp, top = 0, bottom = 100, n.sim = 1000) {
  
  d <- .getData(comp, top = top, bottom = bottom)
  
  
  CS.errors <- .Cstock.errors(d)
  CS.boot <- signif(.Cstock.boot(d, n.sim = n.sim), 2)
  
  # CS.errors.txt <- format(CS.errors, digits = 2, notation = 'plus-minus')
  CS.errors.txt <- sprintf("%s ± %s kg/m^2", signif(CS.errors, 2), signif(errors(CS.errors), 2))
  
  CS.boot.txt <- sprintf("%s [%s\u2014%s] kg/m^2", CS.boot[2], CS.boot[1], CS.boot[3])
  
  res <- list(
    depth.range = c(top, bottom),
    n.comp = length(unique(d$id)),
    stock.errors = CS.errors.txt,
    stock.boot = CS.boot.txt
  )
  
  return(res)
  
}


SOC.stock('Zook')

SOC.stock('Zook', n.sim = 1)

SOC.stock('Pierre')

SOC.stock('Drummer')


SOC.stock('Lucy')

SOC.stock('Sierra')

SOC.stock('Fresno')

SOC.stock('Amador')

SOC.stock('Leon')

SOC.stock('Rindge')

SOC.stock('Clarksville')

SOC.stock('Pentz')

SOC.stock('Ryde')




### still working on this

## fitting distributions
# https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
# https://www.r-project.org/conferences/useR-2009/slides/Delignette-Muller+Pouillot+Denis.pdf

s <- c(soc_l, rep(soc, times = 3), soc_h)
descdist(s, boot = 500)

s.fit <- fitdist(s[s>0], distr = 'lnorm')
plot(s.fit)

Db <- c(db_l, rep(db, times = 3), db_h)
descdist(Db, boot = 500)

db.fit <- fitdist(Db, distr = 'norm')
plot(db.fit)

## this isn't quite right
# descdist(thick)
# thick.fit <- fitdist(thick, distr = 'norm')
# plot(thick.fit)

descdist(fraction, boot = 500)
fraction.fit <- fitdist(fraction, distr = 'lnorm')
plot(fraction.fit)










