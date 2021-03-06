## PREAMBLE ###################################################################
# Script to do bootsrapping of pearson correlation

## LIBRARIES ==================================================================
library(haven)
library(pracma)
library(stringr)
library(foreign)
library(compiler)
library(psych)
library(boot)
library(lmPerm)

## FUNCTIONS ==================================================================
agk.density_p = function(x,value,type = 'smaller') {
  # function to return the probability of value
  # given the density estimated based on x
  # default (one-sided 'smaller')/'greater':
  # one-sided test, that we see value or smaller/greater
  # so useful if you wanna test difference measure "x greater y"
  # two.sided:
  # is the one-sided test, that we see abs(value) or smaller
  # so useful if you wanna test difference measure "x different y"
  if(any(is.na(x))) {return(NA)}
  if (type == 'smaller' | type == 'greater') {
    dens_diffs   = density(x,n = 2^15)
    dens_diffs$y = dens_diffs$y/sum(dens_diffs$y)
    dens_diffs$y = cumsum(dens_diffs$y)
    f <- approxfun(dens_diffs, rule=2)
    if (type == 'smaller') {
      return(f(value))
    } else {
      return(1-f(value))
    }
  } else if (type == 'two.sided') {
    if (value == 0) {
      stop('When testing against value = 0 then you probably do not want to do a two-sided test, cause it will be n.s. anyways!')
    }
    x     = abs(x)
    value = abs(value)
    return(agk.density_p(x,value,type = 'greater'))
  } else {
    stop('Wrong input for "type".')
  }
}

agk.density_p.c = cmpfun(agk.density_p)


## read in data ==================================================================
cur_file = 'data/13V_15P_mit_extrahiertenDaten_preposteinzelstats.sav'
cur_data = as.data.frame(read.spss(cur_file,use.value.labels = F))
  
## Pearson correlation ==================================================================
x = cur_data$Bac_level_Pharma
y = cur_data$desensitization_activation_39_17_8

plot(cur_data$Bac_level_Pharma, cur_data$desensitization_activation_39_17_8)
print(corr.test(x,y))
xy_data = data.frame(x,y)
xy_data = na.omit(xy_data)

## bootstrapping correlation ==================================================================
corr_function = function(xy_data,d) {
  cur_x = xy_data$x[d]
  cur_y = xy_data$y[d]
  cur_r = corr.test(cur_x,cur_y)
  return(cur_r$r)
}

cur_boot = boot(data = xy_data,statistic = corr_function,R = 30000)
agk.density_p.c(x = cur_boot$t,value = 0)

## with permutation test ==================================================================
cur_perm = lmPerm::lmp(x~y,data = xy_data,Ca = 0.0000000001,maxIter = 10000000,nCycle =1)

## Dosis und Serumlevel
d = cur_data$Menge_Hoechstdosis_Pharma
s = cur_data$Bac_level_Pharma

plot(d,s)
print(corr.test(d,s))
ds_data = data.frame(d,s)
ds_data = na.omit(ds_data)

cur_perm = lmPerm::lmp(d~s,data = ds_data,Ca = 0.0000000001,maxIter = 10000000,nCycle =1)
summary(cur_perm)

## bootstrapping correlation ==================================================================
corr_function_ds = function(xy_data,d) {
  cur_x = xy_data$d[d]
  cur_y = xy_data$s[d]
  cur_r = corr.test(cur_x,cur_y)
  return(cur_r$r)
}

cur_boot_ds = boot(data = ds_data,statistic = corr_function_ds,R = 30000)
agk.density_p.c(x = cur_boot_ds$t,value = 0)

## get data for chi-square tests with respect to abstinence ===================================
cur_file = 'data/13V_15P_mit_Alex_permutation.sav'
cur_data = as.data.frame(read.spss(cur_file,use.value.labels = F))

# variables
cur_data$abs_rel_do = factor(cur_data$abs_rel_do, levels = c(1,2),labels = c('abs','rel')) # relapser vs. abstainer
cur_data$abstinenztage_baseline_bis_t2 = as.numeric(cur_data$abstinenztage_baseline_bis_t2) # abstinenztage
cur_data$gruppe_v_p = factor(cur_data$gruppe_v_p, levels = c(1,2),labels = c('VER','PLA')) # placebo, verum

# descriptives
describeBy(cur_data$abstinenztage_baseline_bis_t2, cur_data$gruppe_v_p)
cur_tab = table(cur_data$abs_rel_do,cur_data$gruppe_v_p)
round(prop.table(cur_tab,margin = 2),digits = 2)

# test whether verum leads to more abstinence days
mod00 = lm(abstinenztage_baseline_bis_t2 ~ 1, data = cur_data)
mod01 = lm(abstinenztage_baseline_bis_t2 ~ gruppe_v_p, data = cur_data)
summary(mod01)
anova(mod00,mod01)
BIC(mod00,mod01)

# test whether verum leads to abstinent patients
mod00 = glm(abs_rel_do ~ 1, data = cur_data,family = 'binomial')
mod02 = glm(abs_rel_do ~ gruppe_v_p, data = cur_data,family = 'binomial')
summary(mod02)
anova(mod00,mod02)
BIC(mod00,mod02)


