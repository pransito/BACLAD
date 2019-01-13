## PREAMBLE ===================================================================
# running group stats for sample data BACLAD

## LIBRARIES ==================================================================
library('lmPerm')
library('foreign')
library('xlsx')

## GET DATA ===================================================================
setwd('C:\\Users\\genaucka\\Google Drive\\Promotion\\Baclad')
cur_data = read.spss('13V_15P_für_alex.sav',to.data.frame = T,reencode = 'utf-8')

## analyses ===================================================================
res_p = c()
cur_fun = function(x) return(sum(is.na(x)))
ns = table(cur_data$gruppe_v_p)
for (ii in 1:length(cur_data))  {
  # some tests
  if (names(cur_data)[ii] == 'gruppe_v_p') {
    res_p[ii] = NA
    next
  }
  # try as numeric
  cur_data[[ii]] = as.numeric(cur_data[[ii]])
  if (!is.numeric(cur_data[[ii]])) {
    res_p[ii] = NA
    next
  }
  cur_agg = aggregate(cur_data[[ii]], by=list(cur_data$gruppe_v_p),FUN=cur_fun)
  if (ns[1] == cur_agg$x[1]) {
    res_p[ii] = NA
    next
  }
  if (ns[2] == cur_agg$x[2]) {
    res_p[ii] = NA
    next
  }
  
  cur_form  = as.formula(paste(names(cur_data)[ii],'~','gruppe_v_p'))
  cur_mod   = lmp(cur_form, data=cur_data,maxIter = 1000000,Ca = 0.001,nCycle = 1)
  cur_mods  = summary(cur_mod)
  res_p[ii] = cur_mods$coefficients[2,3]
}

write.xlsx(data.frame(names(cur_data),res_p),file='permuted_p.xlsx',row.names = F)


