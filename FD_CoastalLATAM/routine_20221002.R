##################
# 0. Preamble
##################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach,doParallel, vegan, SYNCSA, tidyr,tseries,strucchange, svMisc, dplyr, mFD)
'%ni%'<- Negate('%in%')

# Load functions
source('functions/FD_df.R')

# Load data (it load a list named df_in with one dataframe with occurrences and a matrix with the proportional occurences per community/column)
load('data_input/data_in.RData')

# Make directory for output

dir.create('data_output')

rep.nulls<-3
method.nulls<- 'r00_samp'  #for other methods see ?commsim

##################
# 1. Observed metrics
##################


df<-FD_df(as.data.frame(df_in[["prop"]]), features = c(4,6,6,3,2,3,3,4))

.fdmetrics<- function(df = NULL, nPC, maxPcoa,nom.features=NULL,ord.features=NULL, weight = NULL, coord=NULL){
  # df = dataframe of relative abundances with dimensions time x ecocodes.
  # nPC = number of PC axes to include in the computation of functional diversity indices
  # maxPcoa = maximum number of axis used for PCoa
  # nom.features = indices of the features that are nominal. Default is NULL
  # ord.features = indices of the features that are ordinal Default is NULL
  # weight = matrix of abundances with dimensions time x ecocodes. If NULL, it is calculated from df.
  # coord = matrix of features/traits with dimensions ecocodes x features. If NULL, it is calculated from df.


  if(is.null(weight) && is.null(coord)){

              weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                                                              names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
              weight<-data.frame(weight, check.names = FALSE)
              rownames(weight)<-weight$time
              weight<-weight[-1]
              weight <- as.matrix(weight)
              rownames(weight)= paste0('t.', rownames(weight))

              coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
              coord<- unique(coord)
              rownames(coord)<- coord[,1]
              coord<-coord[,-1]

              for (n in 1:ncol(coord)) {
                coord[,n] = as.factor(coord[,n])
              }
        }

        coord_cat = data.frame(trait_name= colnames(coord), trait_type=NA, trait_weight=1, fuzzy_name=NA)
        coord_cat[nom.features,'trait_type'] = 'N'
        coord_cat[ord.features,'trait_type'] = 'O'


        sp_dist = mFD::funct.dist(sp_tr         = coord,
                                  tr_cat        = coord_cat,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

        fspaces_quality <- mFD::quality.fspaces(
          sp_dist             = sp_dist,
          maxdim_pcoa         = maxPcoa,
          deviation_weighting = 'absolute',
          fdist_scaling       = FALSE,
          fdendro             = 'average')

        sp_faxes_coord <- fspaces_quality$details_fspaces$sp_pc_coord

        alpha_fd_indices <- mFD::alpha.fd.multidim(
          sp_faxes_coord   = sp_faxes_coord[, 1:nPC],
          asb_sp_w         = weight,
          scaling          = TRUE,
          check_input      = TRUE,
          details_returned = TRUE)

        fd_ind_values <- alpha_fd_indices$functional_diversity_indices

        obsRao<-rao.diversity(traits=coord, weight)

        fd_ind_values$fred <-obsRao$FunRedundancy
        fd_ind_values$raoQ <-obsRao$FunRao
        fd_ind_values$simpson<- obsRao$Simpson
        rownames(fd_ind_values) =substring(rownames(fd_ind_values), 3, 10000L)
        fd_ind_values$time = rownames(fd_ind_values)


return(fd_ind_values)
}
#
obsFD<- .fdmetrics(df = df ,nPC=4, maxPcoa=10,nom.features = 1:8)


write.csv(obsFD, file='data_output/obsFD.csv')

##################
# 2. Breaks
##################
.breaks<- function(x){
  breaks<- list()

  b.TaxRich <- ts(x$sp_richn, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)
  b.raoQ <- ts(x$raoQ, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)
  b.FRic <- ts(x$fric, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)
  b.FRed <- ts(x$fred, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)
  b.FEve <- ts(x$feve, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)
  b.FSpe <- ts(x$fspe, start=c(min(as.numeric(rownames(x))), 1), end=c(max(as.numeric(rownames(x))), 1), frequency=1)

  # Breakpoints
  bp.TaxRich <- breakpoints(b.TaxRich ~ 1)
  bp.raoQ <- breakpoints(b.raoQ  ~ 1)
  bp.FRic <- breakpoints(b.FRic ~ 1)
  bp.FRed <- breakpoints(b.FRed ~ 1)
  bp.FEve<- breakpoints(b.FEve ~ 1)
  bp.FSpe <- breakpoints(b.FSpe ~ 1)

  ## confidence intervals
  ci.TaxRich <- confint(bp.TaxRich)
  ci.raoQ <- confint(bp.raoQ)
  ci.FRic <- confint(bp.FRic)
  ci.FRed <- confint(bp.FRed)
  ci.FEve <- confint(bp.FEve)
  ci.FSpe <- confint(bp.FSpe)

  breaks[['TaxRich']]<- list(bp.TaxRich,ci.TaxRich )
  breaks[['raoQ']]<- list(bp.raoQ,ci.raoQ )
  breaks[['FRic']]<- list(bp.FRic,ci.FRic )
  breaks[['FRed']]<- list(bp.FRed,ci.FRed )
  breaks[['FEve']]<- list(bp.FEve,ci.FEve )
  breaks[['FSpe']]<- list(bp.FSpe,ci.FSpe )

  return(breaks)

}

breaks<- .breaks(x = obsFD)


##################
# 3. Null models
##################

nm<-nullmodel(df_in[["count"]],method.nulls)
null<-simulate(nm, nsim =rep.nulls)

## WARNING: it may generate errors associated with convex hull. Those are removed from the final output.

.nulls<- function(null0, rep.nulls, maxPcoa, nPC, features){
  nullsdf<- data.frame()

  foreach (p=1:rep.nulls) %do% {
    print(paste0('null iteration: ',p, ' out of ',rep.nulls ))
    #progress(p,max.value = rep.nulls )
    m0<- null0[,,p]
    m1<-(t(m0)/rowSums(t(m0)))

    datam<-FD_df(as.data.frame(m1), features)



    df<-datam


    weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                          names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
    weight<-data.frame(weight, check.names = FALSE)
    rownames(weight)<-weight$time
    weight<-weight[-1]
    weight <- as.matrix(weight)
    rownames(weight)= paste0('t.', rownames(weight))

    coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
    coord<- unique(coord)
    rownames(coord)<- coord[,1]
    coord<-coord[,-1]

    for (n in 1:ncol(coord)) {
      coord[,n] = as.factor(coord[,n])
    }

    coord_cat = data.frame(trait_name= colnames(coord), trait_type='N', trait_weight=1, fuzzy_name=NA)


    sp_dist = mFD::funct.dist(sp_tr         = coord,
                              tr_cat        = coord_cat,
                              metric        = "gower",
                              scale_euclid  = "scale_center",
                              ordinal_var   = "classic",
                              weight_type   = "equal",
                              stop_if_NA    = TRUE)

    fspaces_quality <- mFD::quality.fspaces(
      sp_dist             = sp_dist,
      maxdim_pcoa         = maxPcoa,
      deviation_weighting = 'absolute',
      fdist_scaling       = FALSE,
      fdendro             = 'average')

    sp_faxes_coord <- fspaces_quality$details_fspaces$sp_pc_coord

    alpha_fd_indices <- mFD::alpha.fd.multidim(
      sp_faxes_coord   = sp_faxes_coord[, 1:nPC],
      asb_sp_w         = weight,
      scaling          = TRUE,
      check_input      = TRUE,
      details_returned = TRUE)

    fd_ind_values <- alpha_fd_indices$functional_diversity_indices

    obsRao<-rao.diversity(traits=coord, weight)

    fd_ind_values$fred <-obsRao$FunRedundancy
    fd_ind_values$raoQ <-obsRao$FunRao
    fd_ind_values$simpson<- obsRao$Simpson
    rownames(fd_ind_values) =substring(rownames(fd_ind_values), 3, 10000L)
    fd_ind_values$time = rownames(fd_ind_values)
    mdimDF <- fd_ind_values
    nullFunDiv<-data.frame(null=p,mdimDF)


    nullsdf<-rbind(nullsdf, nullFunDiv)

  }
  return(nullsdf)
}

nulldf<- .nulls(null0 = null, rep.nulls, maxPcoa = 10, nPC=4,features =  c(4,6,6,3,2,3,3,4))


.nullsummary<-function(truenull){
  nullmodelsummary<- truenull%>%
    group_by(time)%>%
    mutate(com=as.numeric(time))%>%
    summarise(msimpson= mean(simpson),
              mraoQ= mean(raoQ),
              mfred = mean(fred),
              mfric = mean(fric),
              mfeve= mean(feve),
              mfspe= mean(fspe),
              mfdis = mean(fdis),
              mfdiv= mean(fdiv),
              mfori = mean(fori),
              mTaxRich = mean(sp_richn),
              sd.simpson= sd(simpson, na.rm = TRUE),
              sd.raoQ= sd(raoQ, na.rm = TRUE),
              sd.fred= sd(fred, na.rm = TRUE),
              sd.fric= sd(fric, na.rm = TRUE),
              sd.feve= sd(feve, na.rm = TRUE),
              sd.fspe= sd(fspe, na.rm = TRUE),
              sd.fdis= sd(fdis, na.rm = TRUE),
              sd.fdiv= sd(fdiv, na.rm = TRUE),
              sd.fori= sd(fori, na.rm = TRUE),
              sd.TaxRich= sd(sp_richn, na.rm = TRUE),
              n.fric = n())%>%
    mutate(se.simpson = sd.simpson / sqrt(n.fric),
           se.raoQ = sd.raoQ / sqrt(n.fric),
           se.fred = sd.fred / sqrt(n.fric),
           se.fric = sd.fric / sqrt(n.fric),
           se.feve = sd.feve / sqrt(n.fric),
           se.fspe = sd.fspe / sqrt(n.fric),
           se.fdis = sd.fdis / sqrt(n.fric),
           se.fdiv = sd.fdiv / sqrt(n.fric),
           se.fori = sd.fori / sqrt(n.fric),
           se.TaxRich = sd.TaxRich / sqrt(n.fric),

           lower.ci.simpson = msimpson - qt(1 - (0.05 / 2), n.fric - 1) * se.simpson,
           upper.ci.simpson = msimpson + qt(1 - (0.05 / 2), n.fric - 1) * se.simpson,

           lower.ci.raoQ = mraoQ - qt(1 - (0.05 / 2), n.fric - 1) * se.raoQ,
           upper.ci.raoQ = mraoQ + qt(1 - (0.05 / 2), n.fric - 1) * se.raoQ,

           lower.ci.fred = mfred - qt(1 - (0.05 / 2), n.fric - 1) * se.fred,
           upper.ci.fred = mfred + qt(1 - (0.05 / 2), n.fric - 1) * se.fred,

           lower.ci.fric = mfric - qt(1 - (0.05 / 2), n.fric - 1) * se.fric,
           upper.ci.fric = mfric + qt(1 - (0.05 / 2), n.fric - 1) * se.fric,

           lower.ci.feve = mfeve - qt(1 - (0.05 / 2), n.fric - 1) * se.feve,
           upper.ci.feve = mfeve + qt(1 - (0.05 / 2), n.fric - 1) * se.feve,

           lower.ci.fspe = mfspe - qt(1 - (0.05 / 2), n.fric - 1) * se.fspe,
           upper.ci.fspe = mfspe + qt(1 - (0.05 / 2), n.fric - 1) * se.fspe,

           lower.ci.fdis = mfdis - qt(1 - (0.05 / 2), n.fric - 1) * se.fdis,
           upper.ci.fdis = mfdis + qt(1 - (0.05 / 2), n.fric - 1) * se.fdis,

           lower.ci.fdiv = mfdiv - qt(1 - (0.05 / 2), n.fric - 1) * se.fdiv,
           upper.ci.fdiv = mfdiv + qt(1 - (0.05 / 2), n.fric - 1) * se.fdiv,

           lower.ci.fori = mfori - qt(1 - (0.05 / 2), n.fric - 1) * se.fori,
           upper.ci.fori = mfori + qt(1 - (0.05 / 2), n.fric - 1) * se.fori,

           lower.ci.TaxRich = mTaxRich - qt(1 - (0.05 / 2), n.fric - 1) * se.TaxRich,
           upper.ci.TaxRich = mTaxRich + qt(1 - (0.05 / 2), n.fric - 1) * se.TaxRich)

  return(nullmodelsummary)
}

write.csv(nulldf, file='data_output/allnull.csv')


nullsummary<- .nullsummary(truenull=nulldf)

write.csv(nullsummary, file='data_output/nullsummary.csv')

