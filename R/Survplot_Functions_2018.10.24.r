#Jenny Smith

#march 30, 2017
#Purpose: Create a list of survival analyses and Given survival package Survfit() object, plot the survival curves and calculate the P values for differences in survival.


#' Generate a list of KM and Cox PH objects
#'
#' @param df dataframe with clinical data.
#' @param colNames character vector (length 2) with the colnames of the time and event, in that order.
#' @param group column name as character vector of the explanatory variable to test, or numeric value of 1 for 1 line KM
#' @param rho 1 or 0 or inbetween. 0=logrank test and 1=gehan-wilcoxon test
#' @param time optional if you need to convert days to years for example.
#' @param ref character indicating the group classification that is the reference
#'
#' @returns a list with survfit objects
#' @export
#'
#' @examples
#' TBD <- c("")
SurvObjects <- function(df, colNames, group, rho=0, time=NULL,ref=NULL){
  #df is the dataframe with clinical data.
  #colnames is a character vector (length 2) with the colnames of the time and event, in that order.
  #group is the column name as character vector of the explanatory variable to test
  #time is optional if you need to convert days to years for example.
  #rho is 1 or 0 or inbetween. 0=logrank test and 1=gehan-wilcoxon test

  # require(survival)
  # require(survMisc)

  #optional time unit conversion
  if (is.null(time)){
    time <- (df[[colNames[1]]])
  }else if (time == "DtoY"){
    time <- (df[,colNames[1]]/365)
  }else if (time == "MtoY"){
    time <- (df[,colNames[1]]/12)
  }

  #Kaplan Meier Estimates
  KM <- Surv(time = time, event = df[[colNames[2]]])

  if (group == 1){
    #model fit
    survFit <- survfit(KM ~ 1)
    return(survFit)

  }else{

    if (!is.null(ref)){
      df[,group] <- relevel(as.factor(df[,group]), ref=ref)
    }

    #model fit
    survFit <- survfit(KM ~ df[,group], data=df)

    #cox proportional Hazards
    cph.test <- coxph(KM ~ df[,group])

    #test proportional hazards assumption
    cph.zph <- cox.zph(cph.test, transform = "km")
    rhoPval <- as.data.frame(cph.zph$table)$p #test the PH assumption

    if (length(rhoPval) > 1){
      rhoPval <- rhoPval[length(rhoPval)]
    }

    # survMisc compare stat tests for two curves
    tne <- ten(eval(survFit$call$formula), data=df)
    capture.output(comp(tne)) -> allTests

    if (rhoPval >= 0.05 ){
      attr(tne, "lrt")[1,] -> pVal #log-rank test
      test <- "log.rank"
    }else if (rhoPval < 0.05){
      attr(tne, "lrt")[2, ] -> pVal #Gehan-Breslow-Wilcoxon
      test <- "Gehan.Breslow.Wilcoxon"
    }

    list <- list(survFit, pVal, cph.test, cph.zph, KM)
    names(list) <- c("survFit", test, "CoxPH", "PH_Test","KaplanMeirerEst")
    return(list)
  }
}


#
#' Survival plot without gridlines
#'
#' @param fit the output of the survFit() function.
#' @param LegendTitle character vector which will be the plot legend title.
#' @param timeUnit character vector for the follow-up time in days, years, months etc.
#' @param include.censored boolean whether the add the ticks to indicate censored
#' @param colors character vector the same length as the number of groups.
#' @param pval string (charcter or numeric) of the p-value
#' @param max.year time (singe numeric) in years to end the plot at.
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#'  TBD <- ("")
SurvivalPlot <- function(fit, LegendTitle, timeUnit,include.censored=TRUE, colors,pval=NULL,max.year=NULL){
  #fit the output of the survFit() function.
  #legend title is a character vector which will be the plot legend title.
  #timeUnit is a character vector for the follow-up time in days, years, months etc.
  #colors is a character vector the same length as the number of groups.
  #max.year is the time (singe numeric) in years to end the plot at.
  #pval is string (charcter or numeric) of the p-value

  # require(survival)
  # library(ggplot2)
  # library(GGally)

  if(is.null(pval)){
    pval <- ""
  }else{
    pval <- paste0("p = ",pval)
  }

  if(!is.null(max.year)){
    # if(max(fit$time) < max.year){
    pos.x <- max.year
    # }
  }else{
    pos.x <- max(fit$time)
  }

  if(is.null(fit$strata)){
    lm <- 1 #left margin (lm)
    rm <- 4  #right margin(rm)
    num.cols <- 1
    num.rows <- 1

    #main survival plot
    p <- ggsurv(fit, surv.col = colors,
                cens.col=colors, CI=FALSE,
                lty.est = 1, size.est = 1.0,
                cens.size = 2.0,order.legend=FALSE)

  }else{

    m <- max(nchar(gsub("^.+\\=(.+)", "\\1", names(fit$strata)))) #maximum number of characters
    d <- length(fit$strata) #how many groups
    lm <- case_when(m <= 2 ~ 0.90,
                    m <= 5 ~ 1.25,
                    m > 5 & m <= 9 ~ 1.7,
                    m > 9 & m < 13 ~ 2.6,
                    m >= 13 & m < 17 ~ 3.4,
                    m >= 17 & m < 21 ~ 5.0,
                    m >= 21 & m < 23 ~ 6.0,
                    m >= 23 & m < 25 ~ 6.5,
                    m >= 25 ~ 7.5)#left margin (lm)
    # print(c("km.plot", lm))
    rm <- case_when(d <= 4 ~ 1.5, d > 4 ~ 2.0) #right margin(rm)
    num.cols <- 2
    num.rows <- 2

    if(d > 4 & d < 7){
      num.rows <- 3
    }else if(d > 7 & d <= 18){
      num.cols <- 3
      num.rows <- 6
    }else if(d > 18){
      num.cols <- 5
      num.rows <- 5
    }

    if(num.cols*num.rows < d){
      num.rows <- 10
    }

    #main survival plot
    p <- ggsurv(fit, plot.cens=include.censored,
                surv.col = colors, CI=FALSE,
                lty.est = 1, size.est = 1.0,
                cens.size = 2.0,order.legend=FALSE)
  }

  #customize the plot
  p <- p +
    labs(y= "Fraction Surviving", x = paste(timeUnit, "Follow-up")) +
    scale_y_continuous(limits = c(0,1),
                       breaks = seq(0,1,0.2),
                       minor_breaks = NULL) +
    scale_x_continuous(limits = c(0,pos.x),
                       breaks = seq(0,pos.x, 2)) +

    theme(plot.title = element_text(hjust = 0.5, size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill= NA, size=1.0),

          axis.text = element_text(colour = "black"),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=12),
          axis.title = element_text(size = 12),
          legend.position = "top",
          legend.direction  = "horizontal",
          legend.margin = margin(0,0,0,0, unit = "pt"),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          legend.key=element_blank()) +

    theme(plot.margin = margin(0, rm,.5,lm,unit="cm")) +
    annotate(geom="text", x=1, y=0.05, label=pval, size=5) +
    guides(color=guide_legend(ncol=num.cols, nrow=num.rows))

  if (length(colors) > 1){
    p$scales$scales[[1]]$name <- LegendTitle
    p$scales$scales[[2]]$name <- LegendTitle
  }

  return(p)
}


#risk table dataframe
risk.table <- function(survFit.obj, f.levels=NULL,
                       times=NULL,col="black",
                       ret.data=FALSE, max.year=NULL){
  #survFit.obj is the results  survfit()
  #f.levels is an optional character vector to relevel the order of the groups, if desired. else levels is alphabetical
  #times would be a seperate numeric vector, with years/x-axis breaks points, if desired.

  if(!is.null(max.year)){
    # if (max(survFit.obj$time) < max.year){
      times <- seq(0,max.year, by= 2)
      xlims <- max.year
    # }
  }else{
    times <- seq(0,max(survFit.obj$time), by= 2)
    xlims <- max(survFit.obj$time)
  }

  #Summary dataframe with the selected time points to get % OS/ EFS for plot x-axis
  summ <- summary(survFit.obj, times=times, extend=TRUE)


  if(is.null(f.levels)){
    f.levels <- unique(as.character(summ$strata)) %>%
      gsub("^.+\\=(.+)", "\\1", .)
  }

  #Create a new dataframe to hold the time,KM estimate and strata for plotting
  risk.df <- bind_cols(list(time=summ$time,
                            surv=summ$surv,
                            n.risk=summ$n.risk,
                            Group=as.character(summ$strata))) %>%
    set_colnames(c("time","surv","N.Risk","Group")) %>%
    mutate(Group=gsub("^.+\\=(.+)", "\\1", Group),
           Group=factor(Group, levels = f.levels))


  if(ret.data){
    return(risk.df)
  }

  #variable for deciding length of margins. longer group names would need more space.
  m <- max(nchar(levels(risk.df$Group))) # print(c("length:", m))
  d <- length(levels(risk.df$Group))
  lm <- case_when(m <= 5 ~ 1.75,
                  m > 5 & m <= 9 ~ 1.00,
                  m > 9 & m <= 10 ~ 1.00,
                  m > 10 & m < 13 ~ 1.10,
                  m >= 13 ~ 1.20)

  rm <- case_when(d <= 5 ~ 1.5,
                  m > 5 ~ 2.0)


  #gglot the risk table (risk.df)
  tbl <- ggplot(risk.df, aes(x = time, y = Group, label=N.Risk)) +
    geom_text(size = 5)  +
    theme_bw()  +
    scale_x_continuous("Number at risk",
                       limits = c(0,xlims),
                       breaks = times) +
    xlab(NULL) +
    ylab(NULL) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position="none",
          axis.line = element_line(size=0.7, color="black"),
          axis.title.x = element_text(size=12),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face="bold",
                                     size=12,
                                     color= col[levels(risk.df$Group)],
                                     hjust=1)) +
    #of note: this margin is still not perfect..... the size of window in Rstudio is different than what actually shows up on the tiff...
    theme(plot.margin = margin(-1.5, rm, 0.1, lm, unit="cm"))




  return(tbl)
}



#' Create Kaplan-Meier Plots
#'
#' @param df cleaned CDE with patient IDs as rownames.
#' @param group_vars character vector of a column name to group_by() in dplyr.
#' @param column_names list with colnames OS for OS and EFS for "EFS" in the CDE
#' @param covariate character string for the column in CDE to be the explanatory variable
#' @param cc color codes with one color per group
#' @param riskTable boolean to generate a risk table under the KM plot
#' @param f.levels character vector of factor levels for the risk table
#' @param max.year year to calculate the statistics in the risk table (eg 5yrs after follow-up)
#' @param include.censored boolean to draw the tick marks on the curve to show censored points
#'
#' @returns a grouped data.frame with the surv objects and the KM plots
#' @export
#'
#' @examples
#' my_cols <-  list("OS"=c("OS.time","OS.id"), "EFS"=c("EFS.time","EFS.id"))
#'
KM.plots <- function(df, group_vars=NULL,
                     column_names=list("OS"=c("Overall.Survival.Time.in.Days", "OS.ID")),
                     covariate,
                     cc=NULL,
                     riskTable=TRUE,
                     f.levels=NULL,
                     max.year=NULL,
                     include.censored=TRUE){


  #Define KM estimates for survival analysis
  # if (cohort == "0531" & type=="OS"){
  #   os.est <- "Surv(Overall.Survival.Time.in.Days/365.25, OS.ID)"
  #   efs.est <- "Surv(Event.Free.Survival.Time.in.Days/365.25, Event.ID)"
  # }else if (cohort == "1031" & type=="OS"){
  #   os.est <- "Surv(OS.time..days./365.25, OS.ID)"
  #   efs.est <- "Surv(EFS.time..days./365.25, Event.ID)"
  # }else if (cohort == "TCGA" & type=="OS"){
  #   os.est <- "Surv(OS.months..3.31.12/12, clinData1.vital_status)"
  #   efs.est <- "Surv(EFS.months....3.31.12/12, first.event)"
  # }else if (cohort=="BEAT_AML" & type=="OS"){
  #   os.est <- "Surv(OverallSurvival_VIZOME/365.25, OS.ID)"
  #   #No EFS data available
  # }else if (type=="EOI"){

  os.est <- glue::glue("Surv({column_names$OS[1]}, {column_names$OS[2]})")

  if(!is.null(column_names[["EFS"]])){
    efs.est <- glue::glue("Surv({column_names$EFS[1]}, {column_names$EFS[2]})")
  }

  #function for colors
  colorcodes <- function(df, covariate){
    strata <- unlist(unique(df[,covariate])) #Unique groups (strata)
    strata <- strata[order(strata)] #ensure the alphabetical order
    len <- length(strata)

    if ( ! is.null(cc)){
      return(cc)

    }else{
      colors <- c("dodgerblue4", "brown3", "darkgoldenrod3",
                  "deepskyblue1",
                  "azure4", "darkmagenta",
                  "turquoise4", "green4",
                  "deeppink", "navajowhite2","chartreuse2",
                  "lightcoral",
                  "mediumorchid", "saddlebrown", "#466791",
                  "#60bf37","#953ada",
                  "#4fbe6c","#ce49d3","#a7b43d","#5a51dc",
                  "#d49f36","#552095","#507f2d","#db37aa",
                  "#84b67c","#a06fda","#df462a","#5b83db",
                  "#c76c2d","#4f49a3","#82702d","#dd6bbb",
                  "#334c22","#d83979","#55baad","#dc4555",
                  "#62aad3","#8c3025","#417d61","#862977",
                  "#bba672","#403367","#da8a6d","#a79cd4",
                  "#71482c","#c689d0","#6b2940","#d593a7",
                  "#895c8b","#bd5975")
      cc <- NULL
      for ( i in 1:len){
        c <- colors[i]
        cc <- c(cc,c)
      }
      names(cc) <- strata
      return(cc)
    }
  }


  #Formula for KM curves
  OS.form <- as.formula(paste(os.est,' ~ ',covariate))

  #perform survival analysis for all factor levels
  grouped.df <-  df %>%
      group_by(!!! rlang::syms(group_vars)) %>%
      do(OS.cox=coxph(OS.form, data = .),
         OS.fit=survfit(OS.form, data=.),
         OS.diff=survdiff(OS.form, data=.), #survdiff for log-rank p-value
         OS=SurvivalPlot(survfit(OS.form, data=.), #survival plots
                         LegendTitle = covariate,
                         timeUnit = "Years",
                         max.year = max.year,
                         include.censored=include.censored,
                         colors =  colorcodes(.,covariate = covariate)))

    if(exists("efs.est")){
      #Formula for EFS KM curves
      EFS.form <- as.formula(paste(efs.est,' ~ ', covariate))

      grouped.df.efs <-  df %>%
        group_by(!!! rlang::syms(group_vars)) %>%
        do( EFS.cox=coxph(EFS.form, data = .),
            EFS.fit=survfit(EFS.form, data=.),
            EFS.diff=survdiff(EFS.form, data=.), #survdiff for log-rank p-value
            EFS=DeGSEA::SurvivalPlot(survfit(EFS.form, data=.),
                             LegendTitle= covariate,
                             timeUnit= "Years",
                             max.year = max.year,
                             include.censored=include.censored,
                             colors= colorcodes(., covariate = covariate)))

      grouped.df <- bind_cols(grouped.df, grouped.df.efs)
    }

    #Denote the grouping variable in the results dataframe
    if(is.null(group_vars)){
      grouped.df <- add_column(grouped.df, Group="AML",.before = 1)
    }else{
      col.idx <- which(colnames(grouped.df)=="OS.cox")-1
      grouped.df <- unite(grouped.df,
                          "Group",1:all_of(col.idx), remove = FALSE, sep=" in ")
    }


  #function to calculate a p-value
  pvalue <- function(survdiff.res){
    p <- pchisq(survdiff.res$chisq, df=length(survdiff.res$n) - 1, lower=FALSE)
    p <- ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, digits = 3)))
    return(p)
  }

  #Function to extract N of cohort.
  groupSize <- function(survdiff.res){paste0("N = ", sum(survdiff.res$n))}

  #Update the survival plots with an informative title and p values, and group N
  types <- c("Relapse", "Failure", "OS", "EFS")
  plots <- which(colnames(grouped.df) %in% types)

  #for each row (nrow > 1, when grouping variable given)
  for ( i in 1:nrow(grouped.df)){

    #for each OS/EFS column
    for (plot.idx in plots){
      fit.idx <- plot.idx-2
      fit <- grouped.df[[fit.idx]][[i]]
      summ <- summary(fit)
      median.surv <- summ$table[,"median"]
      x_pos <- ifelse(min(median.surv,na.rm=T) < 1, max(summ$time,na.rm = T), 1)
      y_pos <- ifelse(min(median.surv, na.rm=T) < 1, 0.95, 0.05)

      idx.logrank <- plot.idx-1
      N <- groupSize(grouped.df[[idx.logrank]][[i]])
      pval <- pvalue(grouped.df[[idx.logrank]][[i]]) %>%
        paste(.,N, sep="\n")

      #Add Main Title that is the column name (eg Relapse,OS, EFS, etc) and grouping factors
      grouped.df[[plot.idx]][[i]]$labels$title <- paste(names(grouped.df[plot.idx]),
                                                        grouped.df[["Group"]][i],
                                                        sep = ": ")

      #Add P value and the N
      grouped.df[[plot.idx]][[i]] <- grouped.df[[plot.idx]][[i]] +
        annotate(geom="text", x=x_pos, y=y_pos, label=pval, size=5) +
        theme(plot.title = element_text(face="bold"))

      #Add a risk table (Number of patients at risk at each time )
      if(riskTable){

        risk.tab <- risk.table(survFit.obj = fit,
                               ret.data = TRUE)
        risk.tab.p <- risk.table(survFit.obj = fit,
                                 col = cc,
                                 f.levels = f.levels,
                                 max.year=max.year)

        # Create a blank plot for place-holding
        blank.pic <- ggplot(risk.tab, aes(time, surv)) +
          geom_blank() +
          theme_bw() +
          theme(axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank())


      # grid.arrange(grouped.df[[plot.idx]][[i]], blank.pic, risk.tab.p, clip = FALSE, nrow = 3,
      #                    ncol = 1, heights = unit(c(2, .1, .25), c("null", "null", "null")))

       grouped.df[[plot.idx]][[i]] <-  arrangeGrob(grouped.df[[plot.idx]][[i]],
                                                   blank.pic, risk.tab.p,
                                                   clip = FALSE,
                                                   nrow = 3, ncol = 1,
                                                   heights = unit(c(2, .2, .4),
                                                                  c("null", "null", "null")))

      }
    }
  }

  return(grouped.df)
}


#Function to save the outputs from KM.plots()
saveMultiPlots <- function(KM.plots.res,w=8,h=5){
  #This is for the relapse results from KM.plots()

  N <- nrow(KM.plots.res)

  covar <- function(i){
    names(KM.plots.res[grep("diff",colnames(KM.plots.res))][[i,1]]$n) %>%
      str_split_fixed(.,pattern="=",n=2) %>% .[,1] %>% unique()
  }

  name <- function(i){paste(KM.plots.res[[1]][i],covar(i),col,"KMplot.tiff", sep="_")}
  cols <- grep("diff",colnames(KM.plots.res), value = TRUE, invert = TRUE)[-1]

  for (col in cols){
    lapply(1:N, function(x) ggsave(filename = name(x),
                                   plot = KM.plots.res[[col]][[x]],
                                   device = "tiff",
                                   width = w,
                                   height = 5,
                                   dpi=600))

  }
}




#Added from LSC17  Analysis  on July 12, 2017.
weibull.HR <- function(weibull.mod){
  #weibull.mod is the output of survreg(, dist="wiebull")
  # multivar is a numeric vector idicating the indices of multi varaite analysis to be returned.
  if (is.null(multivar)){
    idx <- 2 #univariate analysis
  }else{
    idx <- 2:length(weibull.mod$coefficients)
  }
  HR <- round(exp(weibull.mod$coefficients*-1*1/weibull.mod$scale)[idx], digits = 3)
}



#Added from LSC17  Analysis  on July 12, 2017.
parametricSumaryTable <- function(survreg.res){
  summ <- summary(survreg.res)
  table <- as.data.frame(summ$table)

  # coef <- summ$coefficients[which(! grepl("Int", names(summ$coefficients)))]
  idx <- which( ! grepl("Int|Log", names(survreg.res$coefficients)))

  UCI_HR2 <- round(exp((survreg.res$coefficients[idx] - 1.96*table$`Std. Error`[idx])*-1*1/survreg.res$scale),
                   digits = 3)
  LCI_HR2 <- round(exp((survreg.res$coefficients[idx]  + 1.96*table$`Std. Error`[idx])*-1*1/survreg.res$scale),
                   digits=3)


  HR <- round(exp(survreg.res$coefficients[idx]*-1*1/survreg.res$scale), digits = 3)
  CI <- paste(LCI_HR2, UCI_HR2, sep="-")
  pVal <- round(table$p[idx], digits=3)


  stats <- data.frame(HR, CI,pVal)
  colnames(stats) <- c("Hazard Ratio", "95% Confidence Interval", "P-value")
  return(stats)
}


#Added from LSC17 Multivariate Analysis  on June 24th, 2017.
coxSummaryTable <- function(coxph.res,Colname="Group"){
  c.mod <- coxph.res
  c.summ <- summary(coxph.res)
  c.table <- as.data.frame(c.summ$coefficients) %>% cbind(.,as.data.frame(c.summ$conf.int))
  # c.table
  c.HR <- round(c.table$`exp(coef)`, digits=3)
  c.CI <- paste(round(c.table$`lower .95`,digits=3),round(c.table$`upper .95`, digits = 3), sep="-")
  c.pVal <- round(c.table$`Pr(>|z|)`, digits=3)

  comp <- gsub("^.+(Yes|No|Unknown|High|Low)", "\\1",names(c.mod$coefficients))
  c.stats <- data.frame(names(c.mod$coefficients) ,comp, c.HR, c.CI, c.pVal)
  colnames(c.stats) <- c("Variable", Colname,"Hazard Ratio", "95% Confidence Interval", "P-value")
  # rownames(c.stats) <- names(c.mod$coefficients)

  return(c.stats)
}



# functions to calculate the p-value for differences in KM curves from survdiff() function
#
#' Tests if there is a difference between two or more survival curves using the G-rho family of tests, or for a single curve against a known alternative.
#'
#' @param survdiff object generated by the survdiff() function from survival package
#' @param digits number of digits in the p-value to report
#'
#' @returns a p-value
#' @export
#'
#' @examples
#' TBD <- c("")
calc_KMcurve_pvalues <- function(survdiff, digits=3){
  p <-  pchisq(survdiff$chisq, length(survdiff$n)-1,
              lower.tail = FALSE) %>%
    round(.,  digits = digits)

  return(p)
}



#
#uses the summary() function on the survfit() object
#
#' create a table with OS/EFS percent at specified time point
#'
#' @param fit the output of the survFit() function
#' @param time the time-point when to calculate the statistics (eg 5 yrs after follow-up)
#' @param pvalues optionally input p-values manually from survdiff() function
#'
#' @returns a data.frame with the outcome statistics
#' @export
#'
#' @examples
#' TBD <- c("")
outcome_table <- function(fit, time=5, pvalues=NULL){

  summ <- summary(fit, times=time)
  type <- ifelse(any(grepl("OS", summ$call)),"OS","EFS")
  summ <- summ[!grepl("table|std.chaz|call", names(summ))]
  summ <- summ[sapply(summ,length)==length(summ$strata)]

  surv_col <- paste("Percent", type)
  n <- length(summ$strata)

  tab <- tibble(as.data.frame(summ)) %>%
    mutate_at(vars(surv,lower,upper,std.err), ~round(.*100, digits = 2)) %>%
    mutate(Group=str_split_fixed(strata, pattern = "=", n=2)[,2])


  if(!is.null(pvalues) & length(pvalues) < n){
    pval_ref <- NA #ref is the first position in the strata
    pvalues=c(pval_ref, pvalues) #must be in same order as strata in the survift()!
    tab <- tab %>%
      add_column(pvalue=pvalues)

    tab <- tab %>%
      select(Group,
             "Number Patients per Group"=n,
             "Time (years)"=time,
             !!surv_col := surv,
             "Lower 95% CI"=lower,
             "Upper 95% CI"=upper,
             "Standard Error"=std.err,
             "p-value"=pvalue)
  } else{

    tab <- tab %>%
      select(Group,
             "Number Patients per Group"=n,
             "Time (years)"=time,
             !!surv_col := surv,
             "Lower 95% CI"=lower,
             "Upper 95% CI"=upper,
             "Standard Error"=std.err)

  }


  return(tab)
}





# introduction This entire function was written by Lindsay Keegan based off
# a tutorial by Edwin Thoen posted here:
# http://www.r-statistics.com/2013/07/creating-good-looking-survival-curves-the-ggsurv-function/
# and slightly adapted by Avery McIntosh, 2015.
#See TCGA AML 05May2017.Rmd for example usage.
#S is the fit object
#NEEDS TO BE TESTED AGAIN

ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                     cens.col = 'red', lty.est = 1, lty.ci = 2,
                     cens.shape = " ", back.white = F, xlab = 'Time (days)',
                     ylab = 'Retention Probability', main = '') {


  strata <- length(s$strata)
  n <- s$strata

  #update JS to deal with the super imposed curves to seq by 4 for 24 and split with comma
  groups <- factor(unlist(lapply(strsplit(names(s$strata),split = "=|,"),
                                 function(x) paste(x[2], x[4], sep="."))))

  gr.name <-  unlist(strsplit(names(s$strata), '='))[1] #LSC17
  gr.df <- vector('list', strata)
  ind <- vector('list', strata)
  n.ind <- c(0,n); n.ind <- cumsum(n.ind)

  for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]

  for(i in 1:strata){
    gr.df[[i]] <- data.frame(
      time = c(0, s$time[ ind[[i]] ]),
      surv = c(1, s$surv[ ind[[i]] ]),
      up = c(1, s$upper[ ind[[i]] ]),
      low = c(1, s$lower[ ind[[i]] ]),
      cens = c(0, s$n.censor[ ind[[i]] ]),
      group = rep(groups[i], n[i] + 1))
  }

  dat <- do.call(rbind, gr.df)
  dat.cens <- subset(dat, cens != 0)

  pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    geom_step(aes(col = group),lty=1, size=1.15) +
    geom_point(data = dat.cens, aes(x=time, y=surv, color=group), shape=3, size=4) +
    labs(y= "Fraction Surviving", x = paste("Follow-up in", "Years")) +
    theme(plot.title = element_text(hjust = 0.5, size=18),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size=16),
          axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "dark grey", fill= NA, size=1.0),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), minor_breaks = NULL) +
    scale_x_continuous(limits = c(0,max(dat$time)), breaks = seq(0,max(dat$time), 2))

  # list <- list(n, groups, gr.name)
  # names(list) <- c("n", "groups", "gr.name")
  return(p1)

}

