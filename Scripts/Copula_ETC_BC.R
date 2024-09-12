




#' Calculates the isoline and relative probability of events on the isoline, given the observational data, for one or more user-specified return periods. Outputs the single "most-likely" design event or an ensemble of possible design events obtained by sampling along the isoline according to these relative probabilities. The design event under the assumption of full dependence is also computed.
#' @param Data Data frame of dimension \code{nx2} containing two co-occurring time series of length \code{n}.
#' @param Data_Con1 Data frame containing the conditional sample (last 30 years) (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column.
#' @param Data_Con2 Data frame containing the conditional sample (last 30 years)(declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column. Can be obtained using the \code{Con_Sampling_2D} function.
#' @param Data_Con1_M Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column.
#' @param Data_Con2_M Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column. Can be obtained using the \code{Con_Sampling_2D} function.
#' @param u1 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the first column was sampled in \code{Data_Con1}.
#' @param u2 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the second column was sampled in \code{Data_Con2}.
#' @param Thres1 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in \code{Data_Con1}. Only one of \code{u1} and \code{Thres1} should be supplied. Default is \code{NA}.
#' @param Thres2 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in \code{Data_Con2}. Only one of \code{u2} and \code{Thres2} should be supplied. Default is \code{NA}.
#' @param Copula_Family1 Numeric vector of length one specifying the copula family used to model the \code{Data_Con1} dataset.
#' @param Copula_Family2 Numeric vector of length one specifying the copula family used to model the \code{Data_Con2} dataset. Best fitting of 40 copulas can be found using the \code{Copula_Threshold_2D} function.
#' @param Marginal_Dist1 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con1}.
#' @param Marginal_Dist2 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con2}.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data}.
#' @param GPD1 Output of \code{GPD_Fit} applied to variable \code{con1} i.e., GPD fit \code{con1}. Default \code{NA}. Only one of \code{u1}, \code{Thres1} and \code{GPD1} is required.
#' @param GPD2 Output of \code{GPD_Fit} applied to variable \code{con2} i.e., GPD fit \code{con2}. Default \code{NA}. Only one of \code{u2}, \code{Thres2} and \code{GPD2} is required.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{Data}. Default is \code{365.25}, daily data.
#' @param GPD_Bayes Logical; indicating whether to use a Bayesian approach to estimate GPD parameters. This involves applying a penalty to the likelihood to aid in the stability of the optimization procedure. Default is \code{FALSE}.
#' @param RP Numeric vector specifying the return periods of interest.
#' @param Decimal_Palace Numeric vector specifying the number of decimal places to which to specifiy the isoline. Defulat is \code{2}.
#' @param Interval Numeric vector specifying the number of equally spaced points comprising the combined isoline.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @param N Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an isoline. Samples are collected from the two joint distribution with proportions consistent with the total number of extreme events conditioned on each variable. Default is \code{10^6}
#' @param N_Ensemble Numeric vector of length one specifying the number of possible design events sampled along the isoline of interest.
#' @param Sim_Max Numeric vector of length one specifying the maximum value, given as a multiple of the largest observation of each variable, permitted in the sample used to estimate the (relative) probabilities along the isoline.
#' @param Plot_Quantile_Isoline Logical; indicating whether to first plot the quantile isoline. Default is \code{FALSE}.
#' @param Isoline_Type Character vector of length one specifying the type of isoline. For isolines obtained using the overlaying method in Bender et al. (2016) use \code{"Combined"} (default). For quantile isoline from the sample conditioned on variable \code{Con1}|(\code{Con2}) use \code{"Con1"}(\code{"Con2"}).
#' @return Plot of all the observations (grey circles) as well as the declustered excesses above Thres1 (blue circles) or Thres2 (blue circles), observations may belong to both conditional samples. Also shown is the isoline associated with \code{RP} contoured according to their relative probability of occurrence on the basis of the sample from the two joint distributions, the "most likely" design event (black diamond), and design event under the assumption of full dependence (black triangle) are also shown in the plot. The function also returns a list comprising the design events assuming full dependence \code{"FullDependence"}, as well as once the dependence between the variables is accounted for the "Most likley" \code{"MostLikelyEvent"} as well as an \code{"Ensemble"} of possible design events and relative probabilities of events on the isoline \code{Contour}. The quantile isolines with \code{Quantile_Isoline_1} and \code{Quantile_Isoline_2}, and GPD thresholds with \code{Threshold_1} and \code{Threshold_2}.
#' @seealso \code{\link{Copula_Threshold_2D}} \code{\link{Diag_Non_Con}} \code{\link{Diag_Non_Con_Trunc}}
#' @export
#' @examples

Copula_ETC_BC<-function(Data, Data_Con1, Data_Con2,Data_Con1_M, Data_Con2_M, u1=NA, u2=NA, Thres1, Thres2, Copula_Family1, Copula_Family2, Marginal_Dist1, Marginal_Dist2, Marginal_Dist1_Par=NA, Marginal_Dist2_Par=NA, Con1="Rainfall",Con2="OsWL", GPD_con1, GPD_con2, mu, GPD_Bayes=FALSE, Decimal_Place=2, RP, Interval=10000, End=F, Resolution="Low", x_lab="Rainfall (mm)",y_lab="O-sWL (mNGVD 29)",x_lim_min = NA,x_lim_max = NA,y_lim_min = NA,y_lim_max = NA,Isoline_Probs="Sample", N=10^6,N_Ensemble=0,Sim_Max=10,Plot_Quantile_Isoline=FALSE,Isoline_Type="Combined",EL_con1,EL_con2){
  
  ###Preliminaries
  
  #Vectors and lists for results
  Quantile_Isoline_1<-vector(mode = "list", length = length(RP)) 
  names(Quantile_Isoline_1)<-RP
  Quantile_Isoline_2<-vector(mode = "list", length = length(RP))
  names(Quantile_Isoline_2)<-RP
  Isoline<-vector(mode = "list", length = length(RP))
  names(Isoline)<-RP
  Contour<-vector(mode = "list", length = length(RP))
  names(Contour)<-RP
  Ensemble<-vector(mode = "list", length = length(RP))
  names(Ensemble)<-RP
  MostLikelyEvent<-vector(mode = "list", length = length(RP))
  names(MostLikelyEvent)<-RP
  FullDependence<-vector(mode = "list", length = length(RP))
  names(FullDependence)<-RP
  
  #Remove 1st column of Data if it is a Date or factor object.
  if(class(Data[,1])=="Date" | class(Data[,1])=="factor" | class(Data[,1])[1]=="POSIXct"){
    Data<-Data[,-1]
  }
  
  #Find the columns in Data (which should be consistent in terms of column order of the other data input objects) of conditioning variable 1 (Con1) and conditioning variable 2 (Con2).
  con1<-which(names(Data)==Con1)
  con2<-which(names(Data)==Con2)
  
  ###Fit the 4 marginal distributions (2 GPD and 2 parametric non-extreme value distributions).
  
  #Fit the GPD to the conditioned variable con1 in Data_Con1.
#  if(is.na(GPD1)==T & is.na(Thres1)==T){
#    Thres1<-quantile(na.omit(Data[,con1]),u1)
#  }
#  
#  if(is.na(GPD1)==T & GPD_Bayes==T){
#    GPD_con1<-evm(Data_Con1[,con1], th = Thres1,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
#  }
#  if(is.na(GPD1)==T & GPD_Bayes==F){
#    GPD_con1<-evm(Data_Con1[,con1], th = Thres1)
#  }
  
  #Fit the specified marginal distribution (Marginal_Dist1) to the non-conditioned variable con2 in Data_Con1.
  if(Marginal_Dist1 == "BS"){
    if(is.na(Marginal_Dist1_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con1_M[,con2])
      marginal_non_con1<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Exp"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2],"exponential")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gam(2)"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2], "gamma")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gam(3)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1_M[,con2])
      marginal_non_con1 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "GamMix(2)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1_M[,con2])
      marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "GamMix(3)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1_M[,con2])
      marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gaus"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2],"normal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "InvG"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdist(Data_Con1_M[,con2], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Logistic"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2], "logistic")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "LogN"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2],"lognormal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "TNorm"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2],"normal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Twe"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-tweedie.profile(Data_Con1_M[,con2] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Weib"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1_M[,con2], "weibull")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  
  #Fit the GPD to the conditioned variable con2 in Data_Con2.
#  if(is.na(GPD2)==T & is.na(Thres2)==T){
#    Thres2<-quantile(na.omit(Data[,con2]),u2)
#  }
#  
#  if(is.na(GPD1)==T & GPD_Bayes==T){
#    GPD_con2<-evm(Data_Con2[,con2], th=Thres2 ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
#  }
#  if(is.na(GPD1)==T & GPD_Bayes==F){
#    GPD_con2<-evm(Data_Con2[,con2], th= Thres2)
#  }
  
  ##Fit the specified marginal distribution (Marginal_Dist2) to the non-conditioned variable con1 in Data_Con2.
  if(Marginal_Dist2 == "BS"){
    if(is.na(Marginal_Dist2_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con2_M[,con1])
      marginal_non_con2<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Exp"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1],"exponential")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gam(2)"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1], "gamma")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gam(3)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2_M[,con1])
      marginal_non_con2 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "GamMix(2)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2_M[,con1])
      marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "GamMix(3)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2_M[,con1])
      marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gaus"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1],"normal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "InvG"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdist(Data_Con2_M[,con1], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Logistic"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1],"logistic")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "LogN"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1],"lognormal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "TNorm"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1],"normal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Twe"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-tweedie.profile(Data_Con2_M[,con1] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Weib"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2_M[,con1], "weibull")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    } 
  }
  
 ###################################################################################################### 
  
  
# The following part should be changed according to the selected marginal distributions (Here it is for Tweedie and Logostic dis.)   
  
  
  # Creating the parrametric space grid
  NTR<- seq(-1,8,0.1)
  RF<- seq(0,300,2)
  Pgrid<-expand.grid(NTR,RF)
  
  # find the cdf for each value
  cdfx1 = pgpd( Pgrid[,1] , loc = Thres1 , scale = exp(GPD_con1$coefficients[1]) , shape = GPD_con1$coefficients[2])
  
  marginal_non_cond_1<-tweedie.profile(Con_NTR[,2] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
  
  cdf_non_con1 = ptweedie(Pgrid[,2], xi=NULL, mean(Con_NTR[,2]), marginal_non_cond_1$phi.max, power=marginal_non_cond_1$p.max)

  
  cdfx2 = pgpd( Pgrid[,2] , loc = Thres2 , scale = exp(GPD_con2$coefficients[1]) , shape = GPD_con2$coefficients[2])
  
  marginal_non_cond_2<-fitdistr(Con_RF[,1],"logistic")
  cdf_non_con2 = plogis(Pgrid[,1], location = marginal_non_cond_2$estimate[1], scale = marginal_non_cond_2$estimate[2], lower.tail = TRUE, log.p = FALSE)
  
  
  
  
  
  
  ###########################################################################################################
  ###Generate samples from the copula models to which a kernel density estimate will be applied to estimate relative probabilities along the isoline.
  
  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj1<-BiCopSelect(pobs(Data_Con1[,1]), pobs(Data_Con1[,2]), familyset=Copula_Family1, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula associated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample<-BiCopSim(round(N*nrow(Data_Con1)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj1)
  #Transform the realizations of the conditioned variable con1 to the original scale using inverse cumulative distribution a.k.a. quantile functions (inverse probability integral transform) of the GPD contained in the u2gpd function.
  #if(is.na(GPD1)==T){
    cop.sample1.con<-u2gpd(sample[,con1], p = 1, th=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2])
  #} else{
  #  cop.sample1.con<-u2gpd(sample[,con1], p = 1, th = GPD1$Threshold, sigma = GPD1$sigma, xi= GPD1$xi)
  #}
  
  #Transform the realizations of the non-conditioned variable con2 to the original scale using the quantile function of the selected parametric (non-extreme value) distribution (Marginal_Dist1).
  if(Marginal_Dist1=="BS"){
    cop.sample1.non.con<-qbisa(sample[,con2], as.numeric(Coef(marginal_non_con1)[1]), as.numeric(Coef(marginal_non_con1)[2]))
  }
  if(Marginal_Dist1=="Exp"){
    cop.sample1.non.con<-qexp(sample[,con2], rate = as.numeric(marginal_non_con1$estimate[1]))
  }
  if(Marginal_Dist1=="Gam(2)"){
    cop.sample1.non.con<-qgamma(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), rate = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Gam(3)"){
    cop.sample1.non.con<-qGG(sample[,con2], mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients)
  }
  if(Marginal_Dist1=="GamMix(2)"){
    xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample1.non.con <- approx(cdf.MX, xx, sample[,con2])$y
  }
  if(Marginal_Dist1=="GamMix(3)"){
    xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- round(marginal_non_con1$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample1.non.con <- approx(cdf.MX, xx, sample[,con2])$y
  }
  if(Marginal_Dist1=="Gaus"){
    cop.sample1.non.con<-qnorm(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="InvG"){
    cop.sample1.non.con<-qinvgauss(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), shape = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Logistic"){
    cop.sample1.non.con<-qlogis(sample[,con2], location = as.numeric(marginal_non_con1$estimate[1]), scale = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="LogN"){
    cop.sample1.non.con<-qlnorm(sample[,con2], meanlog = as.numeric(marginal_non_con1$estimate[1]), sdlog = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="TNorm"){
    cop.sample1.non.con<-qtruncnorm(sample[,con2], a=min(Data_Con1[,con2]), mean = as.numeric(marginal_non_con1$estimate[1]), sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Twe"){
    cop.sample1.non.con<-qtweedie(sample[,con2], power=marginal_non_con1$p.max, mu=mean(Data_Con1[,con2]), phi=marginal_non_con1$phi.max)
  }
  if(Marginal_Dist1=="Weib"){
    cop.sample1.non.con<-qweibull(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), scale=as.numeric(marginal_non_con1$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame
  cop.sample1<-data.frame(cop.sample1.con,cop.sample1.non.con)
  colnames(cop.sample1)<-c("Var1","Var2")
  
  #Fit the specified copula family (Copula_Family2) to the observations in Data_Con2.
  obj2<-BiCopSelect(pobs(Data_Con2[,1]), pobs(Data_Con2[,2]), familyset=Copula_Family2, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula assoicated with Data_Con2 is proportional to the size of Data_Con2 relative to Data_Con1.
  sample<-BiCopSim(round(N*nrow(Data_Con2)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj2)
  
  #Transform the realizations of the conditioned variable con2 to the original scale using the inverse CDF (quantile function) of the GPD contained in the u2gpd function.
 # if(is.na(GPD2)==T){
    cop.sample2.con<-u2gpd(sample[,con2], p = 1, th=Thres2, sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2])
 # } else{
 #   cop.sample2.con<-u2gpd(sample[,con2], p = 1, th = GPD2$Threshold, sigma = GPD2$sigma, xi= GPD2$xi)
#  }
  
  #Transform the realizations of the non-conditioned variable con1 to the original scale using the inverse CDF (quantile function) of the selected parametric (non-extreme value) distribution (Marginal_Dist2).
  if(Marginal_Dist2=="BS"){
    cop.sample2.non.con<-qbisa(sample[,con1], as.numeric(Coef(marginal_non_con2)[1]), as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist2=="Exp"){
    cop.sample2.non.con<-qexp(sample[,con1], rate = as.numeric(marginal_non_con2$estimate[1]))
  }
  if(Marginal_Dist2=="Gam(2)"){
    cop.sample2.non.con<-qgamma(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), rate=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Gam(3)"){
    cop.sample2.non.con<-qGG(sample[,con1], mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients)
  }
  if(Marginal_Dist2=="GamMix(2)"){
    xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample2.non.con <- approx(cdf.MX, xx, sample[,con1])$y
  }
  if(Marginal_Dist2=="GamMix(3)"){
    xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- round(marginal_non_con2$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample2.non.con <- approx(cdf.MX, xx, sample[,con1])$y
  }
  if(Marginal_Dist2=="Gaus"){
    cop.sample2.non.con<-qnorm(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), sd=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="InvG"){
    cop.sample2.non.con<-qinvgauss(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), shape=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="LogN"){
    cop.sample2.non.con<-qlnorm(sample[,con1], meanlog = as.numeric(marginal_non_con2$estimate[1]), sdlog = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Logistic"){
    cop.sample2.non.con<-qlogis(sample[,con1], location = as.numeric(marginal_non_con2$estimate[1]), scale=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="TNorm"){
    cop.sample2.non.con<-qtruncnorm(sample[,con1], a=min(Data_Con2[,con1]), mean = as.numeric(marginal_non_con2$estimate[1]), sd = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Twe"){
    cop.sample2.non.con<-qtweedie(sample[,con1], power=marginal_non_con2$p.max, mu=mean(Data_Con2[,con1]), phi=marginal_non_con2$phi.max)
  }
  if(Marginal_Dist2=="Weib"){
    cop.sample2.non.con<-qweibull(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), scale=as.numeric(marginal_non_con2$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame.
  cop.sample2<-data.frame(cop.sample2.non.con,cop.sample2.con)
  colnames(cop.sample2)<-c("Var1","Var2")
  
  #Combine the data frames containg the samples from two copulas (on the original scale)
  cop.sample<-rbind(cop.sample1,cop.sample2)
  

  
  
  
  ######################################################################################################################
  # Condiioned on var 1
  sample<-BiCopSim(round(N*nrow(Data_Con1)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj1)
  cop.sample1.con<-u2gpd(sample[,1], p = 1, th=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2])
  cop.sample1.non.con<-qtweedie(sample[,2], power=marginal_non_cond_1$p.max, mu=mean(Data_Con1[,2]), phi=marginal_non_cond_1$phi.max)
  cop.sample1<-data.frame(cop.sample1.con,cop.sample1.non.con)
  colnames(cop.sample1)<-c("Var1","Var2")
  
  # Condiioned on var 2
  sample<-BiCopSim(round(N*nrow(Data_Con2)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj2)
  cop.sample2.con<-u2gpd(sample[,2], p = 1, th=Thres2, sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2])
  cop.sample2.non.con<-qlogis(sample[,1], location = as.numeric(marginal_non_cond_2$estimate[1]), scale=as.numeric(marginal_non_cond_2$estimate[2]))
  cop.sample2<-data.frame(cop.sample2.non.con,cop.sample2.con)
  colnames(cop.sample2)<-c("Var1","Var2")
  
  # Saving the copula simulations seperately
  write.csv(cop.sample1, "ETC_Cop_Sample_Con_NTR.csv", row.names=FALSE)
  write.csv(cop.sample2, "ETC_Cop_Sample_Con_RF.csv", row.names=FALSE)
  
  
  #Combine the data frames containg the samples from two copulas (on the original scale)
  cop.sample<-rbind(cop.sample1,cop.sample2)
  write.csv(cop.sample, "ETC_Cop_Sample.csv", row.names=FALSE)
  
  
  
  ##################################################################################################################
  
  
  UU1<-BiCopCDF(cdfx1, cdf_non_con1, obj1)
  UU2<-BiCopCDF(cdf_non_con2, cdfx2, obj2)
  
  #Calculate the time period spanned by the original dataset in terms of mu (only including occasions where both variables are observed).
  time.period<-nrow(Data[which(is.na(Data[,1])==FALSE & is.na(Data[,2])==FALSE),])/mu
  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con1.
  rate<-nrow(Data_Con1)/time.period
  #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con1.
  EL<-1/rate
  
  #ff1<-function(x,y){EL/(1-x-y+UU1[which(cdfx1==x & cdf_non_con1==y)]) }
  ff1<-function(x,y){EL_con1/(1-x-y+UU1[which(cdfx1==x & cdf_non_con1==y)]) }
  #Evaluate the return period at each point on the grid 'u' (the 'outer' function creates the grid internally using the points on the boundary i.e. the x and y we defined earlier).
  RP_con1 <- ff1(cdfx1,cdf_non_con1)
  
  #ff2<-function(x,y){EL/(1-x-y+UU2[which(cdf_non_con2==x & cdfx2==y)]) }
  ff2<-function(x,y){EL_con2/(1-x-y+UU2[which(cdf_non_con2==x & cdfx2==y)]) }
  #Evaluate the return period at each point on the grid 'u' (the 'outer' function creates the grid internally using the points on the boundary i.e. the x and y we defined earlier).
  RP_con2 <- ff2(cdf_non_con2,cdfx2)
  
  RP_ETC <- Pgrid
  RP_ETC[,3] <- RP_con1
  RP_ETC[,4] <- RP_con2
  
  write.csv(RP_ETC, "RP_ETC.csv", row.names=FALSE)
  #########################################################################################
  
  
  

  
  #Create a list of outputs.
  #res<-list("FullDependence" = FullDependence, "MostLikelyEvent" = MostLikelyEvent, "Ensemble"=Ensemble, "Isoline" = Isoline, "Contour"= Contour, "Quantile_Isoline_1" = Quantile_Isoline_1, "Quantile_Isoline_2" = Quantile_Isoline_2, "Threshold_1" = Thres1, "Threshold_2"=Thres2, "Cop_sample"= cop.sample)
  #return(res)
}
