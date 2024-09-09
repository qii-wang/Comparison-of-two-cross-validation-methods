calibrationMetrics = function(data,risk_var,p_weak_cali_title=NULL,weak_cali_ind=TRUE,correlated=TRUE){

  # print(dim(data))
  data$risk_interest = data[,risk_var]
  
  if(correlated==TRUE){
    ## EO ratio:
    E = mean(data$risk_interest)
    O = mean(data$case)
    EORatio = E/O
    library(gee)
    gee_fit = gee(formula = case~1, family = binomial(link="logit"), 
                  data = data, id=as.factor(StudyID_c), corstr="unstructured")
    summary(gee_fit)$coefficients
    beta_hat = summary(gee_fit)$coefficients["(Intercept)","Estimate"]
    var_beta_hat = summary(gee_fit)$coefficients["(Intercept)","Robust S.E."]^2
    se_EORatio = sqrt(E^2/exp(beta_hat)^2 * var_beta_hat)
    EORatio_95CI = c(EORatio - qnorm(0.975)*se_EORatio, EORatio + qnorm(0.975)*se_EORatio)
    
    EORatio_res = paste0(round(EORatio,digits=4)," (",round(EORatio_95CI[1],digits=4),", ", round(EORatio_95CI[2],digits=4),")")
    
    ## Calibration intercept:
    library(gee)
    gee_fit = gee(formula = case~1+offset(I(log(risk_interest/(1-risk_interest)))), family = binomial(link="logit"), data = data, id=as.factor(StudyID_c), corstr="unstructured")
    cox_intercept = summary(gee_fit)$coefficient["(Intercept)","Estimate"] 
    cox_intercept_se = summary(gee_fit)$coefficient["(Intercept)","Robust S.E."]
    cox_intercept_95CI = c(cox_intercept-qnorm(0.975)*cox_intercept_se, cox_intercept+qnorm(0.975)*cox_intercept_se)
    
    cox_intercept_res = paste0(round(cox_intercept,digits=4)," (",round(cox_intercept_95CI[1],digits=4),", ", round(cox_intercept_95CI[2],digits=4),")")
    
    ## Calibration slope:
    library(gee)
    gee_fit = gee(formula = case~I(log(risk_interest/(1-risk_interest))), family = binomial(link="logit"), data = data, id=as.factor(StudyID_c), corstr="unstructured")
    cox_slope = summary(gee_fit)$coefficient["I(log(risk_interest/(1 - risk_interest)))","Estimate"]
    cox_slope_se = summary(gee_fit)$coefficient["I(log(risk_interest/(1 - risk_interest)))","Robust S.E."]
    cox_slope_95CI = c(cox_slope-qnorm(0.975)*cox_slope_se, cox_slope+qnorm(0.975)*cox_slope_se)
    
    cox_slope_res = paste0(round(cox_slope,digits=4)," (",round(cox_slope_95CI[1],digits=4),", ", round(cox_slope_95CI[2],digits=4),")")
    
    ## Weak calibration plot:
    if(weak_cali_ind){
      data$cali_grouping = ceiling(ecdf(data$risk_interest)(data$risk_interest)*10)
      library(dplyr)
      library(gee)
      weak_cali = data %>% group_by(cali_grouping) %>% 
        summarise(n = length(case),
                  E = mean(risk_interest),
                  O = mean(as.integer(case)),
                  beta_hat = summary(gee(formula = case~1, family = binomial(link="logit"),
                                         id=as.factor(StudyID_c),
                                         corstr="unstructured"))$coefficients["(Intercept)","Estimate"],
                  var_beta_hat = summary(gee(formula = case~1, family = binomial(link="logit"),
                                             id=as.factor(StudyID_c),
                                             corstr="unstructured"))$coefficients["(Intercept)","Robust S.E."]^2,
                  se_O = sqrt(exp(beta_hat)^2/(1+exp(beta_hat))^4 * var_beta_hat),
                  se_O_naive = sqrt(O*(1-O)/n),
                  lb = round(O-qnorm(0.975)*se_O,digits=4),
                  ub = round(O+qnorm(0.975)*se_O,digits=4))
      View(weak_cali)
      library(ggplot2)
      p_weak_cali = ggplot(data=weak_cali) +
        geom_point(aes(x=E,y=O),color="dark grey") +
        # xlim(0,0.008) + ylim(0,0.008) + 
        xlab("mean predicted risk") + ylab("observed risk") +
        geom_errorbar(aes(x=E,ymin=lb,ymax=ub),color="dark grey", size=0.5, width=0) +
        geom_abline(intercept = 0, slope = 1, color="white") +
        labs(title=p_weak_cali_title) +
        theme(text=element_text(family="sans",size=12))
      
      return(list(EORatio=EORatio_res,
                  cox_intercept=cox_intercept_res,
                  cox_slope=cox_slope_res,
                  p_weak_cali=p_weak_cali,
                  est_only=c(EORatio=EORatio,
                             cox_intercept=cox_intercept,
                             cox_slope=cox_slope)))
    }
    else{
      
      return(list(EORatio=EORatio_res,
                  cox_intercept=cox_intercept_res,
                  cox_slope=cox_slope_res,
                  est_only=c(EORatio=EORatio,
                             cox_intercept=cox_intercept,
                             cox_slope=cox_slope)))
    }
    
  }
  else{
    ## EO ratio:
    E = mean(data$risk_interest)
    O = mean(data$case)
    EORatio = E/O
    library(gee)
    library(sandwich)
    glm_fit = glm(formula = case~1, family = binomial(link="logit"), 
                  data = data)
    summary(glm_fit)$coefficients
    beta_hat = summary(glm_fit)$coefficients["(Intercept)","Estimate"]
    var_beta_hat = summary(glm_fit)$coefficients["(Intercept)","Std. Error"]^2
    se_EORatio = sqrt(E^2/exp(beta_hat)^2 * var_beta_hat)
    EORatio_95CI = c(EORatio - qnorm(0.975)*se_EORatio, EORatio + qnorm(0.975)*se_EORatio)
    EORatio_res = c("EORatio"=EORatio, "se_EORatio"=se_EORatio, "EORatio_cilb"=EORatio_95CI[1], "EORatio_ciub"=EORatio_95CI[2])
    
    ## Calibration intercept:
    glm_fit = glm(formula = case~1+offset(I(log(risk_interest/(1-risk_interest)))), family = binomial(link="logit"), data = data)
    cox_intercept = summary(glm_fit)$coefficient["(Intercept)","Estimate"] 
    cox_intercept_se = sqrt(diag(vcovHC(glm_fit)))["(Intercept)"][[1]]
    # cox_intercept_se = summary(glm_fit)$coefficient["(Intercept)","Std. Error"]
    cox_intercept_95CI = c(cox_intercept-qnorm(0.975)*cox_intercept_se, cox_intercept+qnorm(0.975)*cox_intercept_se)
    
    cox_intercept_res = c("cox_intercept"=cox_intercept, "cox_intercept_se"=cox_intercept_se, "cox_intercept_cilb"=cox_intercept_95CI[1], "cox_intercept_ciub"=cox_intercept_95CI[2])
    
    ## Calibration slope:
    glm_fit = glm(formula = case~I(log(risk_interest/(1-risk_interest))), family = binomial(link="logit"), data = data)
    cox_slope = summary(glm_fit)$coefficient["I(log(risk_interest/(1 - risk_interest)))","Estimate"]
    cox_slope_se = sqrt(diag(vcovHC(glm_fit)))["I(log(risk_interest/(1 - risk_interest)))"][[1]]
    # cox_slope_se = summary(glm_fit)$coefficient["I(log(risk_interest/(1 - risk_interest)))","Std. Error"]
    cox_slope_95CI = c(cox_slope-qnorm(0.975)*cox_slope_se, cox_slope+qnorm(0.975)*cox_slope_se)
    
    cox_slope_res = c("cox_slope"=cox_slope, "cox_slope_se"=cox_slope_se, "cox_slope_cilb"=cox_slope_95CI[1], "cox_slope_ciub"=cox_slope_95CI[2])
    
    ## Weak calibration plot:
    if(weak_cali_ind){
      data$cali_grouping = ceiling(ecdf(data$risk_interest)(data$risk_interest)*10)
      library(dplyr)
      library(gee)
      weak_cali = data %>% group_by(cali_grouping) %>% 
        summarise(n = length(case),
                  E = mean(risk_interest),
                  O = mean(as.integer(case)),
                  beta_hat = summary(glm(formula = case~1, family = binomial(link="logit")))$coefficients["(Intercept)","Estimate"],
                  var_beta_hat = summary(glm(formula = case~1, family = binomial(link="logit")))$coefficients["(Intercept)","Std. Error"]^2,
                  se_O = sqrt(exp(beta_hat)^2/(1+exp(beta_hat))^4 * var_beta_hat),
                  se_O_naive = sqrt(O*(1-O)/n),
                  lb = round(O-qnorm(0.975)*se_O,digits=4),
                  ub = round(O+qnorm(0.975)*se_O,digits=4))
      View(weak_cali)
      library(ggplot2)
      p_weak_cali = ggplot(data=weak_cali) +
        geom_point(aes(x=E,y=O),color="dark grey") +
        # xlim(0,0.008) + ylim(0,0.008) + 
        xlab("mean predicted risk") + ylab("observed risk") +
        geom_errorbar(aes(x=E,ymin=lb,ymax=ub),color="dark grey", size=0.5, width=0) +
        geom_abline(intercept = 0, slope = 1, color="white") +
        labs(title=p_weak_cali_title) +
        theme(text=element_text(family="sans",size=12))
      
      return(list(EORatio=EORatio_res,
                  cox_intercept=cox_intercept_res,
                  cox_slope=cox_slope_res,
                  p_weak_cali=p_weak_cali,
                  est_only=c(EORatio=EORatio,
                             cox_intercept=cox_intercept,
                             cox_slope=cox_slope)))
    }
    else{
      
      return(list(EORatio=EORatio_res,
                  cox_intercept=cox_intercept_res,
                  cox_slope=cox_slope_res,
                  est_only=c(EORatio=EORatio,
                             cox_intercept=cox_intercept,
                             cox_slope=cox_slope)))
    }
  }
}
