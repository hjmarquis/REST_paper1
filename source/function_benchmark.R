#===============================================#
#                                               #
#     Analysis Ic: RCT vs RWE with LAC cum      #
#                                               #
#===============================================#
# 
# ipcw.grp = proj.grp = grp.all
# max.year=5

  ate.ehr.phecode.bench = function(surv,trt,Z.main, Z.inter, Z.main.inter, 
                                   ipcw.grp, proj.grp, 
                         surv.rct, trt.rct, Z.rct, #Z.match, 
                          max.year=5, or.pen, ps.pen, or.lambda, ps.lambda,
                          pos.exclude, 
                          dr.pen, dr.lambda 
  )
  {
    
    # RCT
    #==========================================
    m1.rct = do.call(data.frame, 
                     summary(survfit(surv.rct~1, subset = trt.rct==1))[c("time","surv")])
    m0.rct = do.call(data.frame, 
                     summary(survfit(surv.rct~1, subset = trt.rct==0))[c("time","surv")])
    
    out = list(mu1.rct = m1.rct, mu0.rct = m0.rct)
    
    
    # EHR project to RCT
    #==========================================
    
    # base.pos = !grepl("17",colnames(Z.main))|grepl
    main.pos = 1+1:ncol(Z.main)
    inter.pos = 1+ncol(Z.main)+1:ncol(Z.inter)
    
    # Z.rct.main = matrix(0,nrow(Z.rct), ncol(Z.main),
    #                      dimnames = list(1:nrow(Z.rct),colnames(Z.main)))
    # Z.rct.main[,colnames(Z.rct)] = Z.rct
    # Z.rct.inter = matrix(0,nrow(Z.rct), ncol(Z.main.inter),
    #                      dimnames = list(1:nrow(Z.rct),colnames(Z.main.inter)))
    # inter.var = colnames(Z.main.inter)[-grep("1017|1417", colnames(Z.main.inter))]
    # Z.rct.inter[, inter.var] = Z.rct[,inter.var]
    
    # Train the four models
    #-------------------------------------
    
    # OR
    if(missing(or.pen))
    {
      or.ridge = cv.glmnet(cbind(trt,Z.main,Z.inter), surv,
                           family = "cox", 
                           alpha = 0)
      
      out$or.pen = or.pen =  1/abs(coef(or.ridge, s = "lambda.min"))
    }
    
    if(missing(or.lambda))
    {
      or.ada = cv.glmnet(cbind(trt,Z.main,Z.inter), surv,
                         family = "cox",
                         alpha = 1, penalty.factor = or.pen)
      or.coef = as.numeric(coef(or.ada, s= "lambda.min"))
      out$or.lambda = or.lambda = or.ada$lambda.min
    }else
    {
      or.ada = glmnet(cbind(trt,Z.main,Z.inter), surv,
                      family = "cox", 
                      alpha = 1, penalty.factor = or.pen)
      or.coef = as.numeric(coef(or.ada, s= or.lambda))
    }
    names(or.coef) = c("LAC", colnames(Z.main), 
                       colnames(Z.inter))
    
    out$rr = pred = exp(trt * or.coef[1]+drop(Z.main %*% or.coef[main.pos]) + 
                 drop(Z.inter %*% or.coef[inter.pos]))
    out$rr1 = exp( or.coef[1]+drop(Z.main %*% or.coef[main.pos]) + 
                   drop(Z.main.inter %*% or.coef[inter.pos]))
    out$rr0 = exp(drop(Z.main %*% or.coef[main.pos]))
    hazfun = cox.breslow(pred,surv)
    
    # Censoring distribution
    n = length(surv)
    ipcw = rep(0,n)
    for(lev in levels(ipcw.grp))
    {
      lev0 = (trt==0)&(ipcw.grp==lev)
      lev1 = (trt==1)&(ipcw.grp==lev)
      G0 = survfit(Surv(surv[,1],1-surv[,2])~1, subset = lev0)
      G1 = survfit(Surv(surv[,1],1-surv[,2])~1, subset = lev1)
      ipcw[lev0] = 1/stepfun(G0$time, c(1,G0$surv))(pmin(surv[lev0,1],max.year))
      ipcw[lev1] = 1/stepfun(G1$time, c(1,G1$surv))(pmin(surv[lev1,1],max.year))
    }
    ipcw = ipcw*surv[,2]
    
    # Calculating the ATE
    #=============================================
    
    time.grid = sort(unique(surv[surv[,2]==1,1]))
    time.grid = time.grid[time.grid <= max.year]
    mu1.all = mu0.all = data.frame(time = time.grid)
    
    # Exclude samples violating positivity
    if(missing(pos.exclude))
    {
      # Initial PS
      if(missing(ps.pen))
      {
        ps.ridge = cv.glmnet(Z.main, trt, family = "binomial",
                             alpha = 0)
        out$ps.pen = ps.pen = 1/abs(coef(ps.ridge, s = "lambda.min")[-1])
      }
      if(missing(ps.lambda))
      {
        ps.ada = cv.glmnet(Z.main, trt, family = "binomial",
                           alpha = 1, penalty.factor = ps.pen)
        ps.coef = as.numeric(coef(ps.ada, s= "lambda.min"))
        out$ps.lambda = ps.lambda = ps.ada$lambda.min
      }else
      {
        ps.ada = glmnet(Z.main, trt, family = "binomial",
                        alpha = 1, penalty.factor = ps.pen)
        ps.coef = as.numeric(coef(ps.ada, s= ps.lambda))
      }
      names(ps.coef) = c("Intercept", colnames(Z.main))
      
      out$ps.old = ps = expit(drop(cbind(1,Z.main) %*% ps.coef))
      out$pos.exclude = pos.exclude = which((ps < 0.1) | (ps > 0.9))
      out$ps.exclude = ps[pos.exclude]
    }
    
    Z.full = Z.main
    if(length(pos.exclude)>0)
    {
      surv = surv[-pos.exclude]
      trt = trt[-pos.exclude]
      Z.main = Z.main[-pos.exclude,]
      Z.inter = Z.inter[-pos.exclude,]
      Z.main.inter = Z.main.inter[-pos.exclude,]
      proj.grp = proj.grp[-pos.exclude]
      ipcw = ipcw[-pos.exclude]

    }
    
    # Refit the PS
    ps.ada = glmnet(Z.main, trt,  
                    family = "binomial",
                    penalty.factor = ps.pen)
    ps.coef = drop(coef(ps.ada, s = ps.lambda))
    ps = expit(ps.coef[1] + drop(Z.main %*% ps.coef[-1]))
    out$ps.new = expit(ps.coef[1] + drop(Z.full %*% ps.coef[-1]))
    
    # Fit density ratio model
    bal.var = colnames(Z.rct)
    
    if(missing(dr.pen))
    {
      dr.ridge = cv.glmnet(rbind(Z.main[,bal.var],
                                 Z.rct), 
                           rep(0:1, c(nrow(Z.main),nrow(Z.rct))),
                           family = "binomial",
                           alpha = 0)
      out$dr.pen = dr.pen = 1/abs(coef(dr.ridge, s = "lambda.min")[-1])
    }
    if(missing(dr.lambda))
    {
      dr.ada = cv.glmnet(rbind(Z.main[,bal.var],
                               Z.rct), 
                         rep(0:1, c(nrow(Z.main),nrow(Z.rct))),
                         family = "binomial",
                         alpha = 1, penalty.factor = dr.pen)
      out$dr.coef = dr.coef = as.numeric(coef(dr.ada, s= "lambda.min"))
      out$dr.lambda = dr.lambda= dr.ada$lambda.min
    }else
    {
      dr.ada = glmnet(rbind(Z.main[,bal.var],
                            Z.rct), 
                      rep(0:1, c(nrow(Z.main),nrow(Z.rct))),
                      family = "binomial",
                      alpha = 1, penalty.factor = dr.pen)
      out$dr.coef = as.numeric(coef(dr.ada, s= dr.lambda))
    }
    names(out$dr.coef) = c("Intercept",colnames(Z.main[,bal.var]))
    dr = exp(-drop(cbind(1,Z.main[,bal.var]) %*% out$dr.coef))
    dr = dr/mean(dr)
    
    
    # Across all years
    #---------------------------------------------
    
    
    # OR prediction
    haz = hazfun(time.grid)
    
    pred1 = exp( or.coef[1]+drop(Z.main %*% or.coef[main.pos]) + 
                   drop(Z.main.inter %*% or.coef[inter.pos]))
    pred0 = exp(drop(Z.main %*% or.coef[main.pos]))
    
    mu1.all$OR = apply(exp(-outer(pred1,haz))*dr,2,mean)
    mu0.all$OR = apply(exp(-outer(pred0,haz))*dr,2,mean)
    
    # IPW prediction
    Y.dich = outer(surv[,1], time.grid, "<=")*ipcw
    iptw1 = trt/ps
    iptw1 = iptw1/mean(iptw1)
    iptw0 = (1-trt)/(1-ps)
    iptw0 = iptw0/mean(iptw0)
    mu1.all$IPW = 1-apply(Y.dich *iptw1*dr ,2,mean)
    mu0.all$IPW = 1-apply(Y.dich *iptw0*dr ,2,mean)
    
    # DR prediction
    mu1.all$DR = mu1.all$OR + mu1.all$IPW - 1+ apply((1-exp(-outer(pred1,haz)))*iptw1*dr,2,mean)
    mu0.all$DR = mu0.all$OR + mu0.all$IPW - 1 + apply((1-exp(-outer(pred0,haz)))*iptw0*dr,2,mean)
    
    # Stratified by years: 2006-2017
    #---------------------------------------------
    
    mu1.year = mu0.year = vector("list",nlevels(proj.grp))
    names(mu1.year) = names(mu0.year) = levels(proj.grp)
    
    # target population: RCT
    for(target in levels(proj.grp))
    {
      # Target year 2014-2017
      lev = target
      mu1.year[[lev]]  = mu0.year[[lev]]  = data.frame(time = time.grid)
      lev.pos = proj.grp == lev
      
      # OR prediction
      haz = hazfun(time.grid)
      
      dr.lev = dr[lev.pos]/mean(dr[lev.pos])
      
      mu1.year[[target]] $OR = apply(exp(-outer(pred1[lev.pos],haz))*dr.lev,2,mean)
      mu0.year[[target]] $OR = apply(exp(-outer(pred0[lev.pos],haz))*dr.lev,2,mean)
      
      # IPW prediction
      iptw1.lev = iptw1[lev.pos]*dr.lev
      iptw0.lev = iptw0[lev.pos]*dr.lev
      
      Y.dich = outer(surv[lev.pos,1], time.grid, "<=")*ipcw[lev.pos]
      mu1.year[[target]] $IPW = 1-apply(Y.dich *iptw1.lev ,2,mean)
      mu0.year[[target]] $IPW = 1-apply(Y.dich *iptw0.lev ,2,mean)
      
      # DR prediction
      mu1.year[[target]] $DR = (mu1.year[[target]]$OR + mu1.year[[target]]$IPW - 1+
                                      apply((1-exp(-outer(pred1[lev.pos],haz)))*iptw1.lev,2,mean))
      mu0.year[[target]] $DR = (mu0.year[[target]]$OR + mu0.year[[target]]$IPW - 1 +
                                      apply((1-exp(-outer(pred0[lev.pos],haz)))*iptw0.lev,2,mean))
    }
    
    
    
    return(c(list(mu1.all = mu1.all, mu0.all = mu0.all, 
                  mu1.year = mu1.year, mu0.year = mu0.year, 
                  or = or.coef, ps = ps.coef),out))
    
    
  }
  
  

