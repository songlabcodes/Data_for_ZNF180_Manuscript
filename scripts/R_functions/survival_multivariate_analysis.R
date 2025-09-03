run_kmplot.wt.group <- function(sdf,group.col,time.col,event.col)
{
  require(ggplot2)
  require(ggkm)
  require(survival)
  
  # check best split of gene expression groups
  # if (is.null(qtiles)) qtiles = seq(0.1,0.9,0.1)
  tcol <- match(time.col,colnames(sdf))
  ecol <- match(event.col,colnames(sdf))
  gcol = match(group.col,colnames(sdf))
  
  lr.test <- survdiff(Surv(sdf[[tcol]],sdf[[ecol]]) ~ sdf[[gcol]])
  logrank.p = 1-pchisq(lr.test$chisq,1);
  
  str <- paste("Logrank P=",format(signif(logrank.p,3),scientific = TRUE),sep = "")
  
  kmobj <- ggplot(data = sdf,
                  aes_string(time = time.col,status = event.col,colour = group.col)) + 
    geom_km() + geom_kmticks() + labs(x = "Time",y = "Survival") + 
    labs(title = str) + 
    theme_bw() 
  output = list(kmplot = kmobj,logrank.p = logrank.p)
  return(output)
}

run_kmplot <- function(sdf,gene.id,time.col,event.col,qtiles = 0.5,thresh = NULL)
{
  require(ggplot2)
  require(ggkm)
  require(survival)
  
  # check best split of gene expression groups
  # if (is.null(qtiles)) qtiles = seq(0.1,0.9,0.1)
  tcol <- match(time.col,colnames(sdf))
  ecol <- match(event.col,colnames(sdf))
  
  gxpr <- sdf[[match(gene.id,colnames(sdf))]]
  
  if (!is.null(qtiles))
  {
    
    logrank.p <- rep(NA,length(qtiles))
    gene.vec <- NA
    minp <- 1
    minq <- NA
    for (i in 1:length(qtiles))
    {
      gvec <- rep(NA,length(gxpr))
      qval <- quantile(gxpr,qtiles[i],na.rm = TRUE)
      gvec[which(gxpr > qval)] <- "high"
      gvec[which(gxpr <= qval)] <- "low"
      
      lr.test <- survdiff(Surv(sdf[[tcol]],sdf[[ecol]]) ~ gvec)
      logrank.p[i] <- 1-pchisq(lr.test$chisq,1);
      
      if (logrank.p[i] < minp) {gene.vec <- gvec;minp <- logrank.p[i];minq = qtiles[i];}
    }
    str <- paste("Logrank P=",format(signif(minp,3),scientific = TRUE),"(",minq*100,"%)",sep = "")
  }else{
    if (!is.null(thresh))
    {
      gene.vec = rep(NA,length(gxpr))
      gene.vec[gxpr > thresh] = "high"
      gene.vec[gxpr <= thresh] = "low"
      
      lr.test <- survdiff(Surv(sdf[[tcol]],sdf[[ecol]]) ~ gene.vec)
      minp <- 1-pchisq(lr.test$chisq,1);
      str <- paste("Logrank P=",format(signif(minp,3),scientific = TRUE),sep = "")
    }else{
      stop("provide threshold value to divide groups!")
    }
  }
  
  sdf$gene.group <- gene.vec
  
  kmobj <- ggplot(data = sdf,
                    aes_string(time = time.col,status = event.col,colour = "gene.group")) + 
    geom_km() + geom_kmticks() + geom_kmband() + labs(x = "Time",y = "Survival") + 
    scale_colour_manual(values = c("high" = "red","low" = "blue")) +
    labs(title = str) + 
    theme_bw() + guides(colour = guide_legend(title = paste(gene.id,"\nexpression group",sep = "")))
  
  if (!is.null(qtiles)) output <- list(kmplot = kmobj,quantile.result = data.frame(quantile = qtiles,logrank.p = logrank.p))
  if (!is.null(thresh)) output <- list(kmplot = kmobj,thresh.result = data.frame(threshold = thresh,logrank.p = minp))
  return(output)
}


run_multivar_survival <- function(cif,data.mat,time.col,event.col,do.univar = TRUE,do.multivar = FALSE,
                                  covariates = NULL,
                                  do.group = TRUE,prob.cut = 0.5)
{
  require(survival)
  
  mult.res = vector("list",nrow(data.mat));
  names(mult.res) = rownames(data.mat)
  for (i in 1:nrow(data.mat))
  {
    if ((i %% 200) == 0) cat(paste("i=",i,"\n",sep = ""))
    gene = rownames(data.mat)[i]
    sdata = cif[,c(1,match(c(time.col,event.col),colnames(cif)))];
    sdata$expression = as.vector(data.mat[i,match(cif[[1]],colnames(data.mat))])
    colnames(sdata)[colnames(sdata) == time.col] = "surv.time"
    colnames(sdata)[colnames(sdata) == event.col] = "surv.outcome"
    
    if (!is.null(covariates)) sdata = cbind.data.frame(sdata,covariates[match(sdata[[1]],covariates[[1]]),-1])
    if (do.group) 
    {
      sdata$expression.group = rep("Low",nrow(sdata))
      sdata$expression.group[sdata$expression > quantile(sdata$expression[!is.nan(sdata$expression) & !is.na(sdata$expression)],probs = prob.cut)] = "High"
    }
    
    # get Surv objct
    time.val = as.numeric(sdata[[which(colnames(sdata) == "surv.time")]])
    event.val = as.numeric(sdata[[which(colnames(sdata) == "surv.outcome")]]);
	geneval = sdata$expression
	
    nas = is.na(time.val) | is.na(event.val) | is.na(geneval) | is.nan(geneval) | is.infinite(geneval)
	
	if (sum(!nas) >= 10 & sum(event.val[!nas]) >= 5)
	{
		sdata = sdata[!nas,]
		time.val = time.val[!nas]
		event.val = event.val[!nas]
		geneval = geneval[!nas]
		
		srv.obj = Surv(time.val,event.val)
	  
		out = c()
		if (do.multivar)
		{
		  dat = sdata[,-c(1,which(colnames(sdata) %in% c("surv.time","surv.outcome")))]
		  # get age and race and gender 
		  res.cox <- coxph(srv.obj ~ .,data = dat);
		  sdat = summary(res.cox);
		  beta = sdat$coef["expression","coef"]
		  HR = sdat$coef["expression","exp(coef)"]
		  HR.lower = sdat$conf.int["expression","lower .95"]
		  HR.upper = sdat$conf.int["expression","upper .95"]
		  p.value.gene = sdat$coef["expression","Pr(>|z|)"]
		  p.value.model = sdat$waldtest["pvalue"]
		  
		  out = c("N" = nrow(sdata),"beta.mult" = beta,"HR.mult" = HR,"HR.mult.lower.95" = HR.lower,"HR.mult.upper.95" = HR.upper,"p.value.gene.mult" = p.value.gene,"p.value.model.mult" = p.value.model)
		}
		
		if (do.univar)
		{
		  #geneval = sdata$expression
		  res.ucox <- coxph(srv.obj ~ geneval);
		  sdat = summary(res.ucox);
		  beta = sdat$coef["geneval","coef"]
		  HR = sdat$coef["geneval","exp(coef)"]
		  HR.lower = sdat$conf.int["geneval","lower .95"]
		  HR.upper = sdat$conf.int["geneval","upper .95"]
		  p.value.gene = sdat$coef["geneval","Pr(>|z|)"]
		  
		  uout = c("beta.uni" = beta,"HR.uni" = HR,"HR.uni.lower.95" = HR.lower,"HR.uni.upper.95" = HR.upper,"p.value.gene.uni" = p.value.gene)
		  out = c(out,uout)
		}
		
		if (do.group)
		{
		  gene.group = factor(sdata$expression.group)
		  if (length(unique(gene.group)) == 2)
		  {
		    diff.cox = survdiff(srv.obj ~ (gene.group))
		    kmpvalue= 1-pchisq(diff.cox$chisq,1)
		    cnt.tbl = table(gene.group,event.val)
		    out = c(out,c("logrank.p" = kmpvalue,"N.High.Event" = cnt.tbl["High","1"],"N.Low.Event" = cnt.tbl["Low","1"]))
		    rm(gene.group,diff.cox,kmpvalue,cnt.tbl)
		  }else{
		    out = c(out,c("logrank.p" = NA,"N.High.Event" = NA,"N.Low.Event" = NA))
		  }
		  
		}
		mult.res[[i]] = out
	}
    
    rm(gene,sdata,time.val,event.val,srv.obj,age,gender,race,histol,geneval,sdat,beta,HR,HR.lower,HR.upper,p.value.gene,p.value.model)
  }
  mult.res = mult.res[!sapply(mult.res,is.null)]
  mult.res.comb = do.call('rbind',mult.res)
  mult.res = data.frame(gene = names(mult.res),as.data.frame(mult.res.comb))
  if (any(colnames(mult.res) == "p.value.gene.mult")) mult.res = mult.res[order(mult.res$p.value.gene.mult),]
  if (any(colnames(mult.res) == "p.value.gene.mult")) mult.res$FDR.gene.mult = p.adjust(mult.res$p.value.gene.mult)
  if (any(colnames(mult.res) == "p.value.model.mult")) mult.res$FDR.model.mult = p.adjust(mult.res$p.value.model.mult)
  if (do.univar) mult.res$FDR.gene.uni = p.adjust(mult.res$p.value.gene.uni)
  if (do.group) mult.res$FDR.logrank = p.adjust(mult.res$logrank.p)
  mult.res$outcome = rep(event.col,nrow(mult.res))
  return(mult.res)
}
