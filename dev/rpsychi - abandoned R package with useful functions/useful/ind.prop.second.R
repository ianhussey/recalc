# calcualte a Z test for the equality of two proportions using x (counts of successes) and n (sample values)

ind.prop.second <-
function(x, n, sig.level=.05, digits=3, ref.ind=1){
     x <- as.numeric(x); n <- as.numeric(n)     
     
     n1 <- n[1]; n2 <- n[2]
     p1 <- x[1]/n1; p2 <- x[2]/n2
   
     
     ##sample statistics
     phi1 <- 2 * asin(sqrt(p1))
     phi2 <- 2 * asin(sqrt(p2))
     cohen.h <- abs(phi1 - phi2)

     samp.stat <- round(c(p1=p1, n1=n1, p2=p2, n2=n2, cohen.h=cohen.h), digits)

     
     ##RD
     if(ref.ind==1){
        rd <- p2 - p1
     }else if(ref.ind==2){
        rd <- p1 - p2
     }
     rd.sem <- sqrt((p1*(1-p1))/n1 + (p2*(1-p2))/n2)
     rd.lower <- rd + qnorm(sig.level/2, lower.tail=TRUE) * rd.sem     
     rd.upper <- rd + qnorm(sig.level/2, lower.tail=FALSE) * rd.sem
     risk.difference <- round(c(es=rd, lower=rd.lower, upper=rd.upper, std=rd.sem), digits)
     
     
     ##RR    
     if(ref.ind==1){
        rr <- p2/p1
     }else if(ref.ind==2){
        rr <- p1/p2
     }     
     log.rr.sem <- sqrt((1-p1)/(n1*p1) + (1-p2)/(n2*p2))
     rr.lower <- exp(log(rr) + qnorm(sig.level/2, lower.tail=TRUE) * log.rr.sem)
     rr.upper <- exp(log(rr) + qnorm(sig.level/2, lower.tail=FALSE) * log.rr.sem)
     risk.ratio <- round(c(es=rr, lower=rr.lower, upper=rr.upper, log.rr.sem=log.rr.sem), digits)
     

     ##OR
     if(ref.ind==1){
        or <- (p2/(1-p2))/(p1/(1-p1))
     }else if(ref.ind==2){
        or <- (p1/(1-p1))/(p2/(1-p2))
     }     
     log.or.sem <- sqrt(1/(n1*p1 *(1-p1)) + 1/(n2*p2*(1-p2)))
     or.lower <- exp(log(or) + qnorm(sig.level/2, lower.tail=TRUE) * log.or.sem)
     or.upper <- exp(log(or) + qnorm(sig.level/2, lower.tail=FALSE) * log.or.sem)
     odds.ratio <- round(c(es=or, lower=or.lower, upper=or.upper, log.or.sem=log.or.sem), digits)
     
    
     #power
     criterion.power <- round(c(small=power.prop(h=.2, n=n, sig.level=sig.level), medium=power.prop(h=.5, n=n, sig.level=sig.level), large=power.prop(h=.8, n=n, sig.level=sig.level)), digits)
      
     output <- list(samp.stat=samp.stat, risk.difference=risk.difference, risk.ratio=risk.ratio, odds.ratio=odds.ratio, power=criterion.power)
     
     return(output)     
}

# ##Kline (2004) Chapter 5
# x1 <- c("relapsed", "not relapsed")
# y1 <- c("control", "treatment")
# 
# dat <- data.frame(y =         
#                     factor(c(rep(x1, c(60, 40)), rep(x1, c(40, 60))), levels=x1),
#                   x = factor(rep(y1, each=100), levels=y1)
# )
# tab <- xtabs(~x+y, data=dat)
# tab
# ind.prop.second(x=tab[,1], n = rowSums(tab))             #Risk for relapse is lower in treatment than control condition.
# ind.prop.second(x=tab[,1], n = rowSums(tab), ref.ind=2)  #Risk for relapse is higher in control than treatment condition