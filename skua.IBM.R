bplist00—_WebMainResource’	
_WebResourceTextEncodingName^WebResourceURL_WebResourceFrameName_WebResourceData_WebResourceMIMETypeUUTF-8_Nhttp://datadryad.org/bitstream/handle/10255/dryad.101319/skua.IBM.R?sequence=1PO,B<html><head></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">###############  Individual based model to simulate sex specific arrival times  ###############
##########                                                                           ##########
##########          This model is based on a code written by Hannah Kokko            ##########
##########                            and described in:                              ##########
##########                 Kokko et al. (2006) Why do female migratory               ##########
##########  birds arrive later than males? Journal of Animal Ecology. 75             ##########
##########                                                                           ##########
###############################################################################################
########### Authors: Simeon Lisovski &amp; Markus Ritz ############################################
############################################################################################### 


## Parameters:
## NM: Number of South Polar Skua males
## NF: Number of South Polar Skua females
## alpha: Male mate choice parameter (1 = no preference, &gt; 1 prevalence towards Brown skua females)
## EPCprob: Extra pair paternity rate
## sigma: ratio of allelic values in hybrids (1 = all hybrids inherit allelic value from South Polar Skua father)
## mutprob: Allel mutation probability
## mortality: mortality rate
## NrT: numbers of Territories
## QT: maximal quality of each terrirory
## plot: logical, if TRUE a plot with the arrival times will be produced



skua.IBM &lt;- function(NM, NF, BS, alpha, EPCprob, sogma, mutprob, mortality, NrT, QT, plot=TRUE) {
  
  # indices
  # index=1; a=2; b=3; arrtime=4; aEPC=5; bEPC=6
  
  # Females and Males contain:
  # index; genes a and b; arrtime
  
  Fe &lt;- matrix(c(1:NF, 40*runif(NF*2)), ncol=3, nrow=NF)
  Ma &lt;- matrix(c(1:NM, 40*runif(NM*2)), ncol=3, nrow=NM)
  BS &lt;- matrix(ncol=2, nrow = BS)
  if(nrow(BS)&gt;0) BS[,1] &lt;- 1:nrow(BS)
  
  Tr &lt;- matrix(c(sort(round(runif(NrT, 1, QT)), decreasing=T), 
                 rep(NA, 2*NrT)), nrow=NrT, ncol=3)
  
  Fe &lt;- cbind(Fe, matrix(ncol=3, nrow=NF))
  Ma &lt;- cbind(Ma, matrix(ncol=1, nrow=NM))
  
  data &lt;- cbind(NA, NA, NA)
  
  generation &lt;- 0
  
  repeat{
    
    Fe[ , 4] &lt;- NA
    Ma[ , 4] &lt;- NA
    Fe[ , c(5,6)] &lt;- NA
    
    if(!is.null(BS)) BS[, 2] &lt;- NA
    
    if(nrow(Fe)&gt;NF) Fe &lt;- Fe[1:NF, ]
    if(nrow(Ma)&gt;NM) Ma &lt;- Ma[1:NM, ]
    
    Tr[ , 2] &lt;- NA
    Tr[ , 3] &lt;- NA
    
    mcurve &lt;- rep(NA, 40)
    fcurve &lt;- rep(NA, 40)
    
    for( t in 0:40) {
      
      # migration: everyone in state 0 moves up if they so wish
      move &lt;- which(is.na(Ma[ , 4]) &amp; runif(nrow(Ma)) &lt; 
                      (1/ (1 + exp( -2 * (t - Ma[ , 2])))))
      if(length(move)&gt;0) Ma[move, 4] &lt;- t
      move &lt;- which(is.na(Fe[ , 4]) &amp; runif(nrow(Fe)) &lt; 
                      (1/ (1 + exp( -2 * (t - Fe[ , 3])))))
      if(length(move)&gt;0) Fe[move, 4] &lt;- t
      
      
      # females who have arrived but not yet chosen an EPC mate do so
      epc &lt;- which(!is.na(Fe[ , 4]) &amp; is.na(Fe[ , 5]))
      if(length(epc)&gt;0) {
        f &lt;- which(!is.na(Ma[ , 4]))
        if(length(f)&gt;0) {
          f &lt;- f[sample(length(f),length(epc), replace=T)]
          Fe[epc, c(5,6)] &lt;- Ma[f, c(2,3)]
        }
      }
      
      
      mcurve[t+1] &lt;- mean(Ma[, 4], na.rm=T)
      fcurve[t+1] &lt;- mean(Fe[, 4], na.rm=T)
      
      
      # some mortality
      malemort &lt;- exp(-mortality*t)
      dead &lt;- which(!is.na(Ma[, 4]) &amp; runif(nrow(Ma))&lt;malemort)
      if(length(dead)&gt;0) {
        Ma &lt;- Ma[-dead, ]
        ind &lt;-  setdiff(Tr[, 2], Ma[, 1])
        if(length(ind)&gt;0) Tr[which(Tr[,2]%in%ind), 2] &lt;- NA
        ind &lt;- setdiff(BS[,2], Ma[, 1])
        if(length(ind)&gt;0) BS[which(BS[,2]%in%ind), 2] &lt;- NA
      }
      
      femalemort &lt;- exp(-mortality*t)
      dead &lt;- which(!is.na(Fe[, 4]) &amp; runif(nrow(Fe))&lt;femalemort)
      if(length(dead)&gt;0) {
        Fe &lt;- Fe[-dead, ]
        ind &lt;- setdiff(Tr[, 3], Fe[, 1])
        if(length(ind)&gt;0) Tr[which(Tr[,3]%in%ind), 3] &lt;- NA
      }
      
      
      # males without a territory look for single BS females
      f1    &lt;-  setdiff(Ma[!is.na(Ma[, 4]), 1], c(Tr[, 2], BS[, 2]))
      f2BS  &lt;-  which(is.na(BS[, 2]))
      f2SPS &lt;-  which(is.na(Tr[, 2]) &amp; !is.na(Tr[, 3]))
      
      if(length(f1)&gt;0 &amp; (length(f2BS)&gt;0 | length(f2SPS)&gt;0)) 
        f1 &lt;- f1[runif(length(f1))&lt;((alpha/2)^(log(nrow(BS)/(nrow(BS)+nrow(Fe)), 0.5)))]
      # if(length(f1)&gt;0 &amp; (length(f2BS)&gt;0 | length(f2SPS)&gt;0)) 
      f1 &lt;- f1[runif(length(f1))&lt;alpha]
      if(length(f1)&gt;0 &amp; length(f2BS)&gt;0){
        if(length(f1)==length(f2BS)) BS[f2BS, 2] &lt;- f1
        if(length(f1) &gt; length(f2BS)) BS[f2BS, 2] &lt;- f1[1:length(f2BS)]
        if(length(f1) &lt; length(f2BS)) BS[f2BS[1:length(f1)], 2] &lt;- f1			
      }
      
      
      # males without a territory look for single territorial females
      f1 &lt;- setdiff(Ma[!is.na(Ma[, 4]), 1], c(Tr[, 2], BS[,2]))
      if(any(is.na(f1))) f1 &lt;- f1[-which(is.na(f1))]
      f2 &lt;- which(!is.na(Tr[, 3]) &amp; is.na(Tr[, 2]))
      if(length(f1)&gt;0 &amp; length(f2)&gt;0) {
        if(length(f1)==length(f2)) Tr[f2, 2] &lt;- f1
        if(length(f1) &gt; length(f2)) Tr[f2, 2] &lt;- f1[1:length(f2)] 
        if(length(f1) &lt; length(f2)) Tr[f2[1:length(f1)], 2] &lt;- f1
      }
      
      
      # those males without a territory look for territories
      f1 &lt;- setdiff(Ma[!is.na(Ma[, 4]), 1], c(Tr[, 2], BS[,2]))
      if(any(is.na(f1))) f1 &lt;- f1[-which(is.na(f1))]
      f2 &lt;- which(is.na(Tr[, 2]))
      if(length(f1)&gt;0 &amp; length(f2)&gt;0) {
        if(length(f1)==length(f2))  Tr[f2, 2] &lt;- f1
        if(length(f1) &gt; length(f2)) Tr[f2, 2] &lt;- f1[1:length(f2)]
        if(length(f1) &lt; length(f2)) Tr[f2[1:length(f1)], 2] &lt;- f1
      }	
      
      # females without a territory look for single territorial males
      f1 &lt;- setdiff(Fe[!is.na(Fe[, 4]), 1], Tr[, 3])
      if(any(is.na(f1))) f1 &lt;- f1[-which(is.na(f1))]
      f2 &lt;- which(!is.na(Tr[, 2]) &amp; is.na(Tr[, 3]))
      if(length(f1)&gt;0 &amp; length(f2)&gt;0) {
        if(length(f1)==length(f2))  Tr[f2, 3] &lt;- f1
        if(length(f1) &gt; length(f2)) Tr[f2, 3] &lt;- f1[1:length(f2)]
        if(length(f1) &lt; length(f2)) Tr[f2[1:length(f1)], 3] &lt;- f1
      }
      
      
      # those females without a territory look for territories
      f1 &lt;- setdiff(Fe[!is.na(Fe[, 4]), 1], Tr[, 3])
      if(any(is.na(f1))) f1 &lt;- f1[-which(is.na(f1))]
      f2 &lt;- which(is.na(Tr[, 3]))
      if(length(f1)&gt;0 &amp; length(f2)&gt;0) {
        if(length(f1)==length(f2))  Tr[f2, 3] &lt;- f1
        if(length(f1) &gt; length(f2)) Tr[f2, 3] &lt;- f1[1:length(f2)]
        if(length(f1) &lt; length(f2)) Tr[f2[1:length(f1)], 3] &lt;- f1
      }
      
    } # end of migration
    
    data &lt;- rbind(data, cbind(mean(Ma[, 2]), mean(Fe[, 3]), sum(!is.na(BS[,2]))))
    
    
    if(plot==TRUE){
      if(generation/2 == floor(generation/2)) {
        plot(0:40, mcurve, ylim=c(0, max(c(mcurve, fcurve), na.rm=T)), 
             type="l", col="blue", xlab="days", ylab="mean arrival")
        lines(0:40, fcurve, col="red")
      }
    }	
    
    newmale &lt;- matrix(ncol=4)
    newfemale &lt;- matrix(ncol=4)
    
    for(i in 1:nrow(Tr)) {
      
      if(!is.na(Tr[i, 2]) &amp; !is.na(Tr[i, 3]) &amp; Tr[i, 1]&gt;0) {
        f1 &lt;- which(Ma[, 1]==Tr[i, 2])
        f2 &lt;- which(Fe[, 1]==Tr[i, 3])
        
        # with certain probability some of male offspring are extra-pair...
        newest &lt;- matrix(rep(Ma[f1,], each=Tr[i, 1]), nrow=Tr[i, 1])
        
        # with 0.5 probability male genes are inherited from mum
        allelfrommom &lt;- matrix(runif(nrow(newest)*2)&lt;0.5, ncol=2)
        newest[allelfrommom[, 1], 2] &lt;- Fe[f2, 2]
        newest[allelfrommom[, 2], 3] &lt;- Fe[f2, 3]
        
        EPC &lt;- which(runif(nrow(newest))&lt;EPCprob)
        if(length(EPC)&gt;0) 
          newest[EPC, c(2,3)] &lt;- matrix(rep(Fe[f2, c(5, 6)], each=length(EPC)), ncol=2)
        
        newmale &lt;- rbind(newmale, newest)
        
        # with certain probability some of female offspring are extra-pair...
        newest &lt;- matrix(rep(Ma[f1,], each=Tr[i, 1]), nrow=Tr[i, 1])
        # EPC &lt;- which(runif(nrow(newest))&lt;EPCprob)
        # if(length(EPC)&gt;0) 
        newest[EPC, c(2,3)] &lt;- matrix(rep(Fe[f2, c(5, 6)], each=length(EPC)), ncol=2)
        
        # with 0.5 probability male genes are inherited from mum
        allelfrommom &lt;- matrix(runif(nrow(newest)*2)&lt;0.5, ncol=2)
        newest[allelfrommom[, 1], 2] &lt;- Fe[f2, 2]
        newest[allelfrommom[, 2], 3] &lt;- Fe[f2, 3]
        
        EPC &lt;- which(runif(nrow(newest))&lt;EPCprob)
        if(length(EPC)&gt;0) 
          newest[EPC, c(2,3)] &lt;- matrix(rep(Fe[f2, c(5, 6)], each=length(EPC)), ncol=2)
        
        newfemale &lt;- rbind(newfemale, newest)
      }
    }	
    
    # hybrids
    if(nrow(BS)&gt;0) {
      for(i in 1:nrow(BS)) {
        
        if(!is.na(BS[i, 2])) {
          
          ratio &lt;- round(QT*BSxSPS_ratio, 0)
          
          if(ratio&gt;0) {
            f1 &lt;- which(Ma[, 1]==BS[i, 2])
            
            newest &lt;- matrix(rep(Ma[f1, ], each=ratio), nrow=ratio)
            newmale &lt;- rbind(newmale, newest)
            
            newest &lt;- matrix(rep(Ma[f1,], each=ratio), nrow=ratio)
            newfemale &lt;- rbind(newfemale, newest)
          }
        }
      }
      
      
      newmale &lt;- newmale[-1, ]
      newfemale &lt;- cbind(newfemale[-1, ], matrix(ncol=2, nrow=nrow(newfemale)-1))
      
      if(nrow(newmale)&gt;0){
        newmale[ , 1] &lt;- seq(max(Ma[,1])+1, max(Ma[,1])+nrow(newmale))
        mut &lt;- matrix(runif(nrow(newmale)*2)&lt;mutprob, ncol=2)
        newmale[ , 2] &lt;- ifelse(mut[,1], runif(sum(mut[,1]))*40, newmale[ , 2])
        newmale[ , 3] &lt;- ifelse(mut[,2], runif(sum(mut[,2]))*40, newmale[ , 3])
      }
      
      if(nrow(newfemale)&gt;0){
        newfemale[ , 1] &lt;- seq(max(Fe[,1])+1, max(Fe[,1])+nrow(newfemale))
        mut &lt;- matrix(runif(nrow(newfemale)*2)&lt;mutprob, ncol=2)
        newfemale[ , 2] &lt;- ifelse(mut[,1], runif(sum(mut[,1]))*40, newfemale[ , 2])
        newfemale[ , 3] &lt;- ifelse(mut[,2], runif(sum(mut[,2]))*40, newfemale[ , 3])
      }
      
      # merge
      Fe &lt;- rbind(Fe, newfemale)
      Ma &lt;- rbind(Ma, newmale)
      
      Fe &lt;- Fe[sample(1:nrow(Fe)), ]
      Ma &lt;- Ma[sample(1:nrow(Ma)), ]
      
      
      generation &lt;- generation + 1
      
      STOP &lt;- FALSE
      if(generation&gt;=50){
        mod &lt;- lm(c(data[(nrow(data)-9):(nrow(data)),2] - 
                    data[(nrow(data)-9):(nrow(data)),1])~c(1:10))
        if(coef(mod)[2]&gt;=(-0.005) &amp; coef(mod)[2]&lt;=0.005) STOP &lt;- TRUE
      }
      if(STOP) break
      
    }
    
    return(data[,2]-data[,1])
  }
  </pre></body></html>Ztext/plain    ( F U l ~ î ö Î Ï-2                           -=