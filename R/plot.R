plotStrata<-function(assess.out,
                      figure.type=c('power','fdr','type1error','tp','fp','taxon'),
                      stratify.by=c('prevalence','abundance'),
                      is.errorbar=TRUE,is.legend=TRUE) {
  
  nsim=dim(assess.out$TP)[2]
  sum.out=summaryAssess(assess.out,assess.type='stratify')
  if(length(sum.out)!=2) stop('assess.type should be set as stratify in summaryAssess!')
  alpha.type=assess.out$alpha.type
  stratify.type=assess.out$stratify.type
  alpha.level=assess.out$alpha.level
  figure.type=match.arg(figure.type)
  
  xlabel=switch(stratify.by,
                prevalence = "Strata (prevalence)",
                abundance = "Strata (abundance)"
  )
  
  ylabel=switch(figure.type,
                power="Power",
                fdr="FDR",
                type1error="Type1 Error",
                tp="TP",
                fp="FP",
                taxon="Taxon")
  
  if(figure.type=='power'){
    TPR.m=sum.out$m[,,'TPR']
    TPR.s=sum.out$s[,,'TPR'] /nsim
    TPR.m=melt(TPR.m)
    TPR.s=melt(TPR.s)
    TPR.m$sd=TPR.s$value
    dat=TPR.m
    #.plotPowerStrata(TPR.m)
  }else if(figure.type=='fdr'){
    FDR.m=sum.out$m[,,'FDR']
    FDR.s=sum.out$s[,,'FDR'] /nsim
    FDR.m=melt(FDR.m)
    FDR.s=melt(FDR.s)
    FDR.m$sd=FDR.s$value
    FDR.m$alpha.level=alpha.level
    dat=FDR.m
    #.plotFDRStrata(FDR.m)
  }else if(figure.type=='type1error'){
    FPR.m=sum.out$m[,,'FPR']
    FPR.s=sum.out$s[,,'FPR'] /nsim
    FPR.m=melt(FPR.m)
    FPR.s=melt(FPR.s)
    FPR.m$sd=FPR.s$value
    FPR.m$alpha.level=alpha.level
    dat=FPR.m
    #.plotType1ErrorStrata(FPR.m)
  }else if(figure.type=='tp'){
    TP.m=sum.out$m[,,'TP']
    TP.s=sum.out$s[,,'TP'] /nsim
    TP.m=melt(TP.m)
    TP.s=melt(TP.s)
    TP.m$sd=TP.s$value
    TP.m$alpha.level=alpha.level
    dat=TP.m
  }else if(figure.type=='fp'){
    FP.m=sum.out$m[,,'FP']
    FP.s=sum.out$s[,,'FP'] /nsim
    FP.m=melt(FP.m)
    FP.s=melt(FP.s)
    FP.m$sd=FP.s$value
    FP.m$alpha.level=alpha.level
    dat=FP.m
  }else if(figure.type=='taxon'){
    taxon.m=sum.out$m[,,'Taxon']
    taxon.s=sum.out$s[,,'Taxon'] /nsim
    taxon.m=melt(taxon.m)
    taxon.s=melt(taxon.s)
    taxon.m$sd=taxon.s$value
    taxon.m$alpha.level=alpha.level
    dat=taxon.m
  }
  
  p=ggplot(dat, aes(x=stratify.name, y=value, group=samplesize, color=samplesize)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(0.05))
  p=p+labs(title="", x=xlabel, y = ylabel)+theme_bw()
  p=p+theme(axis.text.x=element_text(size=12,angle = 30, hjust = 1,face='bold'),
            axis.text.y=element_text(size=12,angle = 0, hjust = 1,face='bold'),
            axis.title=element_text(size=12,face="bold"),
            strip.text = element_text(size = 14,face="bold"),
            legend.text=element_text(size=12,face='bold'),
            legend.title=element_text(size=14,face='bold'))
  if(!is.legend) p=p+theme(legend.position='none') 
  p
  
}



plotStrataAll <- function(assess.out, 
                          figure.type=c('power','fdr','type1error','tp','fp','taxon'),
                          stratify.by=c('prevalence','abundance'),
                          is.errorbar=TRUE,is.legend=TRUE){
  p1=plotStrata(assess.out,'power',stratify.by=stratify.by,is.legend = is.legend)
  p2=plotStrata(assess.out,'fdr',stratify.by=stratify.by,is.legend = is.legend)
  p3=plotStrata(assess.out,'tp',stratify.by=stratify.by,is.legend = is.legend)
  p4=plotStrata(assess.out,'fp',stratify.by=stratify.by,is.legend = is.legend)
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
}


