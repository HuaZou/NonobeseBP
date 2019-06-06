commonTax <- function(x, 
                      occ=0.8,
                      abun=1e-6){
  prf <- x %>% rownames_to_column("tmp") %>%
    filter(tmp!="unclassed") %>% 
    filter(apply(select(., -one_of("tmp")), 1,
                 function(x){sum(x[!is.na(x)] > 0) / 
                     length(x[!is.na(x)])}) > occ) %>%
    data.frame(.) %>%
    column_to_rownames("tmp") %>%
    t() %>% data.frame()
  prf.cln <- apply(prf, 2, median) %>%
    data.frame() %>% setNames("Abundance") %>%
    rownames_to_column("tmp") %>%
    filter(Abundance > abun & nchar(tmp) < 60)
  res <- as.character(prf.cln$tmp)
  return(res)
}

pfun <- function(y,
                 x=phen){
  
  sid <- intersect(as.character(x$SampleID), colnames(y))
  phe <- x %>% filter(SampleID%in%sid)
  tax <- commonTax(y)
  prf <-  y %>% select(as.character(phe$SampleID)) %>%
    rownames_to_column("tmp") %>%
    filter(tmp%in%tax) %>%
    column_to_rownames("tmp") 
  
  fr <- levels(x$Stage)
  wilcoxFun <- function(phs, pro, cfg){
    datprf <- pro %>% select(as.character(phs$SampleID))
    res <- apply(datprf, 1, function(x, grp){
      dat <- as.numeric(x)
      p <- wilcox.test(dat ~ grp, paired = T)$p.value
      md <- median(dat)
      mdn <- tapply(dat, grp, median)
      if ( mdn[1] > mdn[2]) {
        enrich <- levels(grp)[1]
      } else {
        enrich <- levels(grp)[2]
      }
      occ <- tapply(dat, grp, function(x){
        round(sum(x > 0)/length(x), 4)})
      res <- c(p,enrich,occ,md,mdn)
      return(res)
    }, phs$Stage) %>%
      t(.) %>% data.frame(.) %>%
      rownames_to_column("Type") %>%
      varhandle::unfactor(.)
    
    colnames(res)[2:8] <- c("Pvalue", "Enrichment",
                            paste0(fr, "_occurence"), "median", 
                            paste0(fr, "_median"))
    res$Block <- cfg 
    res$Num <- nrow(phs)/2
    res.cln <- res %>% select(c(1, 9:10, 2:8)) %>%
      mutate(Pvalue=as.numeric(Pvalue),
             median=as.numeric(median)) %>%
      arrange(Pvalue, median) %>%
      mutate(FDR = p.adjust(Pvalue, method = "BH"))
    return(res.cln)
  }
  
  pheb <- phe %>% filter(Group%in%grp[1])
  res1 <- wilcoxFun(pheb, prf, grp[1]) %>% 
    arrange(Type, Pvalue) 
  
  phep <- phe %>% filter(Group%in%grp[2])
  res2 <- wilcoxFun(phep, prf, grp[2]) %>% 
    arrange(Type, Pvalue) 
  
  res <- rbind(res1, res2)
  return(res)
}

tfun <- function(y,
                 x=phen){
  
  sid <- intersect(as.character(x$SampleID), colnames(y))
  phe <- x %>% filter(SampleID%in%sid)
  tax <- commonTax(y)
  prf <-  y %>% select(as.character(phe$SampleID)) %>%
    rownames_to_column("tmp") %>%
    filter(tmp%in%tax) %>%
    column_to_rownames("tmp") 
  
  wilcoxFun <- function(tag1, tag2, phs=phe, pro=prf){
    pht <- phs %>% filter(Stage%in%tag1 & Group%in%tag2) %>% 
      mutate(Group=factor(Group, levels=tag2))
    datprf <- pro %>% select(as.character(pht$SampleID))
    pr <- levels(pht$Group)
    res <- apply(datprf, 1, function(x, grp){
      dat <- as.numeric(x)
      p <- wilcox.test(dat ~ grp, paired = F)$p.value	
      md <- median(dat)
      mdn <- tapply(dat, grp, median)
      if ( mdn[1] > mdn[2]) {
        enrich <- levels(grp)[1]
      } else {
        enrich <- levels(grp)[2]
      }
      occ <- tapply(dat, grp, function(x){
        round(sum(x > 0)/length(x), 4)})
      
      res <- c(p,enrich,occ,md,mdn)
      return(res)
    }, pht$Group) %>% 
      t(.) %>% data.frame(.) %>%
      rownames_to_column("Type") %>%
      varhandle::unfactor(.)
    
    colnames(res)[2:8] <- c("Pvalue", "Enrichment",
                            paste0(pr, "_occurence"), "median", 
                            paste0(pr, "_median"))  
    res$Block <- paste0(tag1, " ", pr[1], "_vs_", pr[2])
    number <- as.numeric(table(pht$Group))
    res$Num <- paste0(pr[1], number[1], "_vs_",
                      pr[2], number[2])
    
    res.cln <- res %>% select(c(1,9:10, 2:8)) %>%
      mutate(Pvalue=as.numeric(Pvalue),
             median=as.numeric(median)) %>%      
      arrange(median, Pvalue) %>%
      mutate(FDR=p.adjust(Pvalue, method = "BH"))
    
    return(res.cln)
  }
  res1 <- wilcoxFun("Before", grp)
  res2 <- wilcoxFun("After", grp)
  res <- rbind(res1, res2)
  return(res)
} 
