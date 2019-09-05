#' \code{Evenness}: compute six evenness indices of order q.
#' @param data a data.frame, matrix or list of observed species-by-assemblage frequency vector. For the case
#' of only one assemblage (site), dat should be a numeric vector of species frequency.
#' @param q a vector of diversity orders: user must specify a vector (default  is from 0 to 2 in an increment
#' of 0.1)
#' @return a list of N dataframes, where N is the number of assemblages (sites). Each data.frame shows
#' the profiles of all six classes of evenness indices listed in Chao and Ricotta (2019) Ecology paper.
#' @examples
#' data(Alpine)
#' out1 <- Evenness(data = Alpine)
#' @export
Evenness <- function(data, q = seq(0, 2, 0.1)){
  if(class(data)== "data.frame"||class(data)== "matrix"){
    mydata <- lapply(1:ncol(data), function(i) data[,i])
    names(mydata) <- colnames(data)
  }else if (class(data)=="numeric"){
    mydata <- list(data)
    names(mydata) <- "data"
  }
  output <- lapply(mydata, function(x) cbind(q = q,new_fun(x = x, q.order = q)))
return(output)
}

#' \code{gg_evenness}: use ggplot2 and ggpubr to plot the outcome of function \code{Evenness}.
#' @param x a optput of \code{Evenness}
#' @return return a ggplot2 object
#' @examples
#' data(Alpine)
#' out1 <- Evenness(data = Alpine)
#' gg_evenness(out1)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @export
gg_evenness <- function(x){
  qs <- x[[1]][,1]
  x2 <- sapply(x, function(y) y[,-1], simplify = "array")
  x3 <- aperm(x2, c(1, 3, 2))
  myFUN <- function(x_1,evetype){
    x_1 <- as.data.frame(x_1)
    x_1$q <- qs
    x_1 <- melt(x_1, id.vars = c("q"))
    names(x_1) <- c("q", "Site", "evenness")
    out <- ggplot(x_1, aes("q", "evenness"))+
      geom_line(aes(color = "Site", group = "Site", linetype = "Site"), size = 1.1)+
      scale_linetype_manual(values=c("dashed", "1111", "solid"))+
      theme_bw()+
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text = element_text(size = 14),
            legend.position = "bottom")+
      scale_x_discrete(breaks=seq(0, 2, 0.5))+
      theme(legend.key.width = unit(2,"cm"))+
      theme(plot.title = element_text(size=20, face="bold.italic",hjust = 0.5))+
      ggtitle(evetype)+
      xlab("Diversity order q")
    return(out)
  }
  ggs <- lapply(1:6, function(i){ myFUN(x3[,,i],paste0("E",i)) })
  ggarrange(ggs[[1]],ggs[[2]],ggs[[3]],ggs[[4]],ggs[[5]],ggs[[6]],
            ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
}

new_fun <- function(x,q.order){
  n <- sum(x)
  p <- x/n
  q_profile_evenness <- function(q){
    qDest <- qD(p,q)
    #S <- sum(x>0)
    S <- sum(x>0)
    E1 <- ifelse(q!=1, (1-qDest^(1-q))/(1-S^(1-q)), log(qDest)/log(S))
    E2 <- ifelse(q!=1, (1-qDest^(q-1))/(1-S^(q-1)), log(qDest)/log(S))
    E3 <- (qDest-1)/(S-1)
    E4 <- (1-1/qDest)/(1-1/S)
    E5 <- log(qDest)/log(S)
    if(q==0){
      p <- p[p>0]
      nu <- abs(p - (1/S))
      nu <- nu[nu > 0]
      sub1 <- (sum(log(abs(nu)))/sum(nu>0)-(log(1-1/S)+(1-S)*log(S))/S)
      E6 <- 1-exp(sub1)
    }else{
      p <- p[p>0]
      E6 <- 1-(sum(abs(p-1/S)^q)/((1-1/S)^q+(S-1)*S^(-q)))^(1/q)
    }

    #E6 <- ifelse(q=1, 1-sum(abs(p-1/S)^(1-q))/(abs(1-1/S)^(1-q)+)
    return(c(E1,E2,E3,E4,E5,E6))
  }
  out <- as.matrix(t(sapply(q.order, q_profile_evenness)))
  colnames(out) <- c("E1", "E2", "E3", "E4", "E5", "E6")
  out
}

qD <- function(p,q){
  p <- p[p>0]
  if(q!=1){
    (sum(p^q))^(1/(1-q))
  }else{
    exp(-sum(p*log(p)))
  }
}


#' \code{dis1}: computes the contribution of each species/node to the two types of dissimilarity measures.
#' @param x the species-by-assemblages abundance matrix or data.frame with species names as rownames and the colnames as assemblage/site names.
#' @param q a single value for the diversity order.
#' @param type a string: "tax" (taxonomic) or "phy" (phylogenetic).
#' @param type2 a string: "species" or "k"."species" means the contribution of each species/node to the two types of dissimilarity measures
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity).
#' "k" means the contribution of each assemblage/location/site to the two types of dissimilarity measures
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity). In the worked example, the contribution of each assemblage/stage is not computed.
#' @param tree the pylo object of the phylogenetic tree of all assemblages. Needed only if \code{type="phy"}.
#' @return a data frame tabulating the contribution of each species/node to the two types of dissimilarity measures: Jaccard-type (1-U_qN) and Sorensen-type (1-C_qN)
#' @examples
#' #Taxonomic analysis example
#' data(Alpine)
#' out0 <- dis1(x = Alpine, q = 0, type = "tax")
#' #Phylogenetic analysis example
#' data(Alpine)
#' data(tree_Alpine)
#' out0 <- dis1(x = Alpine, q = 0, type = "phy", tree = tree_Alpine)
#' @importFrom chaoUtility phylo2chaolabphy
#' @export
dis1 <- function(x, q, type = "tax", type2 = "species", tree = NULL){
  if(type2 == "species"){
    FUN <- rowSums
  }else{
    FUN <- colSums
  }
  if(type == "tax"){
    x <- as.matrix(x)
    x <- x[rowSums(x)>0, ]
    N <- ncol(x)
    zbar <- rowSums(x)/N
    x1 <- x[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    if(q==0){
      UqN <- FUN(x==0)/((N-1)*(sum(rowSums(x)>0)))
      CqN <- FUN(x==0)/((N-1)*(sum(apply(x, 2, function(i){sum(i>0)}))))
    }else if(q==2){
      UqN <- FUN((x1-zbar1)^2)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1-zbar1)^2)/((1-N^(1-q))*sum(x1^q))
    }else if(q!=1){
      UqN <- FUN((x1)^q-(zbar1)^q)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1)^q-(zbar1)^q)/((1-N^(1-q))*sum(x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(x1*log(x2), na.rm = T)/((sum(x)*log(N)))
      CqN <- UqN
    }
  }else{
    tree <- phylo2chaolabphy(tree)
    parts_nms <- rev(names(tree$parts))
    tree$parts <- lapply(length(tree$parts):1, function(i){tree$parts[[i]]})
    names(tree$parts) <- parts_nms
    Li <- c(tree$leaves, tree$nodes)
    cumtree = function(a, tree){
      a <- a[names(tree$leaves)]
      for(i in 1:length(tree$parts)){
        a[1+length(a)] <- sum(a[tree$parts[[i]]])
        names(a)[length(a)] <- names(tree$parts)[i]
      }
      a
    }
    ai <- apply(x, 2, cumtree, tree)
    wt <- apply(ai, 1, function(x1)(sum(x1))^q/sum(Li*rowSums(ai, na.rm = T)^q))
    N <- ncol(ai)
    zbar <- rowSums(ai)/N
    x1 <- ai[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    Li <- Li[zbar>0]
    T1 <- sum(rowSums(x1)*Li)
    if(q==0){
      if(type2 == "species"){
        rn <- nrow(x1)
        UqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))})/((N-1)*sum(Li))
        CqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))/((N-1)*sum(Li*rowSums(x1!=0)))})
      }else{
        UqN <- apply(x1, 2, function(x){sum(Li[x==0])})/((N-1)*sum(Li))
        CqN <- apply(x1, 2, function(x){sum(Li[x==0])/((N-1)*sum(Li*colSums(x1!=0)))})
      }

    }else if(q==2){
      UqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else if(q!=1){
      UqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(Li*x1*log(x2), na.rm = T)/(T1*log(N))
      CqN <- UqN
    }
  }

  # c(sum(UqN), sum(CqN))
  t(rbind(UqN, CqN))
}

#' \code{draw_dis_spe} plots the contribution of each species/node to dissimilarity (Jaccard-type dissimilarity and Sorensen-type dissimilarity).
#' @param x a merged table of output values with three columns corresponding to output for q = 0, 1, 2. See example for details.
#' @param title_name the title name of plot (Usually the Jaccard-type dissimilarity or Sorensen-type dissimilarity). See example for details.
#' @param type a string specifying the type of contribution: "tax" for taxonomic and "phy" for phylogenetic
#' @return a ggplot2 object showing the contribution of each species/node.
#' @examples
#' #Taxonomic analysis example
#' data(Alpine)
#' out0 <- dis1(x = Alpine, q = 0, type = "tax")
#' out1 <- dis1(x = Alpine, q = 1, type = "tax")
#' out2 <- dis1(x = Alpine, q = 0, type = "tax")
#' tax_UqN_r <- cbind(out0[, 1], out1[, 1], out2[, 1])
#' tax_CqN_r <- cbind(out0[, 2], out1[, 2], out2[, 2])
#' draw_dis_spe(tax_UqN_r, "Jaccard-type taxonomic dissimilarity")
#' draw_dis_spe(tax_CqN_r, "Sorensen-type taxonomic dissimilarity")
#' #Phylogenetic analysis example
#' data(Alpine)
#' data(tree_Alpine)
#' out0 <- dis1(x = Alpine, q = 0, type = "phy", tree = tree_Alpine)
#' out1 <- dis1(x = Alpine, q = 1, type = "phy", tree = tree_Alpine)
#' out2 <- dis1(x = Alpine, q = 2, type = "phy", tree = tree_Alpine)
#' phy_UqN_r <- cbind(out0[, 1], out1[, 1], out2[, 1])
#' phy_CqN_r <- cbind(out0[, 2], out1[, 2], out2[, 2])
#' draw_dis_spe(phy_UqN_r, "Jaccard-type phylogenetic dissimilarity", type = "phy")
#' draw_dis_spe(phy_CqN_r, "Sorensen-type phylogenetic dissimilarity", type = "phy")
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
draw_dis_spe <- function(x, title_name, type = "tax"){
  colnames(x) <- c("q = 0", "q = 1", "q = 2")
  x <- melt(x)
  g <- ggplot(x, aes(x = as.factor("Var1"), y = "value", fill = "Var2"))+
    geom_col(width = 0.2)+
    facet_grid(Var2~., scales = "free_y")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))+
    guides(fill=FALSE)+
    ggtitle(title_name)

  if(type == "tax"){
    g <- g +
      xlab("Species")+
      ylab("Species contribution")
  }else{
    g <-  g +
      xlab("Species/node")+
      ylab("Species/node contribution")
  }
  return(g)
}

