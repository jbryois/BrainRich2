#Perform linear regression testing whether mean specificity of gene set is higher than average.
reg <- function(d,pathway) {
    
    # Scale data (substracts mean and divide by standard deviation)
    d[-c(1,2,3)] <- apply(d[-c(1,2,3)],2,scale)
    
    # Select column with matching names
    n_genes_symbol <- sum(pathway[[1]]%in%d$Gene)
    n_genes_ensembl <- sum(pathway[[1]]%in%d$gene_id)
    n_genes_entrez<- sum(pathway[[1]]%in%d$entrez_id)
    
    gene_column <- which.max(c(n_genes_symbol,n_genes_ensembl,n_genes_entrez))

    # Only keep genes in Pathway 
    d <- d %>% filter(.[[gene_column]]%in%pathway[[1]])
    
    # Gather data
    d <- gather(d,Lvl5,spe_10k_z,-Gene,-gene_id,-entrez_id)
    
    res <- tidy(lm(spe_10k_z ~0 + Lvl5,data=d)) %>% mutate(term=gsub("Lvl5","",term)) %>%
        rename(p=p.value) %>%
        mutate(p=ifelse(estimate>0,p/2,(1-p/2))) %>%
        mutate(fdr=p.adjust(p,method="fdr")) %>%
        arrange(p) %>%
        mutate(Significance=ifelse(fdr<0.05,"5% FDR", "Not Significant")) %>%
        mutate(Lvl5=factor(term,levels=rev(term)))
}

#Perform Fisher's exact test to test whether the gene set is enriched among the 10% most specific genes
fisher_test <- function(d,pathway) {
    
    d <- d %>% gather(Lvl5,spe_10k,-Gene,-gene_id,-entrez_id) %>% 
        mutate(spe_10k_top10=ifelse(spe_10k>=quantile(spe_10k,0.9),1,0)) %>% 
        select(-spe_10k)
    
    pathway <- pathway  %>% 
        mutate(Pathway=1) %>% 
        select(1,Pathway)
    
    # Select column with matching names
    n_genes_symbol <- sum(pathway[[1]]%in%d$Gene)
    n_genes_ensembl <- sum(pathway[[1]]%in%d$gene_id)
    n_genes_entrez<- sum(pathway[[1]]%in%d$entrez_id)
    
    gene_column <- which.max(c(n_genes_symbol,n_genes_ensembl,n_genes_entrez))
    colnames(pathway)[1] <- colnames(d)[gene_column]
    
    # Add Pathway information to specificity data
    d <- left_join(d,pathway,by=colnames(d)[gene_column])
    d <- mutate(d,Pathway=ifelse(is.na(Pathway),0,1))
    
    fisher <- d %>% group_by(Lvl5) %>% nest() %>%
        mutate(
            fit = map(data, ~ fisher.test(.x$spe_10k_top10, .x$Pathway,alternative="greater")),
            tidied = map(fit, tidy)
        ) %>% ungroup() %>% 
        unnest(tidied)  %>%
        select(Lvl5,estimate,p.value) %>% 
        rename(p=p.value,odds_ratio=estimate) %>%
        mutate(fdr=p.adjust(p,method="fdr")) %>%
        arrange(p) %>%
        mutate(Significance=ifelse(fdr<0.05,"5% FDR", "Not Significant")) %>%
        mutate(Lvl5=factor(Lvl5,levels=rev(Lvl5)))
}

# Uses EWCE to test for cell type enrichment
ewce <- function(d,gene_set,name=NULL,number_of_iteration=9999) {
    d <- d %>% gather(Lvl5,spe_10k,-Gene,-gene_id,-entrez_id)
    
    # Select column with matching names
    n_genes_symbol <- sum(gene_set[[1]]%in%d$Gene)
    n_genes_ensembl <- sum(gene_set[[1]]%in%d$gene_id)
    n_genes_entrez<- sum(gene_set[[1]]%in%d$entrez_id)
    
    gene_column <- which.max(c(n_genes_symbol,n_genes_ensembl,n_genes_entrez))
    gene_column_name <- colnames(d)[gene_column]
    colnames(gene_set)[1] <- gene_column_name
    proportion <- d %>% filter(!is.na((get(gene_column_name)))) %>% select(gene_column,Lvl5,spe_10k) %>% spread(Lvl5,spe_10k)
    
    # Add Pathway information to specificity data
    
    colnames(gene_set)[1] <- gene_column_name
    proportion_gene_set <- proportion[proportion[[gene_column_name]]%in%gene_set[[gene_column_name]],]
    
    cat(c("Number of protein coding genes in dataset",nrow(proportion),"\n"))
    
    #Print length of gene set
    cat(c("Length of gene set original",nrow(gene_set),"\n"))
    
    #Print number of genes in the gene set and also with a 1to1 ortholog
    cat(c("Length of gene set with a match in dataset",nrow(proportion_gene_set),"\n"))
    
    #Average proportion in each cell type
    mean_proportion <- apply(proportion_gene_set[-1],2,mean)
    
    ######
    #Permutation
    ######
    
    mean_proportion_bootstrap_df <- matrix(ncol=ncol(proportion_gene_set)-1,nrow=number_of_iteration)
    for (i in 1:number_of_iteration){
        mean_proportion_bootstrap_df[i,] <- apply(proportion[sample(nrow(proportion),nrow(proportion_gene_set),replace=F),][-1],2,mean)
        if(i%%1000==0){
            cat(i)
            cat("\n")
        }
    }
    
    #Get Pvalue
    
    pvalues <- vector("numeric", ncol(proportion_gene_set)-1)
    for (i in 1:ncol(mean_proportion_bootstrap_df)){
        number_null_more_extreme <- length(which(mean_proportion_bootstrap_df[,i] >= mean_proportion[i]))
        pvalues[i] <- (number_null_more_extreme+1)/(nrow(mean_proportion_bootstrap_df)+1)
    }
    names(pvalues) <- colnames(proportion_gene_set)[-1]
    
    #Get Z-score
    
    sd_boot <- apply(mean_proportion_bootstrap_df,2,sd)
    mean_boot <- apply(mean_proportion_bootstrap_df,2,mean)
    z_scores <- (mean_proportion-mean_boot)/sd_boot
    
    results <- cbind(pvalues,z_scores) %>% as.data.frame() %>% rownames_to_column(var = "Lvl5")  %>% arrange(pvalues,-z_scores) %>%
        rename(p=pvalues) %>% 
        mutate(fdr=p.adjust(p,method="fdr")) %>%
        mutate(Significance=ifelse(fdr<0.05,"5% FDR", "Not Significant")) %>%
        mutate(Lvl5=factor(Lvl5,levels=rev(Lvl5)))
}

#Plot enrichment results
plot_results <- function(results,plot_type){
    if(nrow(results)> 60) {
        results <- arrange(results,p) %>% head(60)
    }
    if(plot_type=="Effect Sizes"){
        if("z_scores"%in%colnames(results)){
            return(ggplot(results,aes(Lvl5,z_scores,fill=Significance)) + geom_col() + coord_flip() + theme_classic() + xlab("") + ylab("Z-scores"))
        }
        if("estimate"%in%colnames(results)){
            return(ggplot(results,aes(Lvl5,estimate,col=Significance)) + geom_point() + geom_errorbar(aes(ymax = estimate + std.error, ymin = estimate - std.error)) + coord_flip() + theme_classic() + xlab("") + ylab("Beta"))
        }
        if("odds_ratio"%in%colnames(results)){
            return(ggplot(results,aes(Lvl5,odds_ratio,fill=Significance)) + geom_col() + coord_flip() + theme_classic() + xlab("") + ylab("Odds ratio") + geom_hline(yintercept = 1))
        }
    }
    if(plot_type=="-log10(pvalue)"){
        return(ggplot(results,aes(Lvl5,-log10(p),fill=Significance)) + geom_col() + coord_flip() + theme_classic() + xlab("") + ylab(expression('-log'[10]*'(pvalue)')))
    }
}

# Parse dataset name
parse_dataset_name <- function(dataset_path){
    name <- gsub("Data/","", dataset_path)
    name <- gsub(".1to1.norm.txt.gz","", name)
    name <- gsub(".norm.txt.gz","", name)
    name <- gsub(".all.norm.txt.gz","", name)
}

# Heatmap
plot_heatmap <- function(d,pathway){
    
    d <- d %>% gather(Lvl5,spe_10k,-Gene,-gene_id,-entrez_id)
    
    pathway <- pathway  %>% 
        mutate(Pathway=1) %>% 
        select(1,Pathway)
    
    # Select column with matching names
    n_genes_symbol <- sum(pathway[[1]]%in%d$Gene)
    n_genes_ensembl <- sum(pathway[[1]]%in%d$gene_id)
    n_genes_entrez<- sum(pathway[[1]]%in%d$entrez_id)
    
    gene_column <- which.max(c(n_genes_symbol,n_genes_ensembl,n_genes_entrez))
    colnames(pathway)[1] <- colnames(d)[gene_column]
    
    # Add Pathway information to specificity data
    d <- left_join(d,pathway,by=colnames(d)[gene_column])
    d <- mutate(d,Pathway=ifelse(is.na(Pathway),0,1))
    
    d <- d %>% group_by(Lvl5) %>% mutate(spe_10k_bin = ntile(spe_10k,10)) %>% filter(Pathway==1) %>% ungroup()
    
    d <- d %>% select(-Pathway,-spe_10k) %>% spread(Lvl5,spe_10k_bin) %>% select(-Gene,-gene_id,-entrez_id) %>% t(.)
    if (nrow(d)<100){
        return(pheatmap::pheatmap(d))
    }
    if (nrow(d)>=100){
        return(pheatmap::pheatmap(d,fontsize=3))
        
    }
}