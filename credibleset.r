##Required functions for credible sets

prior.var <- function()
{
        ##theta's prior: N(0,W); return W
        return(0.04)
}

calc.BFs <- function(data)
{
        data$z = data$effect / data$se
        w = prior.var()

        ## BF = data(P(y | H1) / P(y | H0)
        data$ABF = sqrt(data$se^2 / (data$se^2 + w)) * exp(data$z^2 / 2 * w / (data$se^2 + w))
        return(data)

}

region.BF <- function(data)
{
        k = nrow(data)
        BFregion = 1/k * sum(data$ABF)
        return(BFregion)


}

posterior <- function(data)
{
        k = nrow(data)

        BFregion = region.BF(data)

        data$posterior = data$ABF / (k*BFregion)
        return(data)
}

get.posteriors <- function(data, region_l, region_u, region_chr)
{

        n = which(data$position >= region_l & data$position <= region_u & data$chr == region_chr)
        data = data[n,]

        n=which(data$se==0)
        if(length(n)!=0)
                data=data[-n,]
        data = calc.BFs(data)

        return(posterior(data))

}


credible.sets <- function(data, alpha, region_l, region_u, region_chr)
{
        data = get.posteriors(data, region_l, region_u, region_chr)
        data = data[sort.list(data$posterior, decreasing=TRUE),]

        sum_post=0
        credible=c()
        i = 1
        while(i <= nrow(data) & sum_post <= alpha)
        {
                sum_post = sum_post + data$posterior[i]
                credible = c(credible, as.character(data$rs_number[i]))
                i = i+1
        }

        results = data[which(data$rs_number %in% credible),]
       
        return(results)
}





##To obtain credible sets - full results with locus names and summarise for each locus
determine.credible.sets <- function(alpha, regions_file, input_file, path)
{

         
        ##give list of Locus, Chr, Start, Stop
        regions = read.table(regions_file, as.is=T, header=T)

        regions$n_snps = 0
        regions$distance = 0
        regions$start = 0
        regions$end = 0
        
        ##Main data input file

        data1 = read.table(input_file, as.is=T, header=T)
        ## remove missing stats
    data2<-subset(data1, effect!="NA" & se!="NA" & p!="NA")

        data<- data2[!duplicated(data2[c("rs_number")]),] ## remove any duplicates in the data.

        results<-vector("list", nrow(regions))
        

        for (i in 1:nrow(regions))
        {

                locus = as.character(paste(regions$Locus[i], regions$rs_number[i], sep="_"))
                region_chr = as.numeric(regions$chr[i])
                region_l = as.numeric(regions$region_start[i]) ## use b37 for both regions and summary stats.
                region_u = as.numeric(regions$region_end[i])

                        results[[i]] = credible.sets(data, alpha, region_l, region_u, region_chr)
                        results[[i]]$Locus=locus
                       

 ##summarise in regions
                        n_snps = nrow(results[[i]]) ## includes SNPs that are in multiple loci.
                        distance = max(results[[i]]$position) - min(results[[i]]$position) ## includes SNPs that are in multiple loci.

                
                        regions$n_snps[i] = n_snps

                        regions$distance[i] = distance

                        regions$start[i] = min(results[[i]]$position)
                        regions$end[i] = max(results[[i]]$position)
                        regions$total_posterior[i] = sum(results[[i]]$posterior)

        }

        resultCombined <-  do.call("rbind", results)

                  
        regionsOrdered <- regions[order(regions$chr, regions$position),]

        
        mylist<-list(regionsOrdered, resultCombined) 
        return(mylist)
}




##edit the three lines below and use in R command line for each analysis --> credibleSets=determine.credible.sets([alpha], "Input_with_regions_XXMb.txt", "Input_summary_stats.txt")

##credibleSets=determine.credible.sets(0.99, "GIANT_BMI_p-e5_Locke2015_chrBP_regions_1Mb.txt", "GIANT_2015_BMI_Men_SNP_gwas_mc_merge_nogc.tbl.uniq_chrBP_final.txt")

##write.table(credibleSets[[1]], file="GIANT_2015_BMI_EUR_MEN_CredibleSets0.99_summary_test.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

##write.table(credibleSets[[2]], file="GIANT_2015_BMI_EUR_MEN_CredibleSets0.99_results_test.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


