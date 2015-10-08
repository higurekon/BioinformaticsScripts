for (package in c("survival")) {
	if (!require(package, character.only=T, quietly=T)) {
		install.packages(package)
		library(package, character.only=T)
	}
}

genes_to_search = c("uc001ezg.2","IL17A|3605","IL17F|112744","CCL20|6364","CSF2|1437","DHCR7|1717","DHCR24|1718")
transcript_name = c("rorgt","il17a", "il17f", "ccl20", "csf2", "dhcr7", "dhcr24")
all_transcripts_list = vector("list",length(genes_to_search))

output = "correlation_of_cancer_with_rorgt"
dir.create(output)
for (cancer in c("CRC")) {
	samples_isoforms = read.table(paste(cancer,"allsamples_isoforms.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
	samples_genes = read.table(paste(cancer,"allsamples_genes.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
	
	# isoforms_norm = calcNormFactors(samples_isoforms)
	# genes_norm = calcNormFactors(samples_genes)
	
	# samples_isoforms = sweep(samples_isoforms, 2, isoforms_norm)
	# samples_genes = sweep(samples_genes, 2, genes_norm)
	
	i = 1
	for (gene in genes_to_search) {
		if (length(which(rownames(samples_isoforms)==gene)!=0)) {
			all_transcripts_list[[i]] = c(all_transcripts_list[[i]], as.numeric(t(samples_isoforms[gene,]))) 
		}
		else {
			all_transcripts_list[[i]] = c(all_transcripts_list[[i]], as.numeric(t(samples_genes[gene,])))
		}
		i = i+1
	}

	all_transcripts = do.call(rbind,all_transcripts_list)
	
	for (i in 1:length(transcript_name)) {
		tiff(filename=paste(output,paste(cancer,transcript_name[i],"expression_levels.tiff",sep="_"),sep="/"))
		hist(all_transcripts[i,], xlab = "expression levels", main = transcript_name[i])
		dev.off()
		
		if (transcript_name[i]!="rorgt") {
			nonzero_samples = intersect(which(rorgt>0), which(all_transcripts[i,]>0))
			
			if (length(nonzero_samples) != 0) {
				tiff(filename=paste(output,paste(cancer,"_simple-correlation_rorgt_v_",transcript_name[i],".tiff",sep=""),sep="/"))
				plot(rorgt[nonzero_samples],all_transcripts[i,][nonzero_samples], xlab="rorgt", ylab=transcript_name[i], main=paste("correlation: ",cor(rorgt[nonzero_samples],all_transcripts[i,][nonzero_samples])))
				abline(lm(all_transcripts[i,][nonzero_samples]~rorgt[nonzero_samples]),col="red")
				dev.off()
			}
			
		}
		
	}
	
	print(cancer)
}

