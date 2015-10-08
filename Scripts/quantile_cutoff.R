for (package in c("survival")) {
	if (!require(package, character.only=T, quietly=T)) {
		install.packages(package)
		library(package, character.only=T)
	}
}

genes_to_search = c("uc001ezg.2","IL17A|3605","IL17F|112744","CCL20|6364","CSF2|1437","DHCR7|1717","DHCR24|1718")
transcript_name = c("rorgt","il17a","il17f","ccl20", "csf2", "dhcr7", "dhcr24", "dhcr7-dhcr24")

all_transcripts_list = vector("list",length(genes_to_search))

output = "quantile_cutoff"
dir.create(output)
for (cancer in c("CRC")) {
	samples_isoforms = read.table(paste(cancer,"allsamples_isoforms.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
	samples_genes = read.table(paste(cancer,"allsamples_genes.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)

	samplemap = read.table(paste(cancer,"FILE_SAMPLE_MAP.txt",sep="/"),header=T,fill=T)
    clinicalstat = read.table(paste(cancer,"nationwidechildrens.org_clinical_patient.txt",sep="/"),header=T,fill=T,sep="\t")

    disease_status_vect = mat.or.vec(length(samples_isoforms[1,]),1)
    followup_vect = mat.or.vec(length(samples_isoforms[1,]),1)
    samplenames = colnames(samples_isoforms)

    normals = c()

    for (i in 1:length(samplemap[,1])) {
        barcode = strsplit(as.character(samplemap[i,2]),split="\\-")[[1]]
        regex = paste(barcode[1:3],collapse="-")
        found = grep(regex,clinicalstat[,"bcr_patient_barcode"])
        insert_spot = grep(as.character(samplemap[i,1]),samplenames)
        if (length(found)!=0) {
            followup_val = as.character(clinicalstat[found,"last_contact_days_to"])
            if (followup_val == "[Not Available]") {
                followup_val = as.character(clinicalstat[found,"death_days_to"])
            }
            if (as.numeric(substr(barcode[4],1,2))<=14 && as.numeric(substr(barcode[4],1,2))>=10) {
                if(length(insert_spot)!=0) {normals = rbind(normals, c(insert_spot, as.character(samplemap[i,2])))}
            }
            followup_vect[insert_spot] = followup_val
            if (clinicalstat[found,"vital_status"] == "Alive") {
                disease_status_vect[insert_spot] = 1
            }
            else if (clinicalstat[found,"vital_status"] == "Dead") {
                disease_status_vect[insert_spot] = 0
            }
        }
        else {
            followup_vect[insert_spot] = "[Not Available]"
        }
    }
	
	followup_tweak = which(!is.na(as.numeric(followup_vect)))
	
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

	for (m in 1:length(transcript_name)) {
		transcript = all_transcripts[m,]
		
		if (transcript_name[m] != "rorgt") {
			cutoff_type_x = c("0.25","0.5","0.75")
			cutoffs_x = c(quantile(rorgt,0.25), quantile(rorgt,0.5),quantile(rorgt,0.75))
			cutoff_type_y = c("0.25","0.5","0.75","0.95")
			cutoffs_y = c(quantile(transcript,0.25),quantile(transcript,0.5),quantile(transcript,0.75),quantile(transcript,0.95))
						
			transcript_x = rorgt[followup_tweak]
			transcript_y = transcript[followup_tweak]
			
			survival_input = cbind(as.numeric(disease_status_vect[followup_tweak]),as.numeric(followup_vect[followup_tweak])/30)

			rownames(survival_input) = seq(length(survival_input[,1]))
			colnames(survival_input) = c("disease_status","followup")

			for (i in 1:length(cutoffs_x)) {
				for (j in 1:length(cutoffs_y)) {

					x_Y_idx = intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y > cutoffs_y[j]))
					X_y_idx = intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y <= cutoffs_y[j]))
					x_y_idx = intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y <= cutoffs_y[j]))
					X_Y_idx = intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y > cutoffs_y[j]))
										
					if (length(x_Y_idx)>1 && length(X_y_idx)>1 && length(x_y_idx)>1 && length(X_Y_idx)>1) {
						tiff(filename=paste(output,paste(cancer,"_quadrants_rorgt_vs_",transcript_name[m],"_cutoff_x_",cutoff_type_x[i],"_cutoff_y_",cutoff_type_y[j],".tiff",sep=""),sep="/"))
					
						plot(transcript_x[x_Y_idx],transcript_y[x_Y_idx],xlab="rorgt",ylab=transcript_name[m],col="purple", xlim=c(min(transcript_x),(max(transcript_x))), ylim=c(min(transcript_y),max(transcript_y)))
						points(transcript_x[X_y_idx],transcript_y[X_y_idx],col="green")
						points(transcript_x[x_y_idx],transcript_y[x_y_idx],col="blue")
						points(transcript_x[X_Y_idx],transcript_y[X_Y_idx],col="red")

						dev.off()

						x_low_y_high = as.data.frame(survival_input[intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y > cutoffs_y[j])),])
				
						x_high_y_low = as.data.frame(survival_input[intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y <= cutoffs_y[j])),])
						
						x_low_y_low = as.data.frame(survival_input[intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y <= cutoffs_y[j])),])
						
						x_high_y_high = as.data.frame(survival_input[intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y > cutoffs_y[j])),])

						x_low_y_high.surv = survfit(Surv(followup, disease_status)~1, data = x_low_y_high, type="kaplan-meier", conf.type="none")
						x_high_y_low.surv = survfit(Surv(followup, disease_status)~1, data = x_high_y_low, type="kaplan-meier", conf.type="none")
						x_low_y_low.surv = survfit(Surv(followup, disease_status)~1, data = x_low_y_low, type="kaplan-meier", conf.type="none")
						x_high_y_high.surv = survfit(Surv(followup, disease_status)~1, data = x_high_y_high, type="kaplan-meier", conf.type="none")
											
						tiff(filename=paste(output,paste(cancer,"_survival_expression_rorgt_vs_",transcript_name[m],"_cutoff_x_",cutoff_type_x[i],"_cutoff_y_",cutoff_type_y[j],".tiff",sep=""),sep="/"))

						plot(x_low_y_high.surv, col="purple", xlab=paste("rorgt cutoff: ",cutoff_type_x[i],", ",transcript_name[m]," cutoff: ",cutoff_type_y[j], sep=""), main=paste("rorgt vs ",transcript_name[m]," Survival",sep=""), xlim=c(0,max(survival_input[,2])))
						lines(x_high_y_low.surv, col="green")
						lines(x_low_y_low.surv, col="blue")
						lines(x_high_y_high.surv, col="red")
						
						legend(80,1,inset=0, c(paste("rorgt low, ",transcript_name[m]," high",sep=""),paste("rorgt high, ",transcript_name[m]," low",sep=""),paste("rorgt low, ",transcript_name[m]," low",sep=""),paste("rorgt high, ",transcript_name[m]," high",sep="")), fill = c("purple","green","blue","red"))

						dev.off()
						
					}
					
				}
			}
			
		}
		
	}
	print (cancer)

}
