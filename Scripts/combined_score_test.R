for (package in c("survival")) {
	if (!require(package, character.only=T, quietly=T)) {
		install.packages(package)
		library(package, character.only=T)
	}
}

genes_to_search = c("uc001ezg.2","IL17A|3605","IL17F|112744","CCL20|6364","CSF2|1437","IL22|50616","IL23R|149233")
transcript_name = c("rorgt","il17a","il17f","ccl20", "csf2", "il22", "il23r")

all_transcripts_list = vector("list",length(genes_to_search))

all_followup = c()
all_disease_status = c()

output = "combined_score_test"
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
	
	all_followup = c(all_followup, followup_vect)
	all_disease_status = c(all_disease_status, disease_status_vect)
	
}

followup_tweak = which(!is.na(as.numeric(all_followup)))

all_transcripts = do.call(rbind,all_transcripts_list)

for (i in 1:length(all_transcripts[,1])) {
	tiff(filename=paste(output,paste("log_histogram_of_",transcript_name[i],"_expression.tiff",sep=""),sep="/"))
	log2_all_transcripts = log2(all_transcripts[i,])
	hist(log2_all_transcripts,xlab="Log Expression",main=transcript_name[i])
	dev.off()
}

# manual_cutoffs = c(32, 1, 1, 256, 2, 1024, 4096, 0.25)


all_transcripts_centered = apply(all_transcripts,2,function(x) x/median(x))
score = apply(all_transcripts_centered,2,sum)
all_rorgt = all_transcripts[1,][followup_tweak]
transcript = score[followup_tweak]

transcript_x = all_rorgt
transcript_y = transcript[which(!is.nan(transcript_y))]
		
cutoff_types = c("0.25","0.3","0.4","0.5","0.6","0.7","0.75")
cutoffs_x = c(quantile(transcript_x,0.25), quantile(transcript_x,0.3), quantile(transcript_x,0.4), quantile(transcript_x,0.5), quantile(transcript_x,0.6), quantile(transcript_x,0.7), quantile(transcript_x,0.75))
cutoffs_y = c(quantile(transcript_y,0.25), quantile(transcript_y,0.3), quantile(transcript_y,0.4), quantile(transcript_y,0.5), quantile(transcript_y,0.6), quantile(transcript_y,0.7), quantile(transcript_y,0.75))
					
survival_input = cbind(as.numeric(all_disease_status[followup_tweak]),as.numeric(all_followup[followup_tweak])/30)

rownames(survival_input) = seq(length(survival_input[,1]))
colnames(survival_input) = c("disease_status","followup")

for (i in 1:length(cutoffs_x)) {
	for (j in 1:length(cutoffs_y)) {

		x_Y_idx = intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y > cutoffs_y[j]))
		X_y_idx = intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y <= cutoffs_y[j]))
		x_y_idx = intersect(which(transcript_x <= cutoffs_x[i]),which(transcript_y <= cutoffs_y[j]))
		X_Y_idx = intersect(which(transcript_x > cutoffs_x[i]),which(transcript_y > cutoffs_y[j]))
							
		if (length(x_Y_idx)>1 && length(X_y_idx)>1 && length(x_y_idx)>1 && length(X_Y_idx)>1) {
			tiff(filename=paste(output,paste("quadrants_rorgt_vs_score-of-targets_cutoff_x_",cutoff_types[i],"_cutoff_y_",cutoff_types[j],".tiff",sep=""),sep="/"))
		
			plot(transcript_x[x_Y_idx],transcript_y[x_Y_idx],col="purple", xlab="rorgt", ylab="score-of-targets", xlim=c(min(transcript_x),(max(transcript_x))), ylim=c(min(transcript_y),max(transcript_y)))
			points(transcript_x[X_y_idx],transcript_y[X_y_idx],col="green")
			points(transcript_x[x_y_idx],transcript_y[x_y_idx],col="blue")
			points(transcript_x[X_Y_idx],transcript_y[X_Y_idx],col="red")

			dev.off()

			x_low_y_high = as.data.frame(survival_input[x_Y_idx,])
	
			x_high_y_low = as.data.frame(survival_input[X_y_idx,])
			
			x_low_y_low = as.data.frame(survival_input[x_y_idx,])
			
			x_high_y_high = as.data.frame(survival_input[X_Y_idx,])

			x_low_y_high.surv = survfit(Surv(followup, disease_status)~1, data = x_low_y_high, type="kaplan-meier", conf.type="none")
			x_high_y_low.surv = survfit(Surv(followup, disease_status)~1, data = x_high_y_low, type="kaplan-meier", conf.type="none")
			x_low_y_low.surv = survfit(Surv(followup, disease_status)~1, data = x_low_y_low, type="kaplan-meier", conf.type="none")
			x_high_y_high.surv = survfit(Surv(followup, disease_status)~1, data = x_high_y_high, type="kaplan-meier", conf.type="none")
								
			tiff(filename=paste(output,paste("survival_expression_rorgt_vs_score-of-targets","_cutoff_x_",cutoff_types[i],"_cutoff_y_",cutoff_types[j],".tiff",sep=""),sep="/"))

			plot(x_low_y_high.surv, col="purple", xlab=paste("rorgt cutoff: ",cutoff_types[i],", score-of-targets cutoff: ",cutoff_types[j], sep=""), main="rorgt vs score-of-targets survival", xlim=c(0,max(survival_input[,2])))
			lines(x_high_y_low.surv, col="green")
			lines(x_low_y_low.surv, col="blue")
			lines(x_high_y_high.surv, col="red")
			
			legend(80,1,inset=0, c("rorgt low, score-of-targets high","rorgt high, score-of-targets low","rorgt low, score-of-targets low","rorgt high, score-of-targets high"), fill = c("purple","green","blue","red"))

			dev.off()
			
		}
		
	}
}
