for (package in c("survival")) {
	if (!require(package, character.only=T, quietly=T)) {
		install.packages(package)
		library(package, character.only=T)
	}
}

all_rorgt = c()
all_il17a = c()

output = "correlation_cutoff_analysis"
dir.create(output)
for (cancer in c("CRC")) {
	samples_isoforms = read.table(paste(cancer,"allsamples_isoforms.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
	samples_genes = read.table(paste(cancer,"allsamples_genes.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)

	rorgt = as.numeric(t(samples_isoforms["uc001ezg.2",]))
	il17a = as.numeric(t(samples_genes["IL17A|3605",]))
	il17f = as.numeric(t(samples_genes["IL17F|112744",]))
	ccl20 = as.numeric(t(samples_genes["CCL20|6364",]))
	csf2 = as.numeric(t(samples_genes["CSF2|1437",]))

	all_rorgt = c(all_rorgt, rorgt)
	all_il17a = c(all_il17a, il17a)
}

tiff(filename=paste(output,"simple-correlation_rorgt_v_il17a_allcancers.tiff",sep="/"))
plot(all_rorgt,all_il17a,xlab="all cancers rorgt expression", ylab="all cancers il17 expression", main=paste("correlation: ",cor(all_rorgt,all_il17a)))
abline(lm(all_il17a~all_rorgt),col="red")
dev.off()

resids = residuals(lm(all_il17a~all_rorgt))

cutoff_percentile = 0.1
cutoff_closeness = 0.25 
cutoff_expr = quantile(all_rorgt, cutoff_percentile)
cutoff_actv = quantile(abs(resids), cutoff_closeness)

all_rorgt_active_idx = which(abs(resids) < cutoff_actv)
all_rorgt_inactive_idx = which(abs(resids) >= cutoff_actv)


starter = 0
for (cancer in c("CRC")) {
	samples_isoforms = read.table(paste(cancer,"allsamples_isoforms.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
	samples_genes = read.table(paste(cancer,"allsamples_genes.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)

	samplemap = read.table(paste(cancer,"FILE_SAMPLE_MAP.txt",sep="/"),header=T,fill=T)
    clinicalstat = read.table(paste(cancer,"nationwidechildrens.org_clinical_patient.txt",sep="/"),header=T,fill=T,sep="\t")

    disease_status = mat.or.vec(length(samples_isoforms[1,]),1)
    followup = mat.or.vec(length(samples_isoforms[1,]),1)
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
            followup[insert_spot] = followup_val
            if (clinicalstat[found,"vital_status"] == "Alive") {
                disease_status[insert_spot] = 1
            }
            else if (clinicalstat[found,"vital_status"] == "Dead") {
                disease_status[insert_spot] = 0
            }
        }
        else {
            followup[insert_spot] = "[Not Available]"
        }
    }
	
	followup_tweak = which(!is.na(as.numeric(followup)))
	
	rorgt = as.numeric(t(samples_isoforms["uc001ezg.2",]))
	il17a = as.numeric(t(samples_isoforms["IL17A|3605",]))
	il17f = as.numeric(t(samples_genes["IL17F|112744",]))
	ccl20 = as.numeric(t(samples_genes["CCL20|6364",]))

	all_transcripts = rbind(rorgt, il17a, il17f, ccl20)

	transcript_names = c("rorgt","il17a","il17f","ccl20")
	
	sterol_ratio = as.numeric(t(samples_genes["DHCR7|1717",]))/as.numeric(t(samples_genes["DHCR24|1718",]))
	
	rorgt_active_idx = all_rorgt_active_idx[which(all_rorgt_active_idx <= starter+length(samples_isoforms))]
	rorgt_inactive_idx = all_rorgt_inactive_idx[which(all_rorgt_inactive_idx <= starter+length(samples_isoforms))]
	
	active_tweak = intersect(rorgt_active_idx, followup_tweak)
	inactive_tweak = intersect(rorgt_inactive_idx, followup_tweak)
		
	for (i in 1:4) {

	
		transcript = all_transcripts[i,]

		########################################## Expression
		
		transcript_expr = cbind(as.numeric(as.matrix(transcript)[followup_tweak]),as.numeric(disease_status[followup_tweak]),as.numeric(followup[followup_tweak])/30)

		rownames(transcript_expr) = seq(length(transcript_expr[,1]))
		colnames(transcript_expr) = c("expression","disease_status","followup")

		transcript_high_expr = as.data.frame(transcript_expr[which(transcript_expr[,1] >= cutoff_expr),])
		transcript_low_expr = as.data.frame(transcript_expr[which(transcript_expr[,1] < cutoff_expr),])
		
		high_good = F
		low_good = F
		
		if (length(which(transcript_expr[,1] >= cutoff_expr))>10) {
			high.surv = survfit(Surv(transcript_high_expr[,"followup"], transcript_high_expr[,"disease_status"])~1, data = transcript_high_expr, type="kaplan-meier", conf.type="none")
			high_good = T
		}
		
		if (length(which(transcript_expr[,1] < cutoff_expr))>10) {
			low.surv = survfit(Surv(transcript_low_expr[,"followup"], transcript_low_expr[,"disease_status"])~1, data = transcript_low_expr, type="kaplan-meier", conf.type="none")
			low_good = T
		}


		if (high_good && low_good) {
			firstplot = high.surv
			firstplot_col = "red"
			secondplot = low.surv
			secondplot_col = "blue" 

			tiff(filename=paste(output,paste("survival_expression_",transcript_names[i],"_cutoff_",cutoff_percentile,".tiff",sep=""),sep="/"))

			plot(firstplot, col=firstplot_col, main=transcript_names[i], xlim=c(0,max(transcript_expr[,3])))
			lines(secondplot, col=secondplot_col)

			dev.off()
		}
		
		print ("expression done")
		########################################## Active
		
		transcript_actv = as.data.frame(cbind(as.numeric(as.matrix(transcript)[active_tweak]),as.numeric(disease_status[active_tweak]),as.numeric(followup[active_tweak])/30))

		transcript_inactv = as.data.frame(cbind(as.numeric(as.matrix(transcript)[inactive_tweak]),as.numeric(disease_status[inactive_tweak]),as.numeric(followup[inactive_tweak])/30))
		
		rownames(transcript_actv) = seq(length(transcript_actv[,1]))
		colnames(transcript_actv) = c("expression","disease_status","followup")

		rownames(transcript_inactv) = seq(length(transcript_inactv[,1]))
		colnames(transcript_inactv) = c("expression","disease_status","followup")

		high_good = F
		low_good = F
		if (length(active_tweak)>10) {
			high.surv = survfit(Surv(transcript_actv[,"followup"], transcript_actv[,"disease_status"])~1, data = transcript_actv, type="kaplan-meier", conf.type="none")
			high_good = T
		}
		if (length(inactive_tweak)>10) {
			low.surv = survfit(Surv(transcript_inactv[,"followup"], transcript_inactv[,"disease_status"])~1, data = transcript_inactv, type="kaplan-meier", conf.type="none")
			low_good = T
		}
		
		if (high_good && low_good) {
			firstplot = high.surv
			firstplot_col = "red"
			secondplot = low.surv
			secondplot_col = "blue" 

			tiff(filename=paste(output,paste("survival_active_",transcript_names[i],"_minresidual_",cutoff_actv,".tiff",sep=""),sep="/"))

			plot(firstplot, col=firstplot_col, main=transcript_names[i], xlim=c(0,max(transcript_expr[,3])))
			lines(secondplot, col=secondplot_col)

			dev.off()
		}
		print ("active done")
		
	}
	starter = starter + length(samples_isoforms)
}
