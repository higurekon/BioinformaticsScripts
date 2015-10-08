for (package in c("OptimalCutpoints","survival","gplots")) {
	if (!require(package, character.only=T, quietly=T)) {
		install.packages(package)
		library(package, character.only=T)
	}
}

output = "immune_compendium_clustering"
dir.create(output)
for (cancer in c("CRC")) {

    samples_isoforms = read.table(paste(cancer,"allsamples_isoforms.txt",sep="/"),header=T,check.names=F,stringsAsFactors=F)
    samplemap = read.table(paste(cancer,"FILE_SAMPLE_MAP.txt",sep="/"),header=T,fill=T)
    clinicalstat = read.table(paste(cancer,"nationwidechildrens.org_clinical_patient.txt",sep="/"),header=T,fill=T,sep="\t")
    samples_genes = read.table(paste(cancer,"allsamples_genes.txt",sep="/"),header=T,check.names=F)
    immunegenelist = read.table("genelist.txt",fill=T)

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

   
    ########################################################## Immune Compendium
    
    print("Immune Compendium...")
    
    rorgt = samples_isoforms["uc001ezg.2",]

    followup_tweak = which(!is.na(as.numeric(followup)))

    cutpoint_input = cbind(as.numeric(as.matrix(rorgt)[followup_tweak]),as.numeric(disease_status[followup_tweak]),as.numeric(followup[followup_tweak]))

    rownames(cutpoint_input) = seq(length(cutpoint_input[,1]))
    colnames(cutpoint_input) = c("expression","disease_status","followup")
    
    optimal.cutpoint.geneset <- optimal.cutpoints(X="expression", status="disease_status", tag.healthy=1, methods="MinPvalue", data=as.data.frame(cutpoint_input), pop.prev=NULL, categorical.cov=NULL, ci.fit=T, conf.level=0.95, trace=T)

    cutpoint_output = optimal.cutpoint.geneset$MinPvalue$Global$optimal.cutoff$cutoff
	if (length(samples_genes) < 100) {cutpoint_output = quantile(as.numeric(as.matrix(rorgt)),0.25)}

    immunegenes = c()
    for (i in 1:length(samples_genes[,1])) {
        genename = strsplit(as.character(rownames(samples_genes)[i]),"\\|")[[1]][1]
        index = which(genename == as.character(immunegenelist[,1]))
        if (length(index)!=0) {
            immunegenes = rbind(immunegenes, samples_genes[i,])
        }
    }

    indata = apply(immunegenes,c(1,2),as.numeric)
    samples_mat = as.matrix(samples_genes)
    indata_rorgt_pos = indata[,which(as.numeric(samples_mat["RORC|6097",]) > cutpoint_output),drop=F]
    indata_rorgt_neg = indata[,which(as.numeric(samples_mat["RORC|6097",]) <= cutpoint_output),drop=F]

    breaks = seq(quantile(indata_rorgt_pos, 0.25), quantile(indata_rorgt_pos, 0.75), by=0.2)
    mycol <- colorpanel(n=length(breaks)-1,low="red",mid="black",high="green")

    if (!is.null(dim(indata_rorgt_pos)) && dim(indata_rorgt_pos)[1] > 1 && dim(indata_rorgt_pos)[2] > 1) {
        tiff(filename=paste(output,"immunecompendium_rorgt_pos.tiff",sep="/"), width=25000, height=4000, res=300)
        heatmap.2(t(indata_rorgt_pos),
            key = T,
            symkey = T,
            col = mycol,
            breaks = breaks,
            margins = c(10,10),
            lhei = c(0.3,1),
            lwid = c(0.3,0.97),
            trace = 'none',
            Rowv = T)
        dev.off()
    }
    
    breaks = seq(quantile(indata_rorgt_neg, 0.25), quantile(indata_rorgt_neg, 0.75), by=0.2)
    mycol <- colorpanel(n=length(breaks)-1,low="red",mid="black",high="green")

    if (!is.null(dim(indata_rorgt_neg)) && dim(indata_rorgt_neg)[1] > 1 && dim(indata_rorgt_neg)[2] > 1) {
        tiff(filename=paste(output,"immunecompendium_rorgt_neg.tiff",sep="/"), width=25000, height=4000, res=300)
        heatmap.2(t(indata_rorgt_neg),
            key = T,
            symkey = T,
            col = mycol,
            breaks = breaks,
            margins = c(10,10),
            lhei = c(0.3,1),
            lwid = c(0.3,0.97),
            trace = 'none',
            Rowv = T)
        dev.off()
    }

	print(paste("###",cancer," is done",sep=""))
}
