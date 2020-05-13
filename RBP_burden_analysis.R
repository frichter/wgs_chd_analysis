# Burden analysis function, to calculate enrichments between cases and controls in different genomic subsets.
burden_analysis <- function(desired_tests=NULL, statistical_test="fisher_exact_2s", peak_types="broadPeak", include_individual_RBPs=FALSE) {
    # Initialize some simple string variables.
    disease = "CHD"; control_name = "SSC"; expression_string_token = "HE"; expression_tissue_name = "E14.5 mouse heart"
    # Get case variables for the disease.
    cases <- standardize_colnames(get(paste0(tolower(disease)))); cases <- cases[,c(1:5,which(colnames(cases) %in% c("CADDgt10_Phred", "Chrom_hg19", "Position_hg19", "Chrom_hg38", "Position_hg38")))]
    if (statistical_test != "binomial") { n1 = nrow(cases)
    } else { n1 = get(paste0(toupper(disease), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    # Get control variables for the disease.
    controls <- standardize_colnames(get(paste0(tolower(control_name)))); controls <- controls[,c(1:5,which(colnames(controls) %in% c("CADDgt10_Phred", "Chrom_hg19", "Position_hg19", "Chrom_hg38", "Position_hg38")))]
    if (statistical_test != "binomial") { n0 = nrow(controls)
    } else { n0 = get(paste0(toupper(control_name), "_SAMPLE_COUNT")) } # statistical_test="binomial"
    
    shared_colnames <- intersect(colnames(cases), colnames(controls))
    cases <- cases[,shared_colnames]
    controls <- controls[,shared_colnames]
    
    # Figure out desired tests, if not specified in the function arguments.
    if(is.null(desired_tests)) {
        desired_variant_type_features <- get_features_by_group("variant_type") # is_snv, is_indel, and is_variant
        desired_region_type_features <- c("", "3'UTR", "TSS") # all, 3'UTR, TSS
        desired_gene_set_features <- c("", "constrained") # all, or only haploinsufficient genes
        desired_RBP_features <- c("", "RBP") # combined RBP peaks, with option to include individual RBP features too.
        if(include_individual_RBPs) { desired_RBP_features <- c(desired_RBP_features, get_features_by_group("RBP")) }
        desired_histone_mark_features <- c("", get_features_by_group("H3K36me3")) # H3K97me2 data much sparser and noisier, so is not included by default
        desired_histone_mark_features <- c("", desired_histone_mark_features[rowSums(sapply(peak_types, function(peak_type) grepl(peak_type, desired_histone_mark_features))) > 0])
        desired_CADD_features <- c("", get_features_by_group("CADD")) # variant harmfulness prediction
        
        desired_tests <- apply(expand.grid(desired_variant_type_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features, desired_histone_mark_features), 1, function(x) paste0(x[x != ""], collapse="+")) #desired_gene_set_features
        desired_tests <- unique(c(desired_tests, apply(expand.grid("is_snv", desired_CADD_features, desired_region_type_features, desired_gene_set_features, desired_RBP_features), 1, function(x) paste0(x[x != ""], collapse="+"))))
    }
    num_tests = length(desired_tests)
    all_features <- unique(unlist(strsplit(desired_tests, "\\+")))
    
    # Annotate cases and controls
    print("Annotating cases and controls...")
    cases <- annotate(cases, "nearest_gene"); controls <- annotate(controls, "nearest_gene"); 
    print("Annotating variants with features...")
    all_annotations <- annotate(rbind(cases, controls), all_features)
    cases <- all_annotations[1:nrow(cases),]
    controls <- all_annotations[(nrow(cases)+1):nrow(all_annotations),]
    
    # Function to return the formatted test to append to the tests data frame.
    append_test <- function(test_name, variant_type, d1, d0) { return(c(test_name, variant_type, nrow(d1), nrow(d0)) }
    
    # Add remaining tests to the data frame
    tests <- unfactorize(data.frame(t(sapply(1:num_tests, function(i) {
        print(paste0("Appending test ",desired_tests[i],"...[",i," / ",num_tests,"]"))
        test_elements <- unlist(strsplit(desired_tests[i], "\\+"))
        num_test_elements = length(test_elements)
        annotations_cases <- cases[,which(colnames(cases) %in% test_elements)]; annotations_controls <- controls[,which(colnames(controls) %in% test_elements)]
        if(num_test_elements > 1) { 
            hits_cases <- rowSums(annotations_cases == "Y") == ncol(annotations_cases)
            hits_controls <- rowSums(annotations_controls == "Y") == ncol(annotations_controls)
        } else { 
            hits_cases <- annotations_cases == "Y"
            hits_controls <- annotations_controls == "Y"
        }
        
        # Determine variant type
        if("is_indel" %in% test_elements) { variant_type = "indels"
        } else if("is_snv" %in% test_elements) { 
            variant_type = "SNVs";
            if(subsample_snvs) { 
                hits_cases <- hits_cases & (1:length(hits_cases) %in% sample(1:length(hits_cases), sum(cases$is_indel == "Y"), replace=FALSE))
                hits_controls <- hits_controls & (1:length(hits_controls) %in% sample(1:length(hits_controls), sum(controls$is_indel == "Y"), replace=FALSE)) 
            }
        } else { variant_type = "SNV+indel" }
        
        # Determine test name
        gene_set_feature = intersect(test_elements, desired_gene_set_features)
        region_feature = intersect(test_elements, desired_region_type_features)
        RBP_feature = intersect(test_elements, desired_RBP_features)
        histone_mark_feature = intersect(test_elements, desired_histone_mark_features)
        CADD_feature = intersect(test_elements, desired_CADD_features)
        if(length(gene_set_feature)>0) { test_name = gene_set_feature } else { test_name = NULL }
        if(length(CADD_feature)>0) { test_name = paste(c(test_name, CADD_feature), collapse=" ") }
        if(length(region_feature)>0) { test_name = paste(c(test_name, region_feature), collapse=" ") }
        if(length(RBP_feature)>0) { test_name = paste(c(test_name, RBP_feature), collapse=" ") }
        
        if(length(histone_mark_feature)>0) {
            histone_feature_split <- strsplit(histone_mark_feature, "\\.")[[1]]
            eid = histone_feature_split[1]
            tissue = get_roadmap_epigenome_names(eid); histone_mark = histone_feature_split[2]; peak_type = histone_feature_split[3] }
            test_name = paste(c(test_name, paste0(histone_mark, " in ", tissue, " (", eid, ", ", peak_type, ")")), collapse=" ")
            notes = get_histone_mark_function(histone_mark)
        } else { notes = "" }
        if(is.null(test_name) || test_name == "") {
            test_name = "overall" 
            notes = paste0(round(sum(hits_cases)/get(paste0(toupper(disease),"_SAMPLE_COUNT")),2), " variants/case, ", round(sum(hits_controls)/get(paste0(toupper(control_name),"_SAMPLE_COUNT")),2), " variants/control, two-sided binomial test")
        }
        
        return(append_test(test_name, variant_type, cases[hits_cases,], controls[hits_controls,], notes))
    }))))
    colnames(tests) <- c("test", "variants", "m1", "m0", "notes")
    overall_burden_tests <- tests$test == "overall"
    if(sum(overall_burden_tests)>0) { tests <- tests[c(which(overall_burden_tests),which(!overall_burden_tests)),] }
    
    tests$m1 <- as.numeric(tests$m1); tests$m0 <- as.numeric(tests$m0)
    test_enrichments <- c()
    test_p.values <- c()
    overall_norm_factors <- c(1, 1, 1); names(overall_norm_factors) <- c("SNV+indel", "SNVs", "indels")
    for(i in 1:nrow(tests)) {
        if (statistical_test != "poisson regression" && tests[i,1] == "overall") { 
            test_result <- binomial_test(tests$m1[i], tests$m0[i], get(paste0(toupper(disease),"_SAMPLE_COUNT")), get(paste0(toupper(control_name),"_SAMPLE_COUNT")), alternative="two.sided")
        } else if (statistical_test == "binomial") { norm_factor = overall_norm_factors[tests[i,2]]; test_result <- binomial_test(floor(tests$m1[i]/norm_factor), tests$m0[i], n1, n0, alternative="greater") 
        } else if (statistical_test == "poisson regression") { 
            if (tests[i,1] == "overall") {
                tests[i,6] <- gsub(", two-sided binomial test", "", tests[i,6])
                variants = tolower(tests[i,2]); if(variants == "snp+indel") { variants = "muts" }
                case_samples <- table(get(paste0("cases_",variants))$sample)
                control_samples <- table(get(paste0("controls_",variants))$sample)
            } else {
                case_samples <- table(paste0(sapply(unlist(strsplit(tests$cases_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
                control_samples <- table(paste0(sapply(unlist(strsplit(tests$controls_variant_hits[i],";")), function(x) strsplit(x,"_")[[1]][1])))
            }
            case_samples <- cbind(names(case_samples), case_samples); colnames(case_samples) <- c("sample", "variant_count")
            control_samples <- cbind(names(control_samples), control_samples);  colnames(control_samples) <- c("sample", "variant_count")
            case_paternal_ages <- get("chd_parental_age")
            case_paternal_ages <- case_paternal_ages[case_paternal_ages$Blinded.ID %in% cases_muts$sample,]
            control_paternal_ages <- get("ssc_parental_age")[,1:2]; control_paternal_ages <- control_paternal_ages[control_paternal_ages$Blinded.ID %in% controls_muts$sample,]
            case_dat <- merge(case_paternal_ages, case_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            case_dat$variant_count <- as.numeric(paste0(case_dat$variant_count)); case_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(case_dat$Paternal.Age.at.Proband.Birth))
            case_dat$variant_count[which(is.na(case_dat$variant_count))] <- 0; case_dat$Paternal.Age.at.Proband.Birth[which(is.na(case_dat$Paternal.Age.at.Proband.Birth))] <- mean(case_dat$Paternal.Age.at.Proband.Birth[!is.na(case_dat$Paternal.Age.at.Proband.Birth)])
            control_dat <- merge(control_paternal_ages, control_samples, by.x="Blinded.ID", by.y="sample", all.x=TRUE, all.y=TRUE)
            control_dat$variant_count <- as.numeric(paste0(control_dat$variant_count)); control_dat$Paternal.Age.at.Proband.Birth <- as.numeric(paste0(control_dat$Paternal.Age.at.Proband.Birth))
            control_dat$variant_count[which(is.na(control_dat$variant_count))] <- 0; control_dat$Paternal.Age.at.Proband.Birth[which(is.na(control_dat$Paternal.Age.at.Proband.Birth))] <- mean(control_dat$Paternal.Age.at.Proband.Birth[!is.na(control_dat$Paternal.Age.at.Proband.Birth)])
            variant_count_dat <- rbind(case_dat, control_dat)
            case_control <- rep(0, nrow(variant_count_dat)); case_control[1:nrow(case_dat)] <- 1; variant_count_dat <- cbind(variant_count_dat, case_control)
            
            m1 <- glm(variant_count ~ case_control + Paternal.Age.at.Proband.Birth, family="poisson", data=variant_count_dat)
            test_result <- new.env()
            test_result[["estimate"]] <- coef(summary(m1))[2,1]
            test_result[["p.value"]] <- coef(summary(m1))[2,4]
        } else if (statistical_test == "fisher_exact_2s") {
            if(tests$variants[i] == "SNV+indel") {
                n1_i = n1; n0_i = n0
            } else {
                if(tests$variants[i] == "SNVs") { n1_i = sum(cases$is_snv == "Y"); n0_i = sum(controls$is_snv == "Y")
                } else { n1_i = sum(cases$is_indel == "Y"); n0_i = sum(controls$is_indel == "Y") }
            }
            test_result <- fisher_exact_test(tests$m1[i], tests$m0[i], n1_i, n0_i, alternative="two.sided") # fisher_exact two-sided
        } else { test_result <- fisher_exact_test(tests$m1[i], tests$m0[i], n1, n0, alternative="greater") } # fisher_exact one-sided
        test_enrichments <- c(test_enrichments, test_result[["estimate"]])
        test_p.values <- c(test_p.values, test_result[["p.value"]])
    }
    burdens <- cbind(tests$test, tests$variants, tests$m1, tests$m0, test_enrichments, test_p.values, n1, n0, tests$notes)
    colnames(burdens) <- c("burden", "variants", "m1", "m0", "enrichment", "p.value", "n1", "n0", "notes")
    if(statistical_test == "poisson regression") { colnames(burdens)[5] <- "z_score" }
    
    # Return burden/enrichment results
    return(data.frame(burdens))
}