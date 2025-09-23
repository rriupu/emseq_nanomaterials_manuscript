test_two_treatments <- function(meth_obj, group1, group2, threads, meth_difference = 2, qvalue = 0.5) {
    # Get treatments vector from meth_obj
    treatments <- meth_obj@treatment
    sample_ids <- meth_obj@sample.ids

    # Find indices of samples belonging to either group1 or group2
    selected_samples <- which(treatments %in% c(group1, group2))

    # Subset methylBase object to only these samples using reorganize
    meth_sub <- reorganize(meth_obj, sample.ids = sample_ids[selected_samples], treatment = treatments[selected_samples])

    # Extract donor info from sample IDs of subset
    sample_id_sub <- meth_sub@sample.ids
    donor_sub <- sub(".*_(D[0-9]+)$", "\\1", sample_id_sub)
    donor_sub <- factor(donor_sub)

    # Covariates dataframe with donor factor
    covariates_sub <- data.frame(donor = donor_sub)

    # Run differential methylation test on subset
    myDiff <- calculateDiffMeth(
        meth_sub,
        covariates = covariates_sub,
        overdispersion = "MN",
        mc.cores = threads,
        save.db = FALSE
    )

    # Return the methylDiff object with significant results filtered
    # put very low qvalue threshold to get more results for plots
    sigDiff <- getMethylDiff(
        myDiff,
        difference = meth_difference,
        qvalue = qvalue,
        save.db = FALSE
    )

    return(sigDiff)
}

