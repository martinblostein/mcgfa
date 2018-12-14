create_progress_bar = function(rG,num_G,num_q,num_model,models) {
    CCC_equiv <- c("CCC", "CUC", "UCC", "UUC")
    CCU_equiv <- c("CCU", "CUU", "UCU", "UUU")

    # determine length
    if (1 %in% rG) {
        if (any(models %in% CCC_equiv) && any(models %in% CCU_equiv)) {
            max <- num_model * num_q * (num_G - 1) + 2
        } else {
            max <- num_model * num_q * (num_G - 1) + 1
        }
    } else {
        max <- num_model*num_q*num_G
    }

    # create bar object & return
    txtProgressBar(min = 0, max = max, style = 3)
}
