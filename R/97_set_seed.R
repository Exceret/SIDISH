#' @keywords internal
set_seed <- function(seed, np, torch, random) {
    np$random$seed(seed)
    torch$manual_seed(seed)
    torch$cuda$manual_seed(seed)
    torch$cuda$manual_seed_all(seed) # if you are using multi-GPU.
    random$seed(seed)
    torch$backends$cudnn$deterministic <- TRUE
    torch$backends$cudnn$benchmark <- FALSE
}
