#' Memory Cleanup Function
#'
#' This function performs memory cleanup by running garbage collection and removing unnecessary objects.
#' It is designed to be used as an on.exit handler to ensure memory is properly managed.
#' 
#' @details The function utilizes the on.exit() function in R to automatically execute the cleanup code
#'          when the function scope is exited.
#'
#' @return This function does not return any value.
#'
#' @examples
#' memory_cleanup()  # Call this function to clean up memory
#'
#' @export
#'
memory_cleanup <- function() {
  on.exit({
    gc()  # Clean memory
    rm(list = ls())  # Remove unnecessary objects
  })
}
