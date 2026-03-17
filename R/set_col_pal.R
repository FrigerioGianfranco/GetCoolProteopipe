#' Set the col_pal option.
#'
#' It sets the option GetCoolProteopipe.col_pal, so it will be consistently applied to all the functions that have the col_pal argument.
#'
#' @param col_pal NULL or a character vector with colors. The option GetCoolProteopipe.col_pal will be set to it.
#'
#' @return Nothing. It only sets the GetCoolProteopipe.col_pal option.
#'
#'
#' @examples
#' \dontrun{
#'
#' set_col_pal(c(S = "green", V = "orange"))
#'
#' }
#'
#'
#'
#' @export
set_col_pal <- function(col_pal = NULL) {

  are.colors <- function (vect) {
    map_lgl(vect, ~tryCatch({
      is.matrix(col2rgb(.)) & ncol(col2rgb(.))>0
    }, error = function(e) {
      FALSE
    }))
  }

  if (!is.null(col_pal)) {
    if (is.null(names(col_pal))) {
      cat(paste0('\n --- the col_pal option is now set to be  c("', paste0(col_pal, collapse = '", "'), '") ---\n\n'))
    } else {
      cat('\n --- the col_pal option is now set to be  c(')
      for (i in seq(length(col_pal))) {
        if (names(col_pal)[i]!="") {
          cat(paste0(names(col_pal)[i], ' = "'))
        } else {
          cat('"')
        }
        cat(paste0(col_pal[i], '"'))
        if (i != length(col_pal)) {
          cat(', ')
        }
      }
      cat(') ---\n\n')
    }

  } else {
    cat(paste0('\n --- the col_pal option is now set to be NULL ---\n\n'))
  }

  options(GetCoolProteopipe.col_pal = col_pal)
  invisible(getOption("GetCoolProteopipe.col_pal"))
}
