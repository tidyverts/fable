# TODO: Write quantile.fcdist

#' Create a forecast distribution object
#'  
#' @param f A distribution function producing quantiles (such as `qnorm`)
#' @param ... Arguments for `f` function
#' @param transformation Transformation to be applied to resulting quantiles from `f`
#' @param abbr Abbreviation for display purposes, defaults to the object name of `f`
#' 
#' @examples 
#' mydist <- new_fcdist(qnorm, mean = rep(3, 10), sd = seq(0, 1, length.out=10),
#'  transformation = exp, abbr = "N")
#' mydist
#' hilo(mydist, 95)
#' @export
new_fcdist <- function(f, ..., transformation = ~ .x, abbr = NULL){
  f_quo <- enquo(f)
  t_fn <- as_mapper(transformation)
  pmap(dots_list(...), list) %>%
    enclass("fcdist",
            f = f,
            t = t_fn,
            qname = abbr%||%quo_text(f_quo),
            trans = !is.name(body(t_fn)))
}

#' @export
type_sum.fcdist <- function(x){
  "dist"
}

#' @export
obj_sum.fcdist <- function(x) {
  rep("dist", length(x))
}

#' @export
pillar_shaft.fcdist <- function(x, ...){
  pillar::new_pillar_shaft_simple(format(x), align = "left", min_width = 10)
}

#' @export
print.fcdist <- function(x, ...) {
  print(format(x, ...), quote = FALSE)
  invisible(x)
}

#' @export
format.fcdist <- function(x, ...){
  x %>%
    map_chr(function(qt){
      args <- qt %>%
        imap(~ paste0(ifelse(nchar(.y)>0, paste0(.y, " = "), ""),
                      format(.x, trim = TRUE, digits = 2))) %>%
        invoke("paste", ., sep = ", ")
      out <- paste0(
        attr(x, "qname"),
        "(", args, ")"
      )
      if(attr(x, "trans")){
        paste0("t(", out, ")")
      }
      else{
        out
      }
    })
}

#' @export
`[.fcdist` <- function(x, ...){
  enclass(NextMethod(), "fcdist", 
          !!!attributes(x))
}

#' @export
c.fcdist <- function(...){
  sameAttr <- dots_list(...) %>%
    map(~ if(!inherits(.x, "fcdist")) {abort("Only combining fcdist objects is supported")} else {attributes(.x)}) %>%
    duplicated %>%
    .[-1]
  if(any(!sameAttr)){
    abort("Cannot combine fcdist objects of different types.")
  }
    
  enclass(NextMethod(), "fcdist", 
          !!!attributes(..1))
}

#' @export
length.fcdist <- function(x){
  NextMethod()
}

# #' @importFrom ggplot2 aes_
# autoplot.fcdist <- function(q_fn, q_range = c(0.0001, 0.9999), precision = 0.01){
#   tibble(x = seq(q_range[1], q_range[2], by = precision)) %>%
#     mutate(!!"density":=q_fn(!!sym("x"))) %>%
#     ggplot(aes_(x=~x, y=~density)) + 
#     geom_line()
# }
#' @export
hilo.fcdist <- function(x, level = 95, ...){
  if(length(level)!=1){
    abort("Only one value of 'level' is supported.")
  }
  if (level < 0 || level > 100) {
    abort("'level' can't be negative or greater than 100.")
  }
  list(lower = 50-level/2, upper = 50+level/2) %>%
    map(function(level){
      eval_tidy(quo(attr(x, "t")(attr(x, "f")(level/100, !!!merge_pos_list(!!!x))))) %>%
        unlist(recursive = FALSE, use.names = FALSE)
    }) %>%
    append(list(level = level)) %>%
    invoke("new_hilo", .)
}

#' @importFrom stats quantile
#' @export
quantile.fcdist <- function(x, probs = seq(0, 1, 0.25), ...){
  map(probs, function(prob){
    eval_tidy(quo(attr(x, "t")(attr(x, "f")(prob, !!!merge_pos_list(!!!x))))) %>%
      unlist(recursive = FALSE, use.names = FALSE)
  })
}