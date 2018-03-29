# Lookup table for function inverses
#' @importFrom dplyr first
#' @importFrom purrr compact possibly
inverse_table <- function() {
  table <- new.env(parent = emptyenv())
  list(
    add = function(ns, fn, inv) {
      table[[ns]] <- as.list(table[[ns]])
      table[[ns]][[fn]] <- inv
    },
    
    get = function(ns, fn) {
      ret <- table[[ns]][[fn]]
      if (is.null(ret)) {
        stop("No supported inverse for this function")
      }
      ret
    })
}

#' @importFrom methods findFunction
undo_transformation <- function(operation, target, result){
  fn <- call_name(operation)
  env <- get_env(operation, caller_env())
  ns <- eval_tidy(expr(findFunction(!!fn)), env = env) %>%
    map(possibly(environmentName, otherwise = NULL)) %>% 
    compact %>%
    first
  
  inverse_table$get(ns, fn)(operation, get_expr(target), result)
}

#' @importFrom dplyr last
invert_transformation <- function(transformation){
  transformation <- enquo(transformation)
  
  # Evaluate transformation stack via call traversal along longest argument (hopefully, the data)
  
  transformation_stack <- traverse_call(!!transformation,
                                        f = ~ append(.y, .x[[1]]),
                                        g = ~ .x %>%
                                          get_expr %>%
                                          as.list %>% 
                                          map(new_quosure, env = get_env(.x)) %>%
                                          .[which.max(map(., ~ length(eval_tidy(.x))))])
  
  # Iteratively undo transformation stack
  result <- expr(!!sym("x")) #last(transformation_stack)
  for (i in seq_len(length(transformation_stack) - 1)){
    result <- undo_transformation(transformation_stack[[i]], transformation_stack[[i+1]], result)
  }
  new_function(eval_tidy(expr(alist(x = !!get_expr(last(transformation_stack))))), result)
}

inverse_table <- inverse_table()

map(c("log", "logb"),
    ~ inverse_table$add("base", .x, 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    if(length(args) == 1){
                      call2("exp", !!!replace(args, target_pos, list(result)))
                    }
                    else{
                      call2("^", !!!exprs(!!args[[-target_pos]], !!result))
                    }
                  }
    )
)

inverse_table$add("base", "log10",
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("^", !!!exprs(10, !!result))
                  }
)

inverse_table$add("base", "log2",
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("^", !!!exprs(2, !!result))
                  }
)

inverse_table$add("base", "log1p",
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("expm1", !!!exprs(!!result))
                  }
)

inverse_table$add("base", "exp", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("log", !!!replace(args, target_pos, list(result)))
                  }
)

inverse_table$add("package:forecast", "BoxCox", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("InvBoxCox", !!!replace(args, target_pos, list(result)))
                  }
)

inverse_table$add("package:forecast", "InvBoxCox", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("BoxCox", !!!replace(args, target_pos, list(result)))
                  }
)

inverse_table$add("base", "sqrt", 
                  function(operation, target, result){
                    call2("^", !!!exprs(!!result, 2))
                  }
)

inverse_table$add("base", "^", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("^", !!!exprs(!!result, !!call2("/", !!!exprs(1, !!args[[2]]))))
                  }
)

inverse_table$add("base", "+", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    if(length(args) == 1){
                      call2("+", !!!exprs(!!result))
                    }
                    else{
                      call2("-", !!!exprs(!!result, !!args[[-target_pos]]))
                    }
                  }
)

inverse_table$add("base", "-", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    if(length(args) == 1){
                      call2("-", !!!exprs(!!result))
                    }
                    else{
                      if(target_pos == 1){
                        call2("+", !!!exprs(!!result, !!args[[2]]))
                      }
                      else{
                        call2("-", !!!exprs(!!args[[1]], !!result))
                      }
                    }
                  }
)

inverse_table$add("base", "/", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    if(target_pos == 1){
                      call2("*", !!!exprs(!!args[[2]], !!result))
                    }
                    else{
                      call2("/", !!!exprs(!!args[[1]], !!result))
                    }
                  }
)

inverse_table$add("base", "*", 
                  function(operation, target, result){
                    args <- call_args(operation)
                    target_pos <- match(list(target), args)
                    call2("/", !!!exprs(!!result, !!args[[-target_pos]]))
                  }
)

inverse_table$add("base", "(", 
                  function(operation, target, result){
                    call2("(", !!!exprs(!!result))
                  }
)