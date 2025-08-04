


score_in_minmax <- function(reported_mean, n_items, min_response, max_response, scoring = c("sum", "mean")) {
  
  if(scoring == "sum"){
    
    min_score <- min_response * n_items
    max_score <- max_response * n_items
    
  } else if(scoring == "mean"){
    
    min_score <- min_response
    max_score <- max_response
    
  } else {
    stop("'scoring' must be either 'sum' or 'mean'")
  }
  
  result <- ifelse(min_score <= reported_mean & 
                     max_score >= reported_mean,
                   "consistent",
                   "inconsistent")
  
  res <- data.frame(result = result,
                    reported_mean = reported_mean,
                    min_score = min_score,
                    max_score = max_score,
                    n_items = n_items, 
                    min_response = min_response, 
                    max_response = max_response)
  
  return(res)  
}

score_in_minmax(reported_mean = 28.1, 
                n_items = 4, 
                min_response = 1, 
                max_response = 7,
                scoring = "sum")

score_in_minmax(reported_mean = 2.1, 
                n_items = 5, 
                min_response = -3, 
                max_response = +3,
                scoring = "mean")



