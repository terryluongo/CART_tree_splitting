
  test <- make_random_tibble(range = 30,
                             distance = 0.05,
                             noise = 2,
                             sections = 3,
                             section_ratio = 0.8)
  
  # super weird seems to like cutting more as x increases, has to be a reason
  # the max_cuts really impacts the quality of the tree
  # seems like pruning is where the algorithm works best
  
  make_splits(test,
              max_cuts = 20,
              alpha = 0.5, 
              radius = 2, 
              iterations = 10,
              threshold = 0.01)
  
  
  make_splits <- function(df, max_cuts, alpha, radius, iterations, threshold,
                          round_x_iterations = 5) {
    splits <- c()
    for (i in 1:max_cuts) {
      new <- find_global_split(df,alpha,radius,iterations,splits,round_x_iterations)
      splits <- c(splits,new)
    }
    
    splits <- prune_splits(df,splits,threshold)
    
    plot_parts(df,splits)
  }
  
  
  # recursively prune worst split
  prune_splits <- function(df, splits, threshold) {
    init_rmse <- calculate_rmse(df,splits,NULL)
    rmse_without_each <- c()
    for (i in 1:length(splits)) {
      without_this <- calculate_rmse(df,splits[splits != splits[i]],NULL)
      rmse_without_each <- c(rmse_without_each, without_this - init_rmse)
    }
    sorted_indices <- order(rmse_without_each)
    worst_split <- splits[sorted_indices[1]]
    worst_rmse <- rmse_without_each[sorted_indices[1]]
    
    ratio <- worst_rmse / init_rmse
    print(ratio)
    if (length(splits) < 3 | ratio > threshold) {
      return(splits)
    }
    else {
      prune_splits(df,splits[splits != worst_split],threshold)
    }
    
  }
  
  # calculates midpoints given all previous splits, returns best split from ones starting at each midpoint
  find_global_split <- function(df,alpha,radius,iterations,
                                prev_splits = NULL,
                                round_x_iterations = 5,
                                should_plot = FALSE) {
    splits <- sort(prev_splits)
    splits <- c(min(df$x)-1,splits,max(df$x)+1)
    midpoints <- calculate_midpoints(splits)
    minima <- tibble(split = numeric(), rmse = numeric())  
    print(splits)
    for(midpoint in midpoints) {
      local_min <- step_gradient_descent(df,midpoint,alpha,radius,iterations,
                                         prev_splits,
                                         round_x_iterations,should_plot)
      minima <- bind_rows(minima, local_min)
      
    }
    
    minima %>% arrange(rmse) %>% head(1) %>% select(split) %>% pull()
  
  }
  
  
  # finds local minima, should return in a tibble with RMSE and break point
  # either stops when 0 move reached or after iterations has elapsed rounds last
  # round_x_iterations -> if it is cycling around best point
  step_gradient_descent <- function(df,init,alpha,radius,iterations,
                                    prev_splits = NULL,
                                    round_x_iterations=5,
                                    should_plot = FALSE) {
    zeroed <- FALSE
    iterations_to_round <- c()
    points_visited_x <- c()
    points_visited_y <- c()
    index <- which.min(abs(df$x - init))[[1]]
    
    for (i in 1:iterations) {
      
      value <- df$x[index]
      
      points_visited_x <- c(points_visited_x, value)
      points_visited_y <- c(points_visited_y, df$y[index])
      
      delta <- find_avg_delta(df,radius,prev_splits, value)
      move <- round(delta * alpha)   
      index <- index + move 
      
      if (move == 0) {
        zeroed <- TRUE
        break
      }
      
      if (iterations - i < round_x_iterations) {
        iterations_to_round <- c(iterations_to_round,df$x[index])
      }
    }
    
    if (should_plot) {
      plot(df$x,df$y)
      points(points_visited_x,points_visited_y,col="red",bg="red",pch=16)
    }
    
    split <- ifelse(zeroed,df$x[index],mean(iterations_to_round))
    rmse <- calculate_rmse(df,prev_splits,split)
    minima <- tibble(split = split, rmse = rmse)
    
    return(minima)
  }
  
  
  # calculates average delta around radius of point, by index
  find_avg_delta <- function(df, radius, prev_splits, point) {
    index <- which.min(abs(test$x - point))[[1]]
    
    rmse_vector <- sapply((index - radius):(index + radius), function(i) calculate_rmse(df, prev_splits, df$x[i]))
    -mean(diff(rmse_vector))
  }
  
  
  # calculate rmse for given splits 
  calculate_rmse <- function(df,prev_splits,current) {
    splits <- prev_splits
    if (!is.null(current) && !(current %in% prev_splits)) {
      splits <- c(prev_splits,current)
    }
    splits <- c(min(df$x)-1,splits,max(df$x)+1)
    df <- df %>% mutate(group = cut(x, breaks = splits)) %>%
      group_by(group) %>% 
      mutate(group_mean = mean(y)) %>% 
      mutate(component_rmse = (y-group_mean)^2) %>% 
      ungroup() %>% 
      summarise(total_rmse = sum(component_rmse))
    df %>% select(total_rmse) %>% pull()
  }
  
  # given array of numbers, calculates midpoints returning n-1 numbers
  calculate_midpoints <- function(array) {
    midpoints <- numeric(0)
    # Iterate over pairs of consecutive numbers in the array
    for (i in seq(1, length(array)-1)) {
      start <- array[i]
      end <- array[i + 1]
      # Calculate the midpoint and append it to the result
      midpoint <- (start + end) / 2
      midpoints <- c(midpoints, midpoint)
    }
    
    return(midpoints)
  }
  
  # for n splits in dataframe single feature, calculates group response means for n+1 sides 
  calculate_means <- function(df,splits) {
    splits <- c(min(df$x)-1,splits,max(df$x)+1)
    df <- df %>% mutate(group = cut(x, breaks = splits)) %>%
      group_by(group) %>% 
      summarise(group_mean = mean(y))
    df
  }
  
  calculate_means(test,c(-2,0,2))
  
  # generate synthetic data in tibble
  
  make_random_tibble <- function(range, distance, noise, sections, section_ratio) {
    x <- seq(from = -range, to = range, by = distance)
    y <- runif(length(x)) * noise
    
    chunk_length <- length(y) / sections
    
    for (i in 1:sections) {
      start_index <- ((i - 1) * chunk_length + 1)
      end_index <- (i * chunk_length)
      y[start_index:end_index+1] <- y[start_index:end_index+1] + runif(1) * noise * section_ratio
    }
    tibble(x,y)
  }
  
  # plot what splits look like for dataset 
  plot_parts <- function(df,splits) {
    vlines <- data.frame(xintercept = splits, index = seq_along(splits))
    
    df %>% ggplot(mapping = aes(x = x, y = y)) +
      geom_point() +
      geom_vline(data = vlines, mapping = aes(xintercept = xintercept)) +
      geom_text(data = vlines, aes(x = xintercept, label = index, y = 3),
                vjust = -0.5, hjust = 0.5, color = "red")  # Adjust vjust and hjust as needed
  }
