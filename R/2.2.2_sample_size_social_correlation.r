# 2.4 Sample size 

required_packages <- c("dplyr", "lubridate", "ggplot2",
                       "amt")

for (i in seq_along(required_packages)){
  p <- installed.packages()[,1]
  if(required_packages[i] %in% p){
    library(required_packages[i], character.only = T)
  }else{
    install.packages(required_packages[i])
    library(required_packages[i], character.only = T)
  } 
}

rm(i,p);gc()

# plotting theme
theme_plt <-  ggplot2::theme(panel.background = element_rect(fill = "white",
                                                             inherit.blank = TRUE),
                             panel.grid.major = element_line(linetype = "solid",
                                                             linewidth = 0.5,
                                                             colour = "#efe9e8",
                                                             inherit.blank = FALSE),
                             plot.background = element_rect(color = "black",
                                                            size = 1),
                             axis.text.x = element_text(size = 20, angle = 0),
                             axis.text.y = element_text(size = 20), 
                             axis.line = element_line(colour = "black"),
                             axis.title.x = element_text(size = 22),
                             axis.title.y = element_text(size = 22),
                             plot.margin = ggplot2::margin(t = 1,  # Top margin
                                                           r = 2,  # Right margin
                                                           b = 1,  # Bottom margin
                                                           l = 2,  # Left margin
                                                           unit = "cm"),
                             title = element_text(size = 20),
                             legend.text = element_text(size = 22),
                             legend.title = element_text(size = 24),
                             strip.text = element_text(size = 20))


## Evaluation -----

# get random draws to process
set.seed(1000)
sample_size <- seq(2, 88, 1)

for(n in seq_along(sample_size)) {
  id_list <- list()
  for(i in 1:100){
    id_list[[i]] <- sample(final_ids, # vector collared animals
                           size = sample_size[[n]],
                           replace = FALSE)
  }
  out <- do.call(rbind, id_list)
  saveRDS(out,
          file = paste0(".sample_size/ids_",
                        sample_size[[n]],
                        ".rds"))
}
# ***************************************************************
# ***************************************************************

# ***************************************************************
# extraction to be run on the workstation
# run the first 9 on here
# ***************************************************************
# ***************************************************************

sample_sizes <- seq(2, 88, 1)
sample_files <- paste0("./sample_size/ids_", sample_sizes, ".rds")

ids_list <- lapply(sample_files,
                   readRDS)

# fix_rates <- c(0.5, 1, 2, 4, 12, 24, 168)
dist_scale <- seq(20, 300, 20)

# leave 2 cores for other system processes
library(doMC)
library(doParallel)
cores_to_use <- detectCores() - 2
registerDoMC(cores_to_use)

summary_data <- foreach(u = seq_along(ids_list)) %dopar% {# sample size loop
  result_df <- data.frame()
  for(k in 1:dim(ids_list[[u]])[1]) { # sampling variation loop
    
    tmp_id <- ids_list[[u]][k,]
    
    df_list <- list()
    for(r in seq_along(tmp_id)){
      df_list[[r]] <- df %>%
        filter(ID == tmp_id[r])
    }
    
    out <- data.frame()
    for(i in 1:(length(df_list) - 1)){
      start_ind <- i + 1
      for(j in start_ind:length(df_list)){
        tmp_out <- dyad_dataframe(
          data_1 = df_list[[i]],
          data_2 = df_list[[j]])
        
        out <- rbind(out, tmp_out)
      }
    }
    
    if(nrow(out)==0){
      next
    }else{
      out <- out %>%
        tidyr::drop_na(distance) %>%
        group_by(dyad_id) %>%
        arrange(t) %>%
        mutate(index = seq(1, n(), 1)) %>%
        ungroup()
      
      out <- out %>%
        group_by(dyad_id) %>%
        nest_by() %>%
        ungroup() %>%
        mutate(day_2 = purrr::map(data, function(x) infectious_period_block(x,
                                                                            period = "day",
                                                                            window_size = 2,
                                                                            out_type = "data_frame")),
               day_7 = purrr::map(data, function(x) infectious_period_block(x,
                                                                            period = "day",
                                                                            window_size = 7)),
               day_14 = purrr::map(data, function(x) infectious_period_block(x,
                                                                             period = "day",
                                                                             window_size = 14)),
               day_30 = purrr::map(data, function(x) infectious_period_block(x,
                                                                             period = "day",
                                                                             window_size = 30)),
               day_183 = purrr::map(data, function(x) infectious_period_block(x,
                                                                              period = "day",
                                                                              window_size = 183))) %>%
        select(-data) %>%
        tidyr::unnest(c(day_2, day_7, day_14, day_30, day_183)) %>%
        data.table::data.table() %>%
        rename(day_2 = group)
      
      # get the different def of distance and calculate variance in togetherness 
      contact_out <- do.call(rbind, 
                             lapply(dist_scale,
                                    function(x) contact_only_def_ip(x = out,
                                                                    dist_scale = x))) 
      
      # summarize
      summary_df <- contact_out %>%
        dplyr::group_by(dyad_id, infect_period, contact_dist) %>%
        summarise(N = sum(N),
                  contact_sum = sum(num_contact),
                  p_contact = contact_sum/N,
                  mean = (contact_sum * (1 - p_contact))/p_contact,
                  variance = (contact_sum*(1 - p_contact))/p_contact^2,
                  cv = sqrt(variance)/mean) %>%
        ungroup() %>%
        suppressMessages()
      
      # bind the result
      
      result_df <- rbind(result_df, summary_df)
    }
  }
  
  result_df <- result_df %>%
    mutate(sample_size = sample_sizes[[u]])
  
  saveRDS(result_df,
          file = paste0("./sample_size_",
                        sample_sizes[[u]],
                        ".rds"))
  return(result_df)
}

files_summary <- list.files(path = "./sample_size/",
                            pattern = "sample_size_",
                            full.names = TRUE)


summary_data_rc <- lapply(files_summary, readRDS)

# function to extract metric at sample size level
summary_function <- function(x) {
  out <- x %>%
    group_by(sample_size, infect_period, interaction_dist) %>%
    summarise(mean_mean = mean(mean),
              sd_mean = sd(mean),
              cv_mean = sd_mean/mean_mean,
              mean_cv = mean(cv, na.rm = TRUE),
              sd_cv = sd(cv, na.rm = TRUE),
              mean_var = mean(variance)) %>%
    suppressMessages()
  
  return(out)
}

final_out_rc <- do.call(rbind, lapply(
  summary_data_rc, summary_function
)) %>%
  mutate(sample_size = as.factor(sample_size),
         contact_dist = as.factor(interaction_dist))

plt_data_rc <- final_out_rc %>%
  filter(contact_dist == 40) %>%
  mutate(infect_period = case_when(
    infect_period == "2day" ~ "2 days",
    infect_period == "7day" ~ "7 days",
    infect_period == "14day" ~ "14 days",
    infect_period == "30day" ~ "30 days",
    infect_period == "183day" ~ "183 days",
  ),
  infect_period = factor(infect_period,
                         levels = c("2 days", "7 days", "14 days",
                                    "30 days", "183 days")))

# mean of the cv 
ggplot(plt_data_rc) +
  geom_point(aes(x = sample_size,
                 y = mean_cv),
             size = 2,
             alpha = 0.5) +
  facet_wrap(~infect_period, 
             axis.labels = "margins") +
  scale_x_discrete(breaks = seq(0, 45, 10)) +
  theme_plt +
  labs(x = "Sample size", 
       y = "Mean cv of togetherness") +
  theme(legend.position = "bottom")