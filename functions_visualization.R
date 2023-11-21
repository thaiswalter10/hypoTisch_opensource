draw_rscpline <- function(var, xlab_sr, outc = postop_ischemia, title = "Probability of postop ischemia as a binary variable", df=analyze_df){
  rcspline.plot(
    x = df %>% pull({{var}}), 
                y=df %>% pull({{outc}}), 
                model="logistic", 
                nk=3,
                show="prob",
                plotcl=TRUE, 
                showknots=FALSE, 
                add=FALSE,
                lty=1, 
                noprint=TRUE,
                ylim=c(0,1), 
                xlab=xlab_sr, 
                main=title,
                statloc = "none")  
}

draw_tau_per_bootstrap <- function(tau_per_bootstrap){
  tau_per_bootstrap %>% 
    select(-splits) %>% 
    unnest(cols = c(models)) %>%
    ggplot(aes(x=estimate))+
    geom_histogram(bins = 20, col = "white")+
    facet_wrap(~variable, scales = "free_x")
}