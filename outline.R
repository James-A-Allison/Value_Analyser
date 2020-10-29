library(mixR)
library(tidyverse)
library(ggnewscale)

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

ames_data <- AmesHousing::ames_raw

ggplot(ames_data, aes(x = SalePrice)) + geom_histogram(bins = 50)
ggplot(ames_data, aes(x = SalePrice, fill = `MS SubClass`)) + geom_histogram(bins = 50)
ggplot(ames_data, aes(x = SalePrice, fill = `MS Zoning`)) + geom_histogram(bins = 50)

ggplot(ames_data, aes(x = SalePrice)) + 
  geom_histogram(bins = 50) +
  stat_function(geom = "line",
                fun = plot_mix_comps,
                args = list(mix_model$mu[1], 
                            mix_model$sigma[1],
                            mix_model$lambda[1]),
                colour = "red", lwd = 1.5)


test_fits <- mixR::select(ames_data$SalePrice, 2:9, family = "lnorm")
plot(test_fits)

fit <- mixfit(ames_data$SalePrice, ncomp = 2, family = "lnorm")

ggplot(ames_data, aes(x = SalePrice)) + 
  geom_histogram(bins = 50) +
  stat_function(geom = "line",
                fun = plot_mix_comps,
                args = list(fit$mu[1], 
                            fit$sigma[1],
                            fit$lambda[1]),
                colour = "red", lwd = 1.5)

mix_fn <- function(ncomp, data = ames_data$SalePrice, family = "lnorm"){
  mixfit(ncomp = ncomp, x = data, family = family)
}

fits <-  map(.x = 2:9, ~mixfit(x = ames_data$SalePrice, ncomp = .x, family = "lnorm"))

fits <- list()
for (i in 2:9) {
  fits[[i]] <- mixfit(ames_data$SalePrice, ncomp = i, family = "lnorm")
}

full_data <- lapply(fits[2:9], "[[",  11) %>%
  map(as_tibble) %>%
  map(.f = ~mutate(.data = .x, ID = 1:n())) %>%
  imap(.f = ~mutate(.data = .x, Clusters = .y + 1)) %>%
  map(.f = ~pivot_longer(data = .x,
                         cols = -c(Clusters, ID),
                         names_to = "Component",
                         values_to = "Value")) %>%
  map(.f = ~group_by(.data = .x, ID, Clusters)) %>%
  map(.f = ~slice_max(.data = .x, order_by = Value, n = 1)) %>%
  map(.f = ~ungroup(x = .x)) %>%
  map(.f = ~dplyr::select(.data = .x, Component, Clusters)) %>%
  map(.f = ~bind_cols(.x, ames_data)) %>% 
  bind_rows() %>%
  mutate(Component = gsub("V", "", Component))
  

  
fitted_values <- fit$comp.prob %>% 
  as_tibble() %>%
  mutate(ID = 1:n()) %>%
  pivot_longer(names_to = "Component",
               values_to = "Value",
               -ID) %>%
  group_by(ID) %>%
  slice_max(order_by = Value, n = 1) %>%
  ungroup %>%
  select(Component) %>%
  bind_cols(ames_data)

ggplot(fitted_values, aes(x = SalePrice, fill = Component, ..density..)) + 
  geom_histogram() +
  facet_wrap(facets = vars(`MS Zoning`))

full_data %>%
  filter(Clusters == 4) %>%
  arrange(SalePrice) %>%
  mutate(`Cumulative Revenue` = cumsum(SalePrice),
         `Cumulative Shipments` = 1:n()) %>%
  ggplot(aes(x = `SalePrice`, y = `Cumulative Shipments`)) +
  geom_point(aes(color = `MS Zoning`)) +
  new_scale_color() +
  stat_ellipse(aes(color = Component)) +
  scale_color_brewer(palette = "Blues") +
  facet_wrap(facets = vars(`Yr Sold`))

test_fits_df <- tibble(Clusters = test_fits$ncomp,
                       BIC = test_fits$bic)

test_fits_df %>%
  ggplot(aes(x = Clusters, y = BIC)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = Clusters, y = BIC, color = "red"), data = test_fits_df %>% filter(Clusters == 4), 
             size = 4, shape = 1) +
  geom_point(aes(x = Clusters, y = BIC, color = "black"), data = test_fits_df %>% filter(Clusters == 5), 
             size = 4, shape = 1) +
  scale_color_identity(name = "Cluster Selection: ",
                       breaks = c("red", "black"),
                       labels = c("Optimal", "Selected"),
                       guide = "legend") +
  theme(legend.position = "bottom")

plot(fit, ps = "ggplot2", xlab = "Price", title = "Histogram with clusters") +
  guides(fill = guide_legend(title = "Cluster"))

output <- list()
output$full_data <- full_data
output$test_fits_df <- test_fits_df
output$fits <- fits

saveRDS(output, "preprocessed_data.RDS")
