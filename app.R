#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

data_input <- read_rds("preprocessed_data.RDS")

plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
}

bin <- function(x, brks) {
    k <- length(brks)
    res <- matrix(NA, nrow = k-1, ncol = 3)
    colnames(res) <- c("a", "b", "freq")
    res[,1] <- brks[-k]
    res[,2] <- brks[-1]
    temp <- .bincode(x, brks)
    temp <- temp[!is.na(temp)]
    
    for(i in 1:(k-1)) {
        res[i, 3] <- sum(temp == i)
    }
    res[res[, 3] != 0, ]
}

plot.mixfitEM <- function(x, ps = c("base", "ggplot2"), detail = TRUE, smoothness = 512, ...) {
    # browser()
    ps <- match.arg(ps)
    if(ps == "ggplot2") {
        requireNamespace("ggplot2")
    }
    family <- x$family
    #args <- mget(names(formals()),sys.frame(sys.nframe()))
    mc <- match.call()
    mc$ps <- NULL
    mc$detail <- detail
    mc$smoothness <- smoothness
    
    # fun.name <- ifelse(ps == "base", paste0("plot", family), paste0("ggplot", family))
    if(ps == "base") {
        fun.name <- paste0("plot", family)
        mc[[1]] <- as.name(fun.name)
    } else {
        fun.name <- paste0("ggplot", family)
        mc[[1]] <- as.name(fun.name)
        names(mc)[2] <- "object"
        ggplotlnorm(object = x, detail = mc$detail, smoothness = mc$smoothness)
    }
    # eval(mc, environment())
}

ggplotlnorm <- function(object, detail, smoothness, title = NULL, xlim, ylim, 
                        xlab, ylab, breaks, ..., theme = c("grey", "bw"), trans = 0.5) {
    # browser()
    pi <- object$pi
    mu <- object$mu
    sd <- object$sd
    mulog <- object$mulog
    sdlog <- object$sdlog
    data <- object$data
    family <- object$family
    ncomp <- length(pi)
    
    if(missing(xlab)) {
        xlab <- "Data"
    }
    if(missing(ylab)) {
        ylab <- "Density"
    }
    if(missing(breaks)) {
        breaks <- 30
    }
    if(missing(xlim)) {
        xlim <- c(max(min(mu - 3.5 * sd), 0), max(mu + 4 * sd))
    }
    
    # binwidth = (xlim[2] - xlim[1]) / breaks
    xseq <- seq(xlim[1], xlim[2], length = smoothness)	
    res <- matrix(NA, nrow = length(xseq), ncol = length(pi))
    for(i in 1:length(pi)) {
        res[ ,i] <-  pi[i] * dlnorm(xseq, mulog[i], sdlog[i])
    }
    yt <- apply(res, 1, sum)
    
    if(missing(ylim)) {
        if(is.matrix(data)) {
            count <- data[, 3]
            max_freq <- max(count) / (sum(count) * max(diff(data[, 1])))
            ylim <- c(0, max(c(yt, max_freq)))
        } else {
            brks <- seq(min(data), max(data), length = breaks + 1)
            tmp <- bin(data, brks = brks)
            count <- tmp[, 3]
            max_freq <- max(count) / (sum(count) * (brks[2] - brks[1]))
            ylim <- c(0, max(c(yt, max_freq)))
        }		
    }	
    if(is.matrix(data)) {
        breaks <- sort(unique(c(data[, 1], data[, 2])))
        data <- reinstate(data)	
    } else {
        breaks <- brks
    }
    
    # prepare data frames
    df1 <- data.frame(x = rep(xseq, ncomp), comp = rep(1:ncomp, each = smoothness), 
                      y = as.vector(res))
    df2 <- data.frame(x = xseq, y = yt)
    
    # plot
    if(is.null(title)) title = paste(family,"mixture density")
    if(detail) {
        add <- geom_polygon(data = df1, aes(x, y, fill = as.factor(comp)), alpha = trans)
    } else {
        add <- NULL
    }
    if(theme[1] == "bw") {
        theme <- theme_bw()
    } else {
        theme <- theme_grey()
    }
    ggplot(as.data.frame(data)) + 
        geom_histogram(aes(x = data, y = ..density..),breaks = breaks, color = "black", 
                       fill = "white", size = 0.3) + add + theme + 
        geom_path(data = df2, aes(x, y), ...) +
        scale_fill_discrete(guide = guide_legend(title = "Comp")) +
        labs(title = title, x = xlab, y = ylab) +
        theme(plot.title = element_text(hjust = 0.5))
}

names_vector <- data_input$full_data %>%
    select(-c(Component, Clusters, Order, PID, SalePrice)) %>%
    summarise_all(n_distinct) %>%
    pivot_longer(names_to = "field",
                 values_to = "value",
                 cols = everything()) %>%
    filter(value < 10) %>%
    select(-value) %>%
    pull %>%
    c("", .)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Value Analyser"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("clusters",
                        "Clusters:",
                        min = 2,
                        max = 9,
                        value = 3),
            HTML("This app is an example of a tool that can be used to 
            isolate the different sub-markets that exist within a larger, seemingly uniform market.
                 <br>
                 <br>
                 The dataset is currently house pricing data from AmesHousing"),
            width = 2),

        # Show a plot of the generated distribution
        mainPanel(
        fluidRow(
            column(width = 6, plotOutput("elbow_plot")),
            column(width = 6, plotOutput("hist_plot"))),
        fluidRow(
            column(width = 6, 
                   plotOutput("main_plot_1"),
                   fluidRow(selectizeInput("main_plot_1_rows",
                                  label = "Compare by rows of: ",
                                  choices = names_vector),
                   selectizeInput("main_plot_1_columns",
                                  label = "Compare by columns of: ",
                                  choices = names_vector))
                   ),
            column(width = 6, 
                   plotOutput("main_plot_2"),
                   fluidRow( selectizeInput("main_plot_2_rows",
                                  label = "Compare by rows of: ",
                                  choices = names_vector),
                   selectizeInput("main_plot_2_columns",
                                  label = "Compare by columns of: ",
                                  choices = names_vector))))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$elbow_plot <- renderPlot({
        data_input$test_fits_df %>%
            ggplot(aes(x = Clusters, y = BIC)) +
            geom_point() +
            geom_line() +
            geom_point(aes(x = Clusters, y = BIC, color = "red"), data = data_input$test_fits_df %>% filter(Clusters == 4), 
                       size = 4, shape = 1) +
            geom_point(aes(x = Clusters, y = BIC, color = "blue"), data = data_input$test_fits_df %>% filter(Clusters == input$clusters), 
                       size = 4, shape = 1) +
            scale_color_identity(name = "Cluster Selection: ",
                                 breaks = c("red", "blue"),
                                 labels = c("Optimal", "Selected"),
                                 guide = "legend") +
            theme(legend.position = "bottom") +
            labs(title = "Model Selector",
                 subtitle = "Lower y-axis implies a better fit")
    })
    output$hist_plot <- renderPlot({
        # browser()
        mix_model <- data_input$fits[[input$clusters]]
        
        plot <- plot(mix_model, ps = "ggplot2")
        
        plot <- plot + 
            labs(x = "SalePrice", title = "Histogram with clusters") +
            guides(fill = guide_legend(title = "Cluster"))
                   
        return(plot)    
    })
    output$main_plot_1 <- renderPlot({
        plot <- data_input$full_data %>%
            filter(Clusters == input$clusters) %>%
            arrange(SalePrice) %>%
            mutate(`Cumulative Shipments` = 1:n()) %>%
            ggplot(aes(x = `SalePrice`, y = `Cumulative Shipments`, color = Component)) +
            geom_point() + 
            guides(color = guide_legend(title = "Cluster")) +
            theme(axis.title.y  = element_blank())
        
        
        
        if (isTruthy(input$main_plot_1_rows) && isTruthy(input$main_plot_1_columns)) {
        plot <- plot + facet_grid(rows = vars(.data[[input$main_plot_1_rows]]),
                                  cols = vars(.data[[input$main_plot_1_columns]]))    
        } else if (isTruthy(input$main_plot_1_rows)) {
            plot <- plot + facet_grid(rows = vars(.data[[input$main_plot_1_rows]]))    
        } else if (isTruthy(input$main_plot_1_columns)) {
            plot <- plot + facet_grid(cols = vars(.data[[input$main_plot_1_columns]]))    
            
        }
        return(plot)
    })
    output$main_plot_2 <- renderPlot({
       plot <- data_input$full_data %>%
            filter(Clusters == input$clusters) %>%
            arrange(SalePrice) %>%
            mutate(`Cumulative Shipments` = 1:n()) %>%
            ggplot(aes(x = `SalePrice`, y = `Cumulative Shipments`, color = Component)) +
            geom_point() +
           guides(color = guide_legend(title = "Cluster")) +
           theme(axis.title.y = element_blank())
       
        
        if (isTruthy(input$main_plot_2_rows) && isTruthy(input$main_plot_2_columns)) {
            plot <- plot + facet_grid(rows = vars(.data[[input$main_plot_2_rows]]),
                                      cols = vars(.data[[input$main_plot_2_columns]]))    
        } else if (isTruthy(input$main_plot_2_rows)) {
            plot <- plot + facet_grid(rows = vars(.data[[input$main_plot_2_rows]]))    
        } else if (isTruthy(input$main_plot_2_columns)) {
            plot <- plot + facet_grid(cols = vars(.data[[input$main_plot_2_columns]]))    
        }
        return(plot)
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
