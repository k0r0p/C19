
setwd("~/Dropbox/Research/Corona/C19/")

source(file = "Model/C19.SEIR_Mod.R")

## Plotting functions

plt.db <-  function(a) {
  
  Yvar <- data.frame(c(seq(from = 0, to = 365)))
  for  (i in seq(length(IntL)))      {                                      
    
    
    purrr::map_dfc(IntL[i], ~purrr::map_dfc(.x, ~.x[, a]) ) %>% apply(  1, median) -> p.mat
    
    Yvar <- data.frame(Yvar, p.mat) 
  }
  
  colnames(Yvar) <- c("time", paste("Int", seq(1:15), sep = ''))
                      
  return(Yvar)
  
}


# Data frames for plotting functions

I.db <- tbl_df(plt.db(a = "I")) # Infectious
H.db <- tbl_df(plt.db(a = "H")) # Hospitalized
U.db <- tbl_df(plt.db(a = "U")) # ICU
D.db <- tbl_df(plt.db(a = "D")) # Dead

# setting color palette 

col.pal <- c("#36454f",
             "#1f78b4",
             "#ff6347",
             "#517B58",
             "#b25f8a",
             "#fb9a99",
             "#e31a1c",
             "#fdbf6f",
             "#ff7f00",
             "#cab2d6")



## Plots

# Theme plot

p.main <- function(db) {
  
  ggplot(db, aes(x = seq(from = 0, to = 365))) + 
    
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(colour = "maroon", size = 12, hjust = 0.5),
      plot.title.position = "panel",
      axis.title.x = element_text(colour = "maroon", size = 11),
      axis.text.x  = element_text(colour = "grey20", angle = 45, vjust = 0.7),
      axis.title.y = element_text(colour = "maroon", size = 11, angle = 90),
      axis.text.y  = element_text(colour = "grey20"),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      plot.tag  = element_text(colour = "grey20", size = 10)) +
    
    scale_x_continuous(breaks  = c(seq(from = 7, to = 337, by = 30)),
                       labels = c("April", "May", "June", "July", "August", "September",
                                  "October", "November", "December", "January",
                                  "February", "March")) 
  
}

#  Plots for deaths
# ====================

p2 <- p.main(D.db) + 
  
  geom_line(aes(y = Int1/max(Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int2/max(Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int3/max(Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int4/max(Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int5/max(Int1), color = col.pal[5]), size = 0.75) +
  
  
  scale_color_manual(name = "",
                     values = col.pal[1:5], 
                     breaks = col.pal[1:5],
                     labels = c("No Intervention", "Lockdown \n(Month 1)", "Lockdown \n(Months 1-2)", 
                                "Physical Dist. \n(Months 2-12))", "Physical Dist. \n(Months 2-6)") 
  ) +
  
  scale_y_continuous(breaks  = c(seq(from = 0, to = 1.2, by = 0.2))) +
  
  labs(tag = "(A)",
       x = "Month",
       y = "Annual mortality burden \n(as a proportion of base scenario)",
       #title = "Effects of interventions against COVID-19 on mortality burden",
       subtitle = "Lockdown and Physical Distancing")  +
  
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")


# print(p2)

## Activ case finding

p3 <- p.main(D.db) + 
  
  geom_line(aes(y = Int1/max(Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int6/max(Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int7/max(Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int8/max(Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int9/max(Int1), color = col.pal[5]), size = 0.75) +
  #geom_line(aes(y = Int10, color = col.pal[5]), size = 0.5) +
  
  
  scale_color_manual(name = "",
                     values = col.pal[1:5],
                     breaks = col.pal[1:5],
                     labels = c("Base Scenario", "ACF 5% \n(Months 2-12)", "ACF 5% \n(Months 2-6)", "ACF 5% \n(Months 3-12)",
                                "ACF 10% \n(Months 2-12)")   #, "ACF 10% \n(Months 2-6)"
  ) +
  
  scale_y_continuous(breaks  = c(seq(from = 0, to = 1.2, by = 0.2))) +
  
  labs(tag = "(B)",
       x = "Month",
       y = "Annual mortality burden \n(as a proportion of base scenario)",
       #title = "Effect of interventions on COVID-19 mortality burden",
       subtitle = "Active Case Finding (ACF) Isolation and Quarantine")  +
  
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")


# print(p3)


p4 <- p.main(D.db) + 
  
  geom_line(aes(y = Int1/max(Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int11/max(Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int12/max(Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int13/max(Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int14/max(Int1), color = col.pal[5]), size = 0.75) +
  
  
  
  scale_color_manual(name = "",
                     values = col.pal[1:5],
                     breaks = col.pal[1:5],
                     labels = c("Base \nScenario", "Lockdown + \nPhysical distancing + \nACF 5% (Months 2-12)", 
                                "Lockdown + \nPhysical distancing +\nACF 5% (Months 2-6)", 
                                "Lockdown + \nPhysical distancing +\nACF 5% (Months 3-12)",
                                "Lockdown + \nPhysical distancing +\nACF 10% (Months 2-12)")
  ) +
  
  # scale_y_continuous(breaks  = c(seq(from = 0, to = 500, by = 100))) +
  ylim(0, 0.05) +
  
  labs(tag = "(C)",
       x = "Month",
       y = "Annual mortality burden \n(as a proportion of base scenario)",
       #title = "Effect of interventions on COVID-19 mortality burden",
       subtitle = "Combination of Interventions")  +
  
  theme(legend.position = "bottom",
        #axis.text.x = element_blank(), 
        #axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")  + 
  annotate("text", x = 10, y = 0.04, label = "y-axis\ntruncated", size = 3)

# print(p4)

gridExtra::grid.arrange(p2, p3, p4, nrow = 3)


#==================================== 
## Hospital bed use plots
################################### 

ph2 <- 
  
  ggplot(H.db, aes(x = seq(from = 0, to = 365))) + 
  
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(colour = "maroon", size = 12, hjust = 0.5),
    plot.title.position = "panel",
    axis.title.x = element_text(colour = "maroon", size = 11),
    axis.text.x  = element_text(colour = "grey20", angle = 45, vjust = 0.7),
    axis.title.y = element_text(colour = "maroon", size = 11, angle = 90),
    axis.text.y  = element_text(colour = "grey20"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    plot.tag  = element_text(colour = "grey20", size = 10)) +
  
  scale_x_continuous(breaks  = c(seq(from = 7, to = 337, by = 30)),
                     labels = c("April", "May", "June", "July", "August", "September",
                                "October", "November", "December", "January",
                                "February", "March")) +
  
  geom_line(aes(y = Int1/max(H.db$Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int2/max(H.db$Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int3/max(H.db$Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int4/max(H.db$Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int5/max(H.db$Int1), color = col.pal[5]), size = 0.75) +
  
  geom_area(aes(y = Int1/max(H.db$Int1)), color = col.pal[1], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int2/max(H.db$Int1)), color = col.pal[2], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int3/max(H.db$Int1)), color = col.pal[3], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int4/max(H.db$Int1)), color = col.pal[4], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int5/max(H.db$Int1)), color = col.pal[5], size = 0.75, alpha = 0.2) +
  
  geom_hline(yintercept = 5400/max(H.db$Int1), lty = 2, col = "navy", lwd = 0.5) +
  
  annotate("text", x = 30, y = 7600/max(H.db$Int1), label = "HEALTH SYSTEM \nCAPACITY", size = 3) + 
  
  scale_color_manual(name = "",
                     values = col.pal[1:5], 
                     breaks = col.pal[1:5],
                     labels = c("No Intervention", "Lockdown \n(Month 1)", "Lockdown \n(Months 1-2)", 
                                "Physical Dist. \n(Months 2-12))", "Physical Dist. \n(Months 2-6)") 
  ) +
  
  scale_y_continuous(breaks  = c(seq(from = 0, to = 1.2, by = 0.2))) +
  
  labs(tag = "(A)",
       x = "Month",
       y = "General ward bed demand \n(as a proportion of base scenario)",
       #title = "Effects of interventions against COVID-19 on  health service demand",
       subtitle = "Lockdown and Physical Distancing")  +
  
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")


#print(ph2)

## Activ case finding

ph3 <- p.main(H.db) + 
  
  geom_line(aes(y = Int1/max(Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int6/max(Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int7/max(Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int8/max(Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int9/max(Int1), color = col.pal[5]), size = 0.75) +
  
  geom_area(aes(y = Int1/max(Int1)), color = col.pal[1], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int6/max(Int1)), color = col.pal[2], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int7/max(Int1)), color = col.pal[3], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int8/max(Int1)), color = col.pal[4], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int9/max(Int1)), color = col.pal[5], size = 0.75, alpha = 0.2) +
  
  geom_hline(yintercept = 5400/max(H.db$Int1), lty = 2, col = "navy", lwd = 0.5) +
  
  annotate("text", x = 30, y = 7600/max(H.db$Int1), label = "HEALTH SYSTEM \nCAPACITY", size = 3) + 
  
  
  scale_color_manual(name = "",
                     values = col.pal[1:5],
                     breaks = col.pal[1:5],
                     labels = c("Base Scenario", "ACF 5% \n(Months 2-12)", "ACF 5% \n(Months 2-6)", "ACF 5% \n(Months 3-12)",
                                "ACF 10% \n(Months 2-12)")   #, "ACF 10% \n(Months 2-6)"
  ) +
  
  scale_y_continuous(breaks  = c(seq(from = 0, to = 1.2, by = 0.2))) +
  
  labs(tag = "(B)",
       x = "Month",
       y = "General ward bed demand \n(as a proportion of base scenario)",
       #title = "Effect of interventions on COVID-19 mortality burden",
       subtitle = "Active Case Finding (ACF) Isolation and Quarantine")  +
  
  theme(legend.position = "bottom",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")


# print(ph3)

ph4 <- p.main(H.db) + 
  
  geom_line(aes(y = Int1/max(Int1), color = col.pal[1]), size = 0.75) +
  geom_line(aes(y = Int11/max(Int1), color = col.pal[2]), size = 0.75) +
  geom_line(aes(y = Int12/max(Int1), color = col.pal[3]), size = 0.75) +
  geom_line(aes(y = Int13/max(Int1), color = col.pal[4]), size = 0.75) +
  geom_line(aes(y = Int14/max(Int1), color = col.pal[5]), size = 0.75) +
  
  geom_area(aes(y = Int1/max(Int1)), color = col.pal[1], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int11/max(Int1)), color = col.pal[2], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int12/max(Int1)), color = col.pal[3], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int13/max(Int1)), color = col.pal[4], size = 0.75, alpha = 0.2) +
  geom_area(aes(y = Int14/max(Int1)), color = col.pal[5], size = 0.75, alpha = 0.2) +
  
  geom_hline(yintercept = 5400/max(H.db$Int1), lty = 2, col = "navy", lwd = 0.5) +
  
  annotate("text", x = 300, y = 4900/max(H.db$Int1), label = "HEALTH SYSTEM \nCAPACITY", size = 3) + 
  
  scale_color_manual(name = "",
                     values = col.pal[1:5],
                     breaks = col.pal[1:5],
                     labels = c("Base \nScenario", "Lockdown + \nPhysical distancing + \nACF 5% (Months 2-12)", 
                                "Lockdown + \nPhysical distancing +\nACF 5% (Months 2-6)", 
                                "Lockdown + \nPhysical distancing +\nACF 5% (Months 3-12)",
                                "Lockdown + \nPhysical distancing +\nACF 10% (Months 2-12)")
  ) +
  
  # scale_y_continuous(breaks  = c(seq(from = 0, to = 500, by = 100))) +
  ylim(0, 0.15) +
  
  labs(tag = "(C)",
       x = "Month",
       y = "General ward bed demand \n(as a proportion of base scenario)",
       #title = "Effect of interventions on COVID-19 mortality burden",
       subtitle = "Combination of Interventions")  +
  
  theme(legend.position = "bottom",
        #axis.text.x = element_blank(), 
        #axis.title.x = element_blank(),
        plot.tag = element_text(size = 9),
        plot.tag.position = "bottomleft")   +
  
  annotate("text", x = 0, y = 0.1, label = "y-axis\ntruncated", size = 3)

#print(ph4)

gridExtra::grid.arrange(p2, p3, p4, nrow = 3)

pdf(file = "Figure5.pdf", width = 8.5, height = 11)
gridExtra::grid.arrange(p2, p3, p4, nrow = 3)
#gridExtra::grid.arrange(ph2, ph3, ph4, nrow = 3)
dev.off()



#########################################################
##  FIGURE 3  Chart Plot
#########################################################


gt.df <- data.frame(inter.mat[,-7]) 
gt.df$Interventions <- rownames(gt.df)
rownames(gt.df) <- NULL

gt.df %>%
  reshape::melt(id.vars = "Interventions", 
                variable_name = "Inttype") %>%
  mutate(Intype = substr(Inttype, 1 ,2)) -> gt.dfn 

gt.dfn <- gt.dfn[,-2]

gt.dfn$Interventions <- factor(gt.dfn$Interventions, levels = unique(gt.dfn$Interventions))


print(p.gt)

p.gt <- 
  
  ggplot(gt.dfn, aes(x = value)) +
  geom_line(size = 6, aes(y = Interventions, color = Intype)) +
  
  
  theme(
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(colour = "maroon", size = 12, hjust = 0.5),
    plot.title.position = "panel",
    axis.title.x = element_text(colour = "maroon", size = 12),
    axis.text.x  = element_text(colour = "grey20", angle = 45, vjust = 0.7, size = 11),
    axis.title.y = element_text(colour = "maroon", size = 12, angle = 90),
    axis.text.y  = element_text(colour = "grey20", size = 11),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    
    plot.tag  = element_text(colour = "grey20", size = 10)) +
  
  scale_x_continuous(breaks  = c(seq(from = 7, to = 337, by = 30)),
                     labels = c("April", "May", "June", "July", "August", "September",
                                "October", "November", "December", "January",
                                "February", "March")) +
  
  
  geom_rect(
    fill = "slateblue", alpha = 0.005, 
    xmin = 31, xmax = 366,
    ymin = 10.5, ymax = 15.5) +
  
  
  scale_color_manual(name = "", 
                     values =  col.pal[5:7],
                     breaks = c("ld", "sd", "is"),
                     labels = c("Lockdown",
                                "Physical distancing",
                                "Active Case Finding"
                     )) +
  
  annotate("text", x = -20, y = 13, angle = 90,
           label = "COMBINATION\nINTERVENTIONS", size = 4, col = "slateblue") +
  
  
  labs(x = "Month",
       y = "Interventions")   +
  
  theme(legend.position = "top",
        legend.text = element_text(colour = "gray5", size = 10)) 

pdf(file = "Figure3_cap.pdf", width = 12, height = 8)
print(p.gt)
dev.off()



##########################################3
## FIGURE ONE CONTACT MATRIX PLOT
###########################################

pdf(file = "Figure1_tf.pdf", width = 8, height = 12)   

par(mfrow = c(3,1), oma = c(2,1,2,1))


par(mar = c(0,1,2,1))
plot(log(n.matT), 
     key = list(side = 3, cex = 1), 
     breaks = c(seq(from  = -5, to = 4, by = 9/6)),
     border=NA,
     
     axis.col = list(side = 2, cex = 0.9, pos = 0),
     axis.row = NULL,
     asp = 1,
     
     
     ylab = "", 
     xlab = "", 
     main = "",
     
     cex.axis = 0.8,
     cex.lab = 0.8,
     
     col.lab = "gray5",
     col.axis = "gray5"
)
mtext(side = 1, text = "Normal circumstances", col = "maroon", line = 0, cex = 1.1)
mtext(side = 3, text = "Daily mean contacts (log transformed)", 
      col = "maroon", line = 0, outer = TRUE, cex = 1)



par(mar = c(0,1,2,1))
plot(log(sd.mat), 
     key = NULL, 
     border=NA,
     
     axis.row = NULL,
     breaks = c(seq(from  = -5, to = 4, by = 9/6)),
     asp = 1,
     
     axis.col = list(side = 2, cex = 0.8, pos = 0),
     ylab = "", 
     xlab = "", 
     main = "",
     
     cex.axis = 0.8,
     cex.lab = 0.8,
     
     col.lab = "gray5",
     col.axis = "gray5"
)

mtext(side = 1, text = "Physical distancing", col = "maroon", line = 0, cex = 1.1)


par(mar = c(0,1,2,1))
plot(log(ld.mat), 
     key = NULL, 
     border=NA,
     
     axis.row = NULL,
     breaks = c(seq(from  = -5, to = 4, by = 9/6)),
     asp = 1,
     
     axis.col = list(side = 2, cex = 0.8, pos = 0),
     ylab = "", 
     xlab = "", 
     main = "",
     
     cex.axis = 0.8,
     cex.lab = 0.8,
     
     col.lab = "gray5",
     col.axis = "gray5"
     
)

mtext(side = 1, text = "Lockdown", col = "maroon", line = 0, cex = 1.1)


dev.off()
