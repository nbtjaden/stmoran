library(spdep)
library(stmoran)
library(cluster)
library(spData)
library(dplyr)
library(ggplot2)
library(gridExtra)


#########################################################################
data(votes.repub)
head(votes.repub)

years <- names(votes.repub)

plot(us_states)
str(us_states)
us_states$NAME


## NA's currently still break moran.test_st (and localmoran_st??)
# hack: fill NA with plain 50% voters:
votes.repub[is.na(votes.repub)] <- 50


## build data set

votes <- cbind(NAME=row.names(votes.repub), votes.repub, mean=rowMeans(votes.repub))
# fix names
votes$NAME <- gsub("\\.", " ", votes$NAME)


## spatialize
votes <- inner_join(x = us_states, y = votes, by = "NAME")

str(votes)
plot(votes["X1976"])

##
votes_simple <- st_drop_geometry(votes[years]) %>% t() %>% data.frame()
##

###
neighbours <- spdep::poly2nb(votes, queen=TRUE)
weightslist <- nb2listw(neighbours=neighbours, style="W", zero.policy = TRUE)
###


##########
# test global I
I_1976 <- moran.test(x = votes$X1976, listw = weightslist)
I_1976


I_mean <- moran.test(x = votes$mean, listw = weightslist)
I_mean

I_st <- moran.test_st(x = votes_simple, listw = weightslist)
I_st



############################################
myblue1 <- "#c2d7e4"
myblue2 <- "#3f85ad"
myred1 <- "#ff977c"
myred2 <- "#ff4a2c"
mypurple1 <- "#c2a5cf"
mypurple2 <- "#7b3294"
mygreen1 <- "#a6dba0"
mygreen2 <- "#008837"
NAgrey <- "grey80"

clearmap <- theme_bw()+
theme(panel.grid = element_line(colour = "transparent"),
axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
axis.title.y=element_blank(),
axis.text.y=element_blank(),
axis.ticks.y=element_blank()
)
#########################################


##########
# test local I

Ii_mean <- localmoran(x = votes$mean, listw = weightslist)
Ii_st <- localmoran_st(x = votes_simple, listw = weightslist)

##########################################

Ii_plt_mean <- data.frame(Ii_mean)$Ii
Ii_plt_st <- data.frame(Ii_st)$Ii

plot_Ii <- function(x, title){
  ggplot() +
    geom_sf(data=votes,
            aes(fill = x),
            col="white")+
    scale_fill_viridis_c(name="Ii value",
                         limits = c(-0.1, 5))+
    ggtitle(title)+
    clearmap
}

grid.arrange(plot_Ii(Ii_plt_mean, title="1856-1976 mean data"),
             plot_Ii(Ii_plt_st, title="1856-1976 time series data"))

##########

siglevel=0.05

sig_mean <- Ii_mean[, ncol(Ii_mean)] <= siglevel
sig_st <- Ii_st[, ncol(Ii_st)] <= siglevel

plot_significance <- function(x, title){
  ggplot() +
    geom_sf(data=votes,
            aes(fill = x),
            col="white")+
    scale_fill_manual(values = c(mypurple1, mygreen1),
                      na.value = NAgrey,
                      name = "Ii significant")+
    ggtitle(title)+
    clearmap
}

grid.arrange(plot_significance(sig_mean, title="1856-1976 mean data"),
            plot_significance(sig_st, title="1856-1976 time seriesdata"))


##########
quadr_mean <- attr(Ii_mean, "quadr")
quadr_mean <- ordered(quadr_mean[,"mean"], levels=c("Low-Low", "Low-High", "High-High", "High-Low"))
quadr_st <- attr(Ii_st, "quadr")
quadr_st <- ordered(quadr_st[,"quadr"], levels=c("Low-Low", "Low-High", "High-High", "High-Low"))

plot_quadrants <- function(x, title){
  ggplot() +
    geom_sf(data=votes,
            aes(fill = x),
            col="white")+
    scale_fill_manual(values = c(myblue1, myblue2, myred1, myred2),
                      na.value = NAgrey,
                      name = "quadrants")+
    ggtitle(title)+
    clearmap
}

grid.arrange(plot_quadrants(quadr_mean, title="1856-1976 mean data"),
             plot_quadrants(quadr_st, title="1856-1976 time seriesdata"))


############################################

