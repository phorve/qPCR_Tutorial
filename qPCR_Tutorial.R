qp1 <- read.csv('/path/to/file.csv', header = T,skip = 7, stringsAsFactors = F)

qp1[qp1 == "Undetermined"] <- NA

qp1.1 <- qp1 %>%
  dplyr::group_by(Sample.Name) %>% #, run) %>%
  dplyr::summarise(Ct = mean(as.numeric(as.character(Ct)))) %>%
  as.data.frame()

# Make an object with just our standards
standards1 <- qp1.1[grepl('Standard', qp1.1$Sample), ]

# Remove the standards from the data set with all the other samples    
qp1.1 <- qp1.1[!grepl('Standard', qp1.1$Sample), ]

# Describe our dilution curve and how many gene copies are in each       standard
dil <- c(1,
         .1,
         .01,
         .001,
         .0001, 
         .00001,
         .000001)

standards1 <- cbind(standards1, dil)

copy.1 <- c(354287400,   
            35428740, 
            3542874, 
            354287.4, 
            35428.74,
            3542.874,
            354.2874)

# Now we can combine everything and make our standard curve
dat1 <- data.frame('Sample'= standards1$Sample, 'copy' = copy.1, 'log.copy' = log10(copy.1), 'Ct' = standards1$Ct)

fit1 <- lm(Ct ~ log.copy, data = dat1)

line1 <- lm(Ct ~ log.copy, data = dat1)

ab1 <- coef(line1)

ggplot(dat1, aes(x = log.copy, y = Ct, colour= Sample)) +
  theme_classic() +
  xlab(expression('log'['10']*' Copy No. / µL standard')) +
  ylab(expression('C'['T'])) +
  labs(subtitle = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
                        "Intercept =",signif(fit1$coef[[1]],5 ),
                        " Slope =",signif(fit1$coef[[2]], 5),
                        " P =",signif(summary(fit1)$coef[2,4], 5))) +
  geom_point(size = 2) +
  geom_line(aes(group=as.factor(Sample))) +  
  stat_smooth(method = 'lm', formula = y ~ x, level = 0, size = 0.75, col = "black")

convert <- function(y, b, a){
  x <- 10^((y - b) / a)
  x
}

qp1.1$counts <- convert(y = as.numeric(qp1.1$Ct), b = ab1[1], a = ab1[2])


NTC.1 <- qp1.1[grepl('ntc', qp1.1$Sample), ]

NTC.1 <- as.numeric(NTC.1$counts)

NTC.1

qp1.1
