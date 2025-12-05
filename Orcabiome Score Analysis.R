### Orcabiome Score Analysis ###
#Created by: Brittany Visona-Kelly
#Created: June 3, 2025

###########################
## Image Score Analysis
###########################

score <- read.csv(("Orcabiome score analysis.csv"), stringsAsFactors =TRUE, na.strings = c(""))
str(score)

score$S1.percentlesions.avg <- (as.numeric(score$S1.lesions.pixels.avg) / as.numeric(score$S1.allpixels.avg))*100
score$S2.percentlesions.avg <- (as.numeric(score$S2.lesions.pixels.avg) / as.numeric(score$S2.allpixels.avg))*100
score$S3.percentlesions.avg <- (as.numeric(score$S3.lesions.pixels.avg) / as.numeric(score$S3.allpixels.avg))*100
score$S4.percentlesions.avg <- (as.numeric(score$S4.lesions.pixels.avg) / as.numeric(score$S4.allpixels.avg))*100

### Lesion score test ###
model1 <- aov(S1.percentlesions.avg ~ as.factor(S1.score), data=score)
summary(model1)
#                    Df    Sum Sq   Mean Sq F value Pr(>F)  
#as.factor(S1.score)  2 0.0005361 2.680e-04   4.616 0.0189 *
#  Residuals           27 0.0015679 5.807e-05                           

ls1<-TukeyHSD(model1, conf.level=.95)
ls1
#             diff          lwr          upr     p adj
#2-1 -0.0122582492 -0.023061375 -0.001455124 0.0237507**
#3-1 -0.0125951395 -0.023398265 -0.001792014 0.0198541**
#3-2 -0.0003368903 -0.007747773  0.007073993 0.9930225

model2<-aov((S2.percentlesions.avg~as.factor(S2.score)),data=score)
summary(model2)
#                    Df   Sum Sq   Mean Sq F value Pr(>F)
#as.factor(S2.score)  2 0.000932 0.0004660   2.398   0.11
#Residuals           27 0.005246 0.0001943                     

ls2<-TukeyHSD(model2, conf.level=.95)
ls2
#            diff         lwr         upr     p adj
#2-1 -0.013847867 -0.03248928 0.004793548 0.1753499
#3-1 -0.015565699 -0.03357215 0.002440747 0.0999057
#3-2 -0.001717833 -0.01564330 0.012207638 0.9498346

model3<-aov((S3.percentlesions.avg~as.factor(S3.score)),data=score)
summary(model3)
#                   Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S3.score)  2   27.1   13.57   0.822   0.45
#Residuals           27  445.7   16.51  

ls3<-TukeyHSD(model3, conf.level=.95)
ls3
#        diff       lwr       upr     p adj
#2-1 1.132387 -6.482412  8.747187 0.9280016
#3-1 2.816779 -4.798021 10.431578 0.6343467
#3-2 1.684391 -2.123008  5.491791 0.5241680

model4<-aov((S4.percentlesions.avg~as.factor(S4.score)),data=score)
summary(model4)
#                   Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S4.score)  2    254   127.2   1.028  0.371
#Residuals           27   3342   123.8    

ls4<-TukeyHSD(model4, conf.level=.95)
ls4
#        diff        lwr      upr     p adj
#2-1 1.489643 -12.755091 15.73438 0.9636677
#3-1 7.192657  -7.916165 22.30148 0.4747692
#3-2 5.703014  -5.558437 16.96446 0.4318418

### Scar score test ###
model5 <- aov(as.integer(S1.scar.score.avg)  ~ as.factor(S1.score), data=score)
summary(model5)
#                    Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S1.score)  2  1.828  0.9141   1.407  0.262
#Residuals           27 17.538  0.6496                             

ls5<-TukeyHSD(model5, conf.level=.95)
ls5
#diff        lwr       upr     p adj
#2-1 -0.7692308 -1.9118101 0.3733485 0.2351900
#3-1 -0.5384615 -1.6810408 0.6041178 0.4816915
#3-2  0.2307692 -0.5530337 1.0145722 0.7480309

model6 <- aov(as.integer(S2.scar.score.avg)  ~ as.factor(S2.score), data=score)
summary(model6)
#                    Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S2.score)  2  0.306  0.1532   0.921   0.41
#Residuals           27  4.494  0.1664                           

ls6<-TukeyHSD(model6, conf.level=.95)
ls6
#        diff        lwr       upr     p adj
#2-1 0.1818182 -0.3637376 0.7273739 0.6902170
#3-1 0.2857143 -0.2412586 0.8126872 0.3837337
#3-2 0.1038961 -0.3036438 0.5114360 0.8038790

model7 <- aov(as.integer(S3.scar.score.avg)  ~ as.factor(S3.score), data=score)
summary(model7)
#                    Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S3.score)  2  0.295  0.1476   0.264   0.77
#Residuals           27 15.071  0.5582                        

ls7<-TukeyHSD(model7, conf.level=.95)
ls7
#         diff        lwr       upr     p adj
#2-1 0.2142857 -1.1860301 1.6146015 0.9239384
#3-1 0.3571429 -1.0431730 1.7574587 0.8037293
#3-2 0.1428571 -0.5573008 0.8430151 0.8691095

model8 <- aov(as.integer(S4.scar.score.avg)  ~ as.factor(S4.score), data=score)
summary(model8)
#                     Df Sum Sq Mean Sq F value Pr(>F)
#as.factor(S4.score)  2  0.033  0.0167   0.019  0.981
#Residuals           27 23.433  0.8679                               

ls8<-TukeyHSD(model8, conf.level=.95)
ls8
#          diff        lwr      upr     p adj
#2-1 0.06666667 -1.1261377 1.259471 0.9894730
#3-1 0.10000000 -1.1651601 1.365160 0.9790659
#3-2 0.03333333 -0.9096613 0.976328 0.9957747


### END ###