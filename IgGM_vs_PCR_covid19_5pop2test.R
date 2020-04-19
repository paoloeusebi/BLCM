#install.packages("runjags", "rjags")

library(runjags)
testjags()

covidIgG_M <- "
model
{
  #model specification
  
  y1[1:4] ~ dmulti(p1[1:4], n1)
  
  p1[1] <- pi1*Se[1]*Se[2] + (1-pi1)*(1-Sp[1])*(1-Sp[2])
  p1[2] <- pi1*Se[1]*(1-Se[2]) + (1-pi1)*(1-Sp[1])*Sp[2]
  p1[3] <- pi1*(1-Se[1])*Se[2] + (1-pi1)*Sp[1] *(1-Sp[2])
  p1[4] <- pi1*(1-Se[1])*(1-Se[2]) + (1-pi1)*Sp[1]*Sp[2]
  
  pi1 ~dbeta(1, 1)
  
  y2[1:4] ~ dmulti(p2[1:4], n2)
  
  p2[1] <- pi2*Se[1]*Se[2] + (1-pi2)*(1-Sp[1])*(1-Sp[2])
  p2[2] <- pi2*Se[1]*(1-Se[2]) + (1-pi2)*(1-Sp[1])*Sp[2]
  p2[3] <- pi2*(1-Se[1])*Se[2] + (1-pi2)*Sp[1] *(1-Sp[2])
  p2[4] <- pi2*(1-Se[1])*(1-Se[2]) + (1-pi2)*Sp[1]*Sp[2]
  
  pi2 ~dbeta(1, 1)
  
  y3[1:4] ~ dmulti(p3[1:4], n3)
  
  p3[1] <- pi3*Se[1]*Se[2] + (1-pi3)*(1-Sp[1])*(1-Sp[2])
  p3[2] <- pi3*Se[1]*(1-Se[2]) + (1-pi3)*(1-Sp[1])*Sp[2]
  p3[3] <- pi3*(1-Se[1])*Se[2] + (1-pi3)*Sp[1] *(1-Sp[2])
  p3[4] <- pi3*(1-Se[1])*(1-Se[2]) + (1-pi3)*Sp[1]*Sp[2]
  
  pi3 ~dbeta(1, 1)
  
  y4[1:4] ~ dmulti(p4[1:4], n4)
  
  p4[1] <- pi4*Se[1]*Se[2] + (1-pi4)*(1-Sp[1])*(1-Sp[2])
  p4[2] <- pi4*Se[1]*(1-Se[2]) + (1-pi4)*(1-Sp[1])*Sp[2]
  p4[3] <- pi4*(1-Se[1])*Se[2] + (1-pi4)*Sp[1] *(1-Sp[2])
  p4[4] <- pi4*(1-Se[1])*(1-Se[2]) + (1-pi4)*Sp[1]*Sp[2]
  
  pi4 ~dbeta(1, 1)
  
  y5[1:4] ~ dmulti(p5[1:4], n5)
  
  p5[1] <- pi5*Se[1]*Se[2] + (1-pi5)*(1-Sp[1])*(1-Sp[2])
  p5[2] <- pi5*Se[1]*(1-Se[2]) + (1-pi5)*(1-Sp[1])*Sp[2]
  p5[3] <- pi5*(1-Se[1])*Se[2] + (1-pi5)*Sp[1] *(1-Sp[2])
  p5[4] <- pi5*(1-Se[1])*(1-Se[2]) + (1-pi5)*Sp[1]*Sp[2]
  
  pi5 ~dbeta(1, 1)
  
  Se[1] ~ dbeta(1,1) I(0.1,)
  Sp[1] ~ dbeta(1,1)  I(0.1,)
  Se[2] ~ dbeta(1,1) I(0.1,)
  Sp[2] ~ dbeta(1,1) I(0.1,)

PPVp05[1] <- Se[1]*0.05/(Se[1]*0.05+(1-Sp[1])*(1-0.05))
NPVp05[1] <- Sp[1]*(1-0.05)/((1-Se[1])*0.05+Sp[1]*(1-0.05))

PPVp10[1] <- Se[1]*0.10/(Se[1]*0.10+(1-Sp[1])*(1-0.10))
NPVp10[1] <- Sp[1]*(1-0.10)/((1-Se[1])*0.10+Sp[1]*(1-0.10))

PPVp15[1] <- Se[1]*0.15/(Se[1]*0.15+(1-Sp[1])*(1-0.15))
NPVp15[1] <- Sp[1]*(1-0.15)/((1-Se[1])*0.15+Sp[1]*(1-0.15))
  
}
"

DataIgG_M=list(n1=525, n2=57, n3=190, n4=278, n5=126,
               y1=c(352,	45,	12,	116), # Li 
               y2=c(21, 3, 24, 9), # Jia
               y3=c(21,49,13,107), # Paradiso
               y4=c(91, 7, 6, 174), # Cellex
               y5=c(103,	9,	0,	14)) #Aytu

resultsIgG_M <- run.jags(covidIgG_M, 
                         data=DataIgG_M,
                         monitor = c("Se","Sp",
                                     "PPVp05[1]", "NPVp05[1]",
                                     "PPVp10[1]", "NPVp10[1]",
                                     "PPVp15[1]", "NPVp15[1]"),
                         n.chains = 3)

print(resultsIgG_M)
plot(resultsIgG_M, vars = c("Se", "Sp"),
     plot.type = c("density", "trace"))


# Gold standard

covidIgG_M_gold <- "
model
{
  #model specification
  
  y1[1:4] ~ dmulti(p1[1:4], n1)
  
  p1[1] <- pi1*Se[1] + (1-pi1)*(1-Sp[1])
  p1[2] <- (1-pi1)*(1-Sp[1])
  p1[3] <- pi1*(1-Se[1]) + (1-pi1)*Sp[1] 
  p1[4] <- (1-pi1)*Sp[1]
  
  pi1 ~dbeta(1, 1)
  
  y2[1:4] ~ dmulti(p2[1:4], n2)
  
  p2[1] <- pi2*Se[1] + (1-pi2)*(1-Sp[1])
  p2[2] <- (1-pi2)*(1-Sp[1])
  p2[3] <- pi2*(1-Se[1]) + (1-pi2)*Sp[1] 
  p2[4] <- (1-pi2)*Sp[1]
  
  pi2 ~dbeta(1, 1)
  
  y3[1:4] ~ dmulti(p3[1:4], n3)
  
  p3[1] <- pi3*Se[1] + (1-pi3)*(1-Sp[1])
  p3[2] <- (1-pi3)*(1-Sp[1])
  p3[3] <- pi3*(1-Se[1]) + (1-pi3)*Sp[1] 
  p3[4] <- (1-pi3)*Sp[1]
  
  pi3 ~dbeta(1, 1)
  
  y4[1:4] ~ dmulti(p4[1:4], n4)
  
  p4[1] <- pi4*Se[1] + (1-pi4)*(1-Sp[1])
  p4[2] <- (1-pi4)*(1-Sp[1])
  p4[3] <- pi4*(1-Se[1]) + (1-pi4)*Sp[1] 
  p4[4] <- (1-pi4)*Sp[1]
  
  pi4 ~dbeta(1, 1)
  
  y5[1:4] ~ dmulti(p5[1:4], n5)
  
  p5[1] <- pi5*Se[1] + (1-pi5)*(1-Sp[1])
  p5[2] <- (1-pi5)*(1-Sp[1])
  p5[3] <- pi5*(1-Se[1]) + (1-pi5)*Sp[1] 
  p5[4] <- (1-pi5)*Sp[1]
  
  pi5 ~dbeta(1, 1)
  
  Se[1] ~ dbeta(1,1) I(0.1,)
  Sp[1] ~ dbeta(1,1)  I(0.1,)
  
  PPVp05[1] <- Se[1]*0.05/(Se[1]*0.05+(1-Sp[1])*(1-0.05))
  NPVp05[1] <- Sp[1]*(1-0.05)/((1-Se[1])*0.05+Sp[1]*(1-0.05))
  
  PPVp10[1] <- Se[1]*0.10/(Se[1]*0.10+(1-Sp[1])*(1-0.10))
  NPVp10[1] <- Sp[1]*(1-0.10)/((1-Se[1])*0.10+Sp[1]*(1-0.10))
  
}
"

DataIgG_M=list(n1=525, n2=57, n3=190, n4=278, n5=126,
               y1=c(352,	45,	12,	116), # Li 
               y2=c(21, 3, 24, 9), # Jia
               y3=c(21,49,13,107), # Paradiso
               y4=c(91, 7, 6, 174), # Cellex
               y5=c(103,	9,	0,	14)) #Aytu

resultsIgG_M_gold <- run.jags(covidIgG_M_gold, 
                         data=DataIgG_M,
                         monitor = c("Se","Sp",
                                     "PPVp05[1]", "NPVp05[1]",
                                     "PPVp10[1]", "NPVp10[1]"),
                         n.chains = 3)

print(resultsIgG_M_gold)
plot(resultsIgG_M_gold, vars = c("Se", "Sp"),
     plot.type = c("density", "trace"))