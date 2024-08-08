
# Load modules
library("BacArena")

args = commandArgs()

# Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

setwd(args[6])  #<--CHANGE ACCORDINGLY !!!


C1 <- readRDS(args[7])
C2 <- readRDS(args[8])
C3 <- readRDS(args[9])
C4 <- readRDS(args[10])
C5 <- readRDS(args[11])
C6 <- readRDS(args[12])
C7 <- readRDS(args[13])




# Convert to BacArena models

# Add the organisms
bac1 <- BacArena::Bac(C1,type=args[14])
bac2 <- BacArena::Bac(C2,type=args[15])
bac3 <- BacArena::Bac(C3,type=args[16])
bac4 <- BacArena::Bac(C4,type=args[17])
bac5 <- BacArena::Bac(C5,type=args[18])
bac6 <- BacArena::Bac(C6,type=args[19])
bac7 <- BacArena::Bac(C7,type=args[20])

# Set size
arena <- BacArena::Arena(n=100, m=100)

#Add organisms
arena <- BacArena::addOrg(arena, bac1, amount=10)
arena <- BacArena::addOrg(arena, bac2, amount=10)

# Add media for all 
arena <- BacArena::addDefaultMed(arena, bac1) # add default media
arena <- BacArena::addDefaultMed(arena, bac2) # add default media
arena <- BacArena::addDefaultMed(arena, bac3) # add default media
arena <- BacArena::addDefaultMed(arena, bac4) # add default media
arena <- BacArena::addDefaultMed(arena, bac5) # add default media
arena <- BacArena::addDefaultMed(arena, bac6) # add default media
arena <- BacArena::addDefaultMed(arena, bac7) # add default media

# Simulate
sim <- BacArena::simEnv(arena, time=7)

p<-plotGrowthCurve(sim)

data = p[[1]]["data"]

write.csv(data, file = sprintf("statistics.csv"))

