
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
C8 <- readRDS(args[14])
C9 <- readRDS(args[15])
C10 <- readRDS(args[16])
C11 <- readRDS(args[17])
C12 <- readRDS(args[18])
C13 <- readRDS(args[19])
C14 <- readRDS(args[20])
C15 <- readRDS(args[21])




# Convert to BacArena models

# Add the organisms
bac1 <- BacArena::Bac(C1,type=args[22])
bac2 <- BacArena::Bac(C2,type=args[23])
bac3 <- BacArena::Bac(C3,type=args[24])
bac4 <- BacArena::Bac(C4,type=args[25])
bac5 <- BacArena::Bac(C5,type=args[26])
bac6 <- BacArena::Bac(C6,type=args[27])
bac7 <- BacArena::Bac(C7,type=args[28])
bac8 <- BacArena::Bac(C8,type=args[29])
bac9 <- BacArena::Bac(C9,type=args[30])
bac10 <- BacArena::Bac(C10,type=args[31])
bac11 <- BacArena::Bac(C11,type=args[32])
bac12 <- BacArena::Bac(C12,type=args[33])
bac13 <- BacArena::Bac(C13,type=args[34])
bac14 <- BacArena::Bac(C14,type=args[35])
bac15 <- BacArena::Bac(C15,type=args[36])

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
arena <- BacArena::addDefaultMed(arena, bac8) # add default media
arena <- BacArena::addDefaultMed(arena, bac9) # add default media
arena <- BacArena::addDefaultMed(arena, bac10) # add default media
arena <- BacArena::addDefaultMed(arena, bac11) # add default media
arena <- BacArena::addDefaultMed(arena, bac12) # add default media
arena <- BacArena::addDefaultMed(arena, bac13) # add default media
arena <- BacArena::addDefaultMed(arena, bac14) # add default media
arena <- BacArena::addDefaultMed(arena, bac15) # add default media

# Simulate
sim <- BacArena::simEnv(arena, time=7)

p<-plotGrowthCurve(sim)

data = p[[1]]["data"]

write.csv(data, file = sprintf("statistics.csv"))

