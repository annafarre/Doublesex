# install.packages("UpSetR")
library(UpSetR)

Ame_diffMatches.caste <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_caste.txt", header = T)
Ame_diffMatches_2_4days.caste <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_2_4_days_caste.txt", header = T)
Ame_diffMatches_2days.caste <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_2days_caste.txt", header = T)
Ame_diffMatches_4days.caste <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_4days_caste.txt", header = T)
Wau_diffMatches.caste <-  read.table("~/Unimelb/Doublesex/R/Wau_diffMatches_caste.txt", header = T)
Mph_diffMatches.caste <-  read.table("~/Unimelb/Doublesex/R/Mph_diffMatches_DRP002877_caste.txt", header = T)
Aec_diffMatches.caste <-  read.table("~/Unimelb/Doublesex/R/Aec_diffMatches_casteGyneWorker.txt", header = T)


listInput <- list(Ame.prepupa=Ame_diffMatches.caste$unlist.Ame_diffMatches.caste., 
                  Ame.2day=Ame_diffMatches_2days.caste$unlist.Ame_diffMatches.caste.,
                  Ame.4day=Ame_diffMatches_4days.caste$unlist.Ame_diffMatches.caste.,
                  Wau=Wau_diffMatches.caste$unlist.Wau_diffMatches.caste., 
                  Aec=Aec_diffMatches.caste$unlist.Aec_diffMatches.sex., 
                  Mph=Mph_diffMatches.caste$unlist.Mph_diffMatches.DRP002877.caste.)

upset(fromList(listInput), order.by = "freq", empty.intersections = T)

Ame_diffMatches.sex <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_sex.txt", header = T)
Ame_diffMatches_2_4days.sex <- read.table("~/Unimelb/Doublesex/R/Ame_diffMatches_2_4_days_sex.txt", header = T)
Wau_diffMatches.sex <-  read.table("~/Unimelb/Doublesex/R/Wau_diffMatches_sex.txt", header = T)

listInput.sex <- list(Ame.prepupa=Ame_diffMatches.sex$unlist.Ame_diffMatches.sex., 
                  Ame.larva=Ame_diffMatches_2_4days.sex$unlist.Ame_diffMatches.sex.,
                  Wau=Wau_diffMatches.sex$unlist.Wau_diffMatches.sex.)

upset(fromList(listInput.sex), order.by = "freq")

length(intersect.Vector(Ame_diffMatches.sex$unlist.Ame_diffMatches.sex., Ame_diffMatches_2_4days.sex$unlist.Ame_diffMatches.sex.))
