# demog
demog <- read.csv("/mnt/nfs/class/exampledata/fmri/adhd200/command/demog.txt", sep="\t")
## peking_1, peking_2, peking_3
peking1 <- list.files("/mnt/nfs/class/exampledata/fmri/adhd200/data/athena/training/Peking_1", pattern="^[0-9]+")
peking2 <- list.files("/mnt/nfs/class/exampledata/fmri/adhd200/data/athena/training/Peking_2", pattern="^[0-9]+")
peking3 <- list.files("/mnt/nfs/class/exampledata/fmri/adhd200/data/athena/training/Peking_3", pattern="^[0-9]+")
## test release 
csvfiles <- Sys.glob("/mnt/nfs/class/exampledata/fmri/adhd200/data/athena/test/*/*.csv")
site_choices <- c("Peking_1", "Brown", "KKI", "NI", "NYU", "OHSU", "Pitt", "WashU")
gender_choices <- c("Female", "Male")
handedness_choices <- c("Left", "Right", "Ambidextrous")
testdfs <- ldply(csvfiles, read.csv)
testdfs[c(71,87),]$Handedness <- 0    # ghetto fix
testdfs$Site <- sapply(testdfs$Site, function(i) site_choices[i])
testdfs$Gender <- sapply(testdfs$Gender, function(i) gender_choices[i+1])
testdfs$Handedness <- sapply(testdfs$Handedness, function(i) handedness_choices[i+1])
## combine
tmp1 <- cbind(demog, type=rep("training", nrow(demog)))
sites <- as.character(tmp1$Site)
sites[sites=="Peking"] <- sapply(tmp1$ScanDir.ID[tmp1$Site=="Peking"], function(id) {
    if (any(grepl(id, peking1)))
        return("Peking_1")
    else if ((any(grepl(id, peking2))))
        return("Peking_2")
    else if ((any(grepl(id, peking3))))
        return("Peking_3")
    else
        stop("Unrecognized id ", id)
})
tmp1$Site <- factor(sites)
tmp2 <- cbind(testdfs[,1:17], DX2=testdfs$DX, type=rep("test", nrow(testdfs)))
colnames(tmp2)[7] <- "Secondary.Dx"
d <- rbind(tmp1, tmp2)
## fixes
colnames(d)[1] <- "ID"
d$ID <- sapply(d$ID, function(x) sprintf("%07i", as.integer(x)))
d$Full4.IQ[is.na(d$Full4.IQ)] <- 1
d$Gender[is.na(d$Gender)] <- sample(levels(d$Gender), sum(is.na(d$Gender)))
## remove WashU
d <- d[d$Site!="WashU",]
d$Site <- factor(d$Site)
## save
write.table(d, sep="\t", file="demog_run1_all.txt", row.names=F, quote=T)

# setup subjects
fnames <- Sys.glob("/mnt/nfs/class/exampledata/fmri/adhd200/data/niak/rois_3000/*/*/tseries*mat")
exclude <- c("0010003", "0010015", "3662296", "6568351", "0016017", "0026041", "0021017", "0021029", "0021036")
## remove bad subs
w <- unlist(lapply(exclude, function(x) grep(x, fnames)))
fnames <- fnames[-w]
## restrict to 1st run
w <- grep("run1", fnames)
fnames <- fnames[w]
## save
write.table(fnames, file="tseries_fnames.txt", row.names=F, col.names=F, quote=F)

# mismatch?
## missing
tmp <- sapply(d$ID, function(id) any(grepl(id, fnames)))
w1 <- which(!tmp)
## exclude
w2 <- sapply(exclude, function(e) which(d$ID==e))
## save
w <- c(w1, w2)
badsubs <- d$ID[w]
write.table(badsubs, file="excluded_subs.txt", row.names=F, col.names=F, quote=F)
write.table(d[w,], file="demog_run1_all_excluded_subs.txt", sep="\t", row.names=F)
## filter
dfilt <- d[-w,]
write.table(dfilt, sep="\t", file="demog_run1_all_included.txt", row.names=F, quote=T)

# design matrix
mf <- model.frame(ID ~ Site + Gender + Age + Full4.IQ + type, dfilt)
X <- model.matrix(attr(mf, "terms"), data=mf)
if (nrow(X) != nrow(dfilt))
    warning("Size mismatch in create of design matrix", immediate. = TRUE)
write.table(as.matrix(X), file="design.txt", row.names=F, col.names=T, quote=F)
