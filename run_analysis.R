## Prior to running this script the data was downloaded from the web using the
## url: 
## "https://d396qusza40orc.cloudfront.net/getdata%2Fprojectfiles%2FUCI%20HAR%20Dataset.zip"
## This downloaded a directory named "UCI HAR Data" onto the desktop which 
## contained informative text files and directories named "train" and "test"
## These directories contained the data that were used by this script to obtain 
## the "tidydata.txt" that was uploaded to the Coursera site
# First download and unzip the dataset
target_url <- "https://d396qusza40orc.cloudfront.net/getdata%2Fprojectfiles%2FUCI%20HAR%20Dataset.zip"
target_localfile = "dataset.zip"
if (!file.exists(target_localfile)) {
        # We wave also metadata, thanks to David Hood for his contribution on this useful piece of code.
        download.file(target_url, method="curl", target_localfile)
        library(tools)       # for md5 checksum
        sink("download_metadata.txt")
        print("Download date:")
        print(Sys.time() )
        print("Download URL:")
        print(target_url)
        print("Downloaded file Information")
        print(file.info(target_localfile))
        print("Downloaded file md5 Checksum")
        print(md5sum(target_localfile))
        sink()
}  ## from student 2 eval

library(data.table)
unzip("dataset.zip")

## reading files needed for test and training data to be merge into one data set
features <- read.table("./UCI HAR Dataset/features.txt")
subject_test <- read.table("./UCI HAR Dataset/test/subject_test.txt")
x_test <-  read.table("./UCI HAR Dataset/test/x_test.txt")
y_test <-  read.table("./UCI HAR Dataset/test/y_test.txt")
subject_train <- read.table("./UCI HAR Dataset/train/subject_train.txt")
x_train <-  read.table("./UCI HAR Dataset/train/x_train.txt")
y_train <-  read.table("./UCI HAR Dataset/train/y_train.txt")

## Merge all components of test and taining data into one data set
## Merge test data into one data set
sub_y_test <- cbind(subject_test, y_test)
colnames(sub_y_test) <- c("subject", "activity") # names first two columns
sub_yx_test <- cbind(sub_y_test, x_test)
## Merge train data into one data set
sub_y_train <- cbind(subject_train, y_train)
colnames(sub_y_train) <- c("subject", "activity") # names first two columns
sub_yx_train <- cbind(sub_y_train, x_train)
## Merge test and training data sets into one large data set
test_train_data <- rbind(sub_yx_train, sub_yx_test)

## Subset the columns containing the mean and standard deviation for each 
## measurement
test_train_subset <- test_train_data[, c(1:8,43:48,83:88,123:128,163:168,
                                        203:204,216:217,229:230,242:243,
                                        255:256,268:273,347:352,426:431,
                                        505:506,518:519,531:532,544:545)]

## Order the data sets by subject and activity
test_train_ordered <- test_train_subset[order(test_train_subset[,1], 
                                        test_train_subset[,2], 
                                        decreasing=FALSE),]

## Replace the activity numbers with descriptive text
test_train_ordered$activity[test_train_ordered$activity==1] <- "walking"
test_train_ordered$activity[test_train_ordered$activity==2] <- "walkingupstairs"
test_train_ordered$activity[test_train_ordered$activity==3] <- "walkingdownstairs"
test_train_ordered$activity[test_train_ordered$activity==4] <- "sitting"
test_train_ordered$activity[test_train_ordered$activity==5] <- "standing"
test_train_ordered$activity[test_train_ordered$activity==6] <- "laying"

## Convert data to a data table and calculate the mean for each measurement for 
## each subject and activity
dt <- data.table(test_train_ordered)
setkey(dt, subject, activity)
dt.meanstdtable <- dt[, list(tBodyAccmeanx=mean(V1),
                                  tBodyAccmeany=mean(V2),
                                  tBodyAccmeanz=mean(V3), 
                                  tBodyAccstdx=mean(V4), 
                                  tBodyAccstdy=mean(V5), 
                                  tBodyAccstdz=mean(V6),
                                  tGravityAccmeanx=mean(V41),
                                  tGravityAccmeany=mean(V42), 
                                  tGravityAccmeanz=mean(V43),
                                  tGravityAccstdx=mean(V44), 
                                  tGravityAccstdy=mean(V45), 
                                  tGravityAccstdz=mean(V46), 
                                  tBodyAccJerkmeanx=mean(V81), 
                                  tBodyAccJerkmeany=mean(V82), 
                                  tBodyAccJerkmeanz=mean(V83), 
                                  tBodyAccJerkstdx=mean(V84), 
                                  tBodyAccJerkstdy=mean(V85), 
                                  tBodyAccJerkstdz=mean(V86), 
                                  tBodyGyromeanx=mean(V121),
                                  tBodyGyromeany=mean(V122),
                                  tBodyGyromeanz=mean(V123), 
                                  tBodyGyrostdx=mean(V124),
                                  tBodyGyrostdy=mean(V125), 
                                  tBodyGyrostdz=mean(V126), 
                                  tBodyGyroJerkmeanx=mean(V161), 
                                  tBodyGyroJerkmeany=mean(V162),
                                  tBodyGyroJerkmeanz=mean(V163), 
                                  tBodyGyroJerkstdx=mean(V164),
                                  tBodyGyroJerkstdy=mean(V165), 
                                  tBodyGyroJerkstdz=mean(V166), 
                                  tBodyAccMagmean=mean(V201), 
                                  tBodyAccMagstd=mean(V202),
                                  tGravityAccMagmean=mean(V214),
                                  tGravityAccMagstd=mean(V215), 
                                  tBodyAccJerkMagmean=mean(V227),
                                  tBodyAccJerkMagstd=mean(V228), 
                                  tBodyGyroMagmean=mean(V240), 
                                  tBodyGyroMagstd=mean(V241), 
                                  tBodyGyroJerkMagmean=mean(V253), 
                                  tBodyGyroJerkMagstd=mean(V254), 
                                  fBodyAccmeanx=mean(V266), 
                                  fBodyAccmeany=mean(V267), 
                                  fBodyAccmeanz=mean(V268), 
                                  fBodyAccstdx=mean(V269), 
                                  fBodyAccstdy=mean(V270),
                                  fBodyAccstdz=mean(V271),
                                  fBodyAccJerkmeanx=mean(V345), 
                                  fBodyAccJerkmeany=mean(V346),
                                  fBodyAccJerkmeanz=mean(V347), 
                                  fBodyAccJerkstdx=mean(V348), 
                                  fBodyAccJerkstdy=mean(V349), 
                                  fBodyAccJerkstdz=mean(V350),
                                  fBodyGyromeanx=mean(V424),
                                  fBodyGyromeany=mean(V425), 
                                  fBodyGyromeanz=mean(V426),
                                  fBodyGyrostdx=mean(V427), 
                                  fBodyGyrostdy=mean(V428), 
                                  fBodyGyrostdz=mean(V429), 
                                  fBodyAccMagmean=mean(V503), 
                                  fBodyAccMagstd=mean(V504), 
                                  fBodyBodyAccJerkMagmean=mean(V516), 
                                  fBodyBodyAccJerkMagstd=mean(V517), 
                                  fBodyBodyGyroMagmean=mean(V529), 
                                  fBodyBodyGyroMagstd=mean(V530), 
                                  fBodyBodyGyroJerkMagmean=mean(V542),
                                  fBodyBodyGyroJerkMagstd=mean(V543)),
                        by=list(subject, activity)]

## Write this data table to a new "tidydata.txt" file 
write.table(dt.meanstdtable, file="tidydata.txt", row.name=FALSE)