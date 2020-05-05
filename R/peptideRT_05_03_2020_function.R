fracpred <- function(){

# This R script allows you to predict the retention times of undetected peptides in a sample based on the retentions times 
# generated from another sample. This method can help select the fractions containing the target peptides during a 2D LC-MS/MS analysis.

library(tidyverse)
cat("Please choose the Excel file that contains your target peptides.","\n")
filecontainingtargets <- file.choose()

pepRT_HLMplasma_1 <- read.csv(filecontainingtargets)

pepRT_HLMplasma_1 <- pepRT_HLMplasma_1 %>% select("Peptide.Retention.Time", "Peptide.Sequence")
pepRT_HLMplasma_1 <- pepRT_HLMplasma_1 %>% rename(RT_HLMplasma_1 = Peptide.Retention.Time, 
                                                  Pepseq = Peptide.Sequence)

cat("Please choose the Excel file that you want to predict the retentiontimes of your target peptides.","\n")
filepredict <- file.choose()
pepRT_plasma_1 <- read.csv(filepredict)

pepRT_plasma_1 <- pepRT_plasma_1 %>% select("Peptide.Retention.Time", "Peptide.Sequence")
pepRT_plasma_1 <- pepRT_plasma_1 %>% rename(RT_plasma_1 = Peptide.Retention.Time, 
                                            Pepseq = Peptide.Sequence)

pepRT_merged <- merge(pepRT_HLMplasma_1, pepRT_plasma_1)
pepRT_merged$RT_HLMplasma_1 <- as.numeric(as.character(pepRT_merged$RT_HLMplasma_1))
pepRT_merged$RT_plasma_1 <- as.numeric(as.character(pepRT_merged$RT_plasma_1))
pepRT_merged_ordered <- pepRT_merged[order(pepRT_merged$RT_HLMplasma_1), ]
pepRT_merged_ordered_rev <- pepRT_merged[order(pepRT_merged$RT_HLMplasma_1), ]

#find the differences between adjacent rows in the plasma column
pepRT_merged_ordered <- pepRT_merged_ordered %>%
  mutate(Diff = RT_plasma_1 - lag(RT_plasma_1))
pepRT_merged_ordered [1, "Diff"] <- 0
count_row_before <- 1
count_row_after <- 2

while (count_row_after != count_row_before){                                   #remove peptides that are out of orders
  count_row_before <- nrow(pepRT_merged_ordered)
  pepRT_merged_ordered <- pepRT_merged_ordered %>% subset(Diff > 0)
  count_row_after <- nrow(pepRT_merged_ordered)
  pepRT_merged_ordered <- pepRT_merged_ordered %>%
    mutate(Diff = RT_plasma_1 - lag(RT_plasma_1))
  pepRT_merged_ordered [1, "Diff"] <- 0.1
} 

#find the differences between adjacent rows in the plasma column in the reversed direction
pepRT_merged_ordered_rev <- pepRT_merged_ordered_rev %>%  mutate(Diff_rev = RT_plasma_1 - lead(RT_plasma_1))
pepRT_merged_ordered_rev [nrow(pepRT_merged_ordered_rev), "Diff_rev"] <- 0
count_row_before_rev <- 1
count_row_after_rev <- 2

while (count_row_after_rev != count_row_before_rev){                                   #remove peptides that are out of orders
  count_row_before_rev <- nrow(pepRT_merged_ordered_rev)
  pepRT_merged_ordered_rev <- pepRT_merged_ordered_rev %>% subset(Diff_rev < 0)
  count_row_after_rev <- nrow(pepRT_merged_ordered_rev)
  pepRT_merged_ordered_rev <- pepRT_merged_ordered_rev %>%
    mutate(Diff_rev = lag(RT_plasma_1) - RT_plasma_1)
  pepRT_merged_ordered_rev [1, "Diff_rev"] <- -0.1
} 

#choose the ordered peptide file with more peptides

if (nrow(pepRT_merged_ordered) > nrow(pepRT_merged_ordered_rev)){
  pepRT_merged_ordered_selected <- pepRT_merged_ordered
} else {pepRT_merged_ordered_selected <- pepRT_merged_ordered_rev}

#keep the peptide sequence and RT columns
pepRT_HLMplasma_2 <- pepRT_HLMplasma_1 %>% select(1, 2)
pepRT_HLMplasma_2 <- pepRT_HLMplasma_2 %>% rename(RT = RT_HLMplasma_1)
pepRT_plasma_2 <- pepRT_plasma_1 %>% select(1, 2)
pepRT_plasma_2 <- pepRT_plasma_2 %>% rename(RT = RT_plasma_1)

#remove duplicate rows
pepDup <- duplicated(pepRT_HLMplasma_2$Pepseq)
pepRT_HLMplasma_2 <- cbind(pepRT_HLMplasma_2, pepDup)
pepRT_HLMplasma_2 <- pepRT_HLMplasma_2  %>% subset(pepDup == FALSE)

#input the seqeunces of the target peptide, test CES1 peptide: AISESGVALTSVLVK
candidatepep <- readline(prompt="Please input the seqeunces of your target peptide: ")

if (mean(grepl(candidatepep, pepRT_HLMplasma_1$Pepseq)) == 0){
  cat("The peptide was not found in your file!", "\n")
  predicted_RT_results <- data.frame("Peptide" = candidatepep, "RT" = 0)
} else {cat("The peptide was found in your file!", "\n")
  
#extract the target peptide from the mix sample and made the dataframe consistent with the dataframe to be compared

targetPep <- pepRT_HLMplasma_2 %>% filter (Pepseq == candidatepep) %>% select(1,2) %>% rename(RT_HLMplasma_1 = RT)
targetPep$RT_plasma_1 <- 0

pepRT_merged_ordered_selected <- pepRT_merged_ordered_selected %>% select(Pepseq, RT_HLMplasma_1, RT_plasma_1)
targetPep <- rbind(pepRT_merged_ordered_selected, targetPep)

targetPep$RT_HLMplasma_1 <- as.numeric(as.character(targetPep$RT_HLMplasma_1))
targetPep$RT_plasma_1 <- as.numeric(as.character(targetPep$RT_plasma_1))
targetPep <- targetPep[order(targetPep$RT_HLMplasma_1),]

#find the row # n from the merged dataframe
n <- which(grepl(candidatepep, targetPep$Pepseq))

#calculate the predicted RT of the target peptide
if (n == 1) {predicted_RT <- (targetPep[n+1, "RT_plasma_1"]-
                               ((targetPep[n+2, "RT_plasma_1"]-targetPep[n+1, "RT_plasma_1"])*(targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n, "RT_HLMplasma_1"])/
                                  (targetPep[n+2, "RT_HLMplasma_1"]-targetPep[n+1, "RT_HLMplasma_1"])))}
else if (n == nrow(targetPep)) {predicted_RT <- (targetPep[n-1, "RT_plasma_1"]+
                                                 ((targetPep[n-1, "RT_plasma_1"]-targetPep[n-2, "RT_plasma_1"])*
                                                    (targetPep[n, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])/
                                                    (targetPep[n-1, "RT_HLMplasma_1"]-targetPep[n-2, "RT_HLMplasma_1"])))}
else {predicted_RT <- ((targetPep[n, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])*(targetPep[n+1, "RT_plasma_1"]-targetPep[n-1, "RT_plasma_1"])+
    (targetPep[n-1, "RT_plasma_1"])*(targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"]))/
  (targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])}

cat("The predicted retention time of the peptide", candidatepep, "is", predicted_RT, "minutes.", "\n")


#create a dataframe of the results
predicted_RT_results <- data.frame("Peptide" = candidatepep, "RT" = predicted_RT)
}

#input the seqeunces of MORE target peptide, test CES1 peptide: AISESGVALTSVLVK
candidatepep_more <- readline(prompt="Do you want to predict the retention times of more peptides? How many?")

if (candidatepep_more == 0){
  print("Thank you very much!")
} else {
  for (candidatepep_more in 1:candidatepep_more){
    
    #input the seqeunces of the target peptide, test CES1 peptides: NLSVEDAAR, RSTVAQLVK, AISESGVALTSVLVK, 
    # FTPPQPAEPWSFVK, ESQPLLGTVIDGMLLLK, NFLAFIQHLR
    candidatepep <- readline(prompt="Please input the seqeunces of your target peptide: ")
  
    if (mean(grepl(candidatepep, pepRT_HLMplasma_1$Pepseq)) == 0){
      cat("The peptide was not found in your file!", "\n")
      predicted_RT_results_more2 <- data.frame("Peptide" = candidatepep, "RT" = 0)
      predicted_RT_results <- rbind(predicted_RT_results, predicted_RT_results_more2)
    } else {cat("The peptide was found in your file!", "\n")
      
      #extract the target peptide from the mix sample and made the dataframe concsistent the dataframe to be compared
      
      targetPep <- pepRT_HLMplasma_2 %>% filter (Pepseq == candidatepep) %>% select(1,2) %>% rename(RT_HLMplasma_1 = RT)
      targetPep$RT_plasma_1 <- 0
      
      pepRT_merged_ordered_selected <- pepRT_merged_ordered_selected %>% select(Pepseq, RT_HLMplasma_1, RT_plasma_1)
      targetPep <- rbind(pepRT_merged_ordered_selected, targetPep)
      
      targetPep$RT_HLMplasma_1 <- as.numeric(as.character(targetPep$RT_HLMplasma_1))
      targetPep$RT_plasma_1 <- as.numeric(as.character(targetPep$RT_plasma_1))
      targetPep <- targetPep[order(targetPep$RT_HLMplasma_1),]
      
      
      #find the row # n from the merged dataframe
      n <- which(grepl(candidatepep, targetPep$Pepseq))
      
      #calculate the predicted RT of the target peptide
      if (n == 1) {predicted_RT <- (targetPep[n+1, "RT_plasma_1"]-
                                      ((targetPep[n+2, "RT_plasma_1"]-targetPep[n+1, "RT_plasma_1"])*(targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n, "RT_HLMplasma_1"])/
                                         (targetPep[n+2, "RT_HLMplasma_1"]-targetPep[n+1, "RT_HLMplasma_1"])))}
      else if (n == nrow(targetPep)) {predicted_RT <- (targetPep[n-1, "RT_plasma_1"]+
                                                         ((targetPep[n-1, "RT_plasma_1"]-targetPep[n-2, "RT_plasma_1"])*
                                                            (targetPep[n, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])/
                                                            (targetPep[n-1, "RT_HLMplasma_1"]-targetPep[n-2, "RT_HLMplasma_1"])))}
      else {predicted_RT <- ((targetPep[n, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])*(targetPep[n+1, "RT_plasma_1"]-targetPep[n-1, "RT_plasma_1"])+
                               (targetPep[n-1, "RT_plasma_1"])*(targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"]))/
        (targetPep[n+1, "RT_HLMplasma_1"]-targetPep[n-1, "RT_HLMplasma_1"])}
      
      cat("The predicted retention time of the peptide", candidatepep, "is", predicted_RT, "minutes.", "\n")
     
    
      predicted_RT_results_more <- data.frame("Peptide" = candidatepep, "RT" = predicted_RT)
      predicted_RT_results <- rbind(predicted_RT_results, predicted_RT_results_more)
    }
    
      }
  
}

print(predicted_RT_results)
cat("The results have been saved to predicted_RT.csv.", "\n")
write.csv(predicted_RT_results, "predicted_RT.csv")
}

