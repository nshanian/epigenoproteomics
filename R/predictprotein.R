#' Predict protein abundance
#'
#' Calculate protein predictive significance for all genes
#'
#' @param Peak_annotation A tibble data frame that lists assigned genes to peaks and their status as differential features
#' @param PeakCount A tibble data frame that lists assigned genes to peaks and their accessibility counts for all samples. Sample columns have the same order as TranscriptCount and ProteinCount. Make sure to use same columns for annotations.
#' @param TranscriptCount A tibble data frame that lists genes and their RNA expression in log transformation of transcript per million units for all samples. Sample columns have the same order as PeakCount and ProteinCount
#' @param ProteinCount A tibble data frame that lists genes and their protein expression in log transformation of normalized PSM counts for all samples. Sample columns have the same order as PeakCount and TranscriptCount
#' @param orderSamples A vector of sample IDs in the order of the sample columns in Count data frames (i.e. PeakCount, TranscriptCount, and ProteinCount)
#' @param graphSamples A vector of sample IDs randomized to test group, which was output from the epigenoproteomics::randomize(orderSamples)
#' @return A data frame with Peak annotations and p-value of protein predicitive significance
#' @examples
#' predictprotein(data/Peak_annotation.rda,data/PeakCount_NonPromoter.rda,data/TranscriptCount.rda,data/ProteinCount.rda,data/orderSamples.rda,c("1","2"))
predictprotein <- function(Peak_annotation,PeakCount,TranscriptCount,ProteinCount,orderSamples,testSamples){
  PeakGene_lmModel_predict <- c()
  for(p in 1:nrow(Peak_annotation)){
    gene <- unlist(unname(Peak_annotation[p,1]))
    rowName <- unlist(unname(Peak_annotation[p,2]))
    status <- unlist(unname(Peak_annotation[p,3]))
    transcript_count <- c()
    protein_count <- c()
    peak_count <- c()
    transcript_count <- unlist(unname(dplyr::filter(TranscriptCount,Gene==gene)[1,2:ncol(TranscriptCount)]))
    protein_count <- unlist(unname(dplyr::filter(ProteinCount,Gene==gene)[1,2:ncol(ProteinCount)]))
    peak_count <- unlist(unname(dplyr::filter(PeakCount,Gene==gene,peakName==rowName)[1,3:ncol(PeakCount)]))
    cor_transProt_Peak <- dplyr::inner_join(dplyr::bind_cols(dplyr::tibble(transcript=transcript_count),dplyr::tibble(peak=peak_count),dplyr::tibble(protein=protein_count,samples=orderSamples)),dplyr::tibble(samples=testSamples))
    PeakGene_lmModel_predict[p] <- unlist(unname(stats::pf(summary(stats::lm(protein ~ transcript+transcript*peak,data=cor_transProt_Peak))$fstatistic[1],summary(stats::lm(protein ~ transcript+transcript*peak,data=cor_transProt_Peak))$fstatistic[2],summary(stats::lm(protein ~ transcript+transcript*peak,data=cor_transProt_Peak))$fstatistic[3],lower.tail=F)))
  }
  PeakGene_lmModel_predict_genes <- dplyr::mutate(Peak_annotation,lm_Ftest=PeakGene_lmModel_predict)
  PeakGene_lmModel_predict_genes
}




