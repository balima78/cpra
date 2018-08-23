#' cPRA value for a single candidate
#'
#' This function computes a cPRA value for one candidate with their HLA sensitization
#' and given HLA allelic and haplotypic frequencies
#'
#' @param id a numeric valic with patients identification
#' @param sensib a 2 columns data frame with patients identification 'ID'
#' and HLA antibodies 'acs'
#' @param A_freq a 2 columns dataframe with HLA-A alelic frequencies.
#' A column (named 'AA') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param B_freq a 2 columns dataframe with HLA-B alelic frequencies.
#' A column (named 'BB') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param DR_freq a 2 columns dataframe with HLA-DR alelic frequencies.
#' A column (named 'DDR') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param AB_freq a 3 columns dataframe with HLA-A-B haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles, another (named 'BB') with HLA-B alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param ADR_freq a 3 columns dataframe with HLA-A-DR haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles, another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param BDR_freq a 3 columns dataframe with HLA-B-DR haplotypic frequencies.
#' A column (named 'BB') with HLA-B alleles, another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param ABDR_freq a 4 columns dataframe with HLA-A-B-DR haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles,
#' another column (named 'BB') with HLA-B alleles,
#' another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @return a percentual value corresponding to cPRA value
#' @author Bruno A Lima
#' @details This function returns a single percentage for cPRA of a single candidate
#' considering their HLA immunization and HLA allelic and haplotypic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors.
#' @export
#' @import tidyverse

cpra1<-function(id = 1, sensib = acHLA,
                A_freq = HLA_A,
                B_freq = HLA_B,
                DR_freq = HLA_DR,
                AB_freq = HLA_AB,
                ADR_freq = HLA_ADR,
                BDR_freq = HLA_BDR,
                ABDR_freq = HLA_ABDR){


  ac<-sensib %>% filter(ID == id) %>% select(acs)

  Sa<-if(sum(ac$acs %in% A_freq$AA)>=1)
  {sum(A_freq %>% filter(AA %in% t(ac)) %>% select(freq))}
  else Sa<-0

  Sb<-if(sum(ac$acs %in% B_freq$BB)>=1)
  {sum(B_freq %>% filter(BB %in% t(ac)) %>% select(freq))}
  else Sb<-0

  Sdr<-if(sum(ac$acs %in% DR_freq$DDR)>=1)
  {sum(DR_freq %>% filter(DDR %in% t(ac)) %>% select(freq))}
  else Sdr<-0

  Sab<-if(sum(ac$acs %in% AB_freq$AA) >=1 & sum(ac$acs %in% AB_freq$BB) >=1)
  {sum(AB_freq %>% filter(AA %in% t(ac) & BB %in% t(ac)) %>% select(freq))}
  else Sab<-0

  Sadr<-if(sum(ac$acs %in% ADR_freq$AA) >=1 & sum(ac$acs %in% ADR_freq$DDR) >=1)
  {sum(ADR_freq %>% filter(AA %in% t(ac) & DDR %in% t(ac)) %>% select(freq))}
  else Sadr<-0

  Sbdr<-if(sum(ac$acs %in% BDR_freq$BB) >=1 & sum(ac$acs %in% BDR_freq$DDR) >=1)
  {sum(BDR_freq %>% filter(BB %in% t(ac) & DDR %in% t(ac)) %>% select(freq))}
  else Sbdr<-0

  Sabdr<-if(sum(ac$acs %in% ABDR_freq$AA && ac$acs %in% ABDR_freq$BB && ac$acs %in% ABDR_freq$DDR) >= 1)
  {sum(ABDR_freq %>% filter(AA %in% t(ac) & BB %in% t(ac) & DDR %in% t(ac)) %>% select(freq))}
  else Sabdr<-0

  l<-list(Sa,Sb,Sdr,Sab,Sadr,Sbdr,Sabdr)

  round(
    (1-(1-Sa-Sb-Sdr+Sab+Sadr+Sbdr-Sabdr)^2) * 100,
    2)
}


#' cPRA value for multiple candidates
#'
#' This function computes cPRA values with \code{cpra1} for multiple candidates with their own HLA sensitization
#' and given HLA allelic and haplotypic frequencies
#'
#' @param sensib a 2 columns data frame with patients identifications 'ID'
#' and respective HLA antibodies 'acs'
#' @param A_freq a 2 columns dataframe with HLA-A alelic frequencies.
#' A column (named 'AA') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param B_freq a 2 columns dataframe with HLA-B alelic frequencies.
#' A column (named 'BB') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param DR_freq a 2 columns dataframe with HLA-DR alelic frequencies.
#' A column (named 'DDR') with HLA-A alleles and a column (named 'freq') with respective frequencies
#' @param AB_freq a 3 columns dataframe with HLA-A-B haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles, another (named 'BB') with HLA-B alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param ADR_freq a 3 columns dataframe with HLA-A-DR haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles, another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param BDR_freq a 3 columns dataframe with HLA-B-DR haplotypic frequencies.
#' A column (named 'BB') with HLA-B alleles, another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @param ABDR_freq a 4 columns dataframe with HLA-A-B-DR haplotypic frequencies.
#' A column (named 'AA') with HLA-A alleles,
#' another column (named 'BB') with HLA-B alleles,
#' another (named 'DDR') with HLA-DR alleles
#' and a column (named 'freq') with respective haplotypic frequencies
#' @return a data frame with 2 columns, one with patients identification ('ids')
#' and the other  with percentual values corresponding to respective cPRA values ('cpra')
#' @author Bruno A Lima
#' @details This function returns percentages for cPRAs of multiple candidates
#' considering their respective HLA immunization and HLA allelic and haplotypic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors.
#' @export
#' @import tidyverse

cpras<-function(sensib = acHLA,
                A_freq = HLA_A,
                B_freq = HLA_B,
                DR_freq = HLA_DR,
                AB_freq = HLA_AB,
                ADR_freq = HLA_ADR,
                BDR_freq = HLA_BDR,
                ABDR_freq = HLA_ABDR){
  ids<-unique(sensib$ID)
  cpra<-NA
  for(i in 1:length(ids)){
    cpra[i]<-cpra1(id = ids[i], sensib = sensib,
                   A_freq = A_freq, B_freq = B_freq, DR_freq = DR_freq,
                   AB_freq = AB_freq, ADR_freq = ADR_freq, BDR_freq = BDR_freq,
                   ABDR_freq = ABDR_freq)
  }
  data.frame(ids,cpra)
}
