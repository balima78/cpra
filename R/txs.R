#' Transplantability percentage for a single candidate, AB0 **compatible**
#'
#' This function computes an AB0 compatible transplantability percentual score for a single cnadidate
#' given their own HLA sensitization, it's AB0 blood group, and HLA and AB0 genetic frequencies
#'
#' @param id a numeric valic with patient's identification
#' @param sensib a 2 columns data frame with patient's identification 'ID'
#' and HLA antibodies 'acs'
#' @param grupos a 2 columns data frame with patients identification 'ID'
#' and AB0 blod type 'abo'
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
#' @param ABO_freq a 2 columns data frame with AB0 blood types frequencies;
#' A column (named 'grupo') with AB0 blood types,
#' another column (named 'freq') with respective frequencies.
#' @return a percentual value corresponding to compatible transplantability percentage
#' @author Bruno A Lima
#' @details This function returns a single percentage transplantability of a single candidate
#' considering AB0 compatible, candidate's own HLA immunization and HLA and AB0 genetic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors;
#' and AB0 frequencies are from Portuguese blood donors.
#' @export
#' @import tidyverse
#'
txComp1<-function(id = 1, sensib = acHLA, grupos = ABO,

                  A_freq = read.csv2("data/hla_A.csv"),
                  B_freq = read.csv2("data/hla_B.csv"),
                  DR_freq = read.csv2("data/hla_DRB1.csv"),
                  AB_freq = read.csv2("data/hla_AB.csv"),
                  ADR_freq = read.csv2("data/hla_ADR.csv"),
                  BDR_freq = read.csv2("data/hla_BDR.csv"),
                  ABDR_freq = read.csv2("data/hla_ABDR.csv"),

                  ABO_freq = read.csv2("data/AB0.csv")){

  if(filter(grupos,ID == id)$abo == "A")
  {probAB0 <- (1-filter(ABO_freq,grupo == "B")$freq)}
  else if(filter(grupos,ID == id)$abo == "B")
  {probAB0 <- (1-filter(ABO_freq,grupo == "A")$freq)}
  else if(filter(grupos,ID == id)$abo == "AB")
  {probAB0 <- 1}
  else probAB0 <- (filter(ABO_freq,grupo == "O")$freq)

  round(
    (1-cpra1(id = id ,sensib = sensib,
             A_freq = A_freq,
             B_freq = B_freq,
             DR_freq = DR_freq,
             AB_freq = AB_freq,
             ADR_freq = ADR_freq,
             BDR_freq = BDR_freq,
             ABDR_freq = ABDR_freq)/100) * probAB0^2 * 100,
    2)

}

#' Transplantability percentage for multiple candidates, AB0 **compatible**
#'
#' This function computes with \code{txComp,1} AB0 compatible transplantability percentual scores for multiple cnadidates
#' given their own HLA sensitization, it's AB0 blood group, and HLA and AB0 genetic frequencies
#'
#' @param sensib a 2 columns data frame with patients' identification 'ID'
#' and HLA antibodies 'acs'
#' @param grupos a 2 columns data frame with patients' identification 'ID'
#' and AB0 blod type 'abo'
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
#' @param ABO_freq a 2 columns data frame with AB0 blood types frequencies;
#' A column (named 'grupo') with AB0 blood types,
#' another column (named 'freq') with respective frequencies.
#' @return a 4 columns data frame with patients' identification ('ID'),
#' respective cPRA value computed withn \code{cpras} ('cPRA'), compatible transplantability percentage ('TxComp'),
#' and AB0 blood type ('abo')
#' @author Bruno A Lima
#' @details This function returns percentage transplantability of multiple candidates
#' considering AB0 compatible, candidates' own HLA immunization and HLA and AB0 genetic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors;
#' and AB0 frequencies are from Portuguese blood donors.
#' @export
#' @import tidyverse
#'
txComps<-function(sensib = acHLA, grupos = ABO,

                  A_freq = read.csv2("data/hla_A.csv"),
                  B_freq = read.csv2("data/hla_B.csv"),
                  DR_freq = read.csv2("data/hla_DRB1.csv"),
                  AB_freq = read.csv2("data/hla_AB.csv"),
                  ADR_freq = read.csv2("data/hla_ADR.csv"),
                  BDR_freq = read.csv2("data/hla_BDR.csv"),
                  ABDR_freq = read.csv2("data/hla_ABDR.csv"),

                  ABO_freq = read.csv2("data/AB0.csv")) {

  dados<-sensib %>% inner_join(grupos,by="ID")
  ids<-unique(dados$ID)
  cpra<-NA
  txb<-NA
  for(i in 1:length(ids)){

    cpra[i]<-cpra1(id = ids[i], sensib = sensib,
                   A_freq = A_freq, B_freq = B_freq, DR_freq = DR_freq,
                   AB_freq = AB_freq, ADR_freq = ADR_freq, BDR_freq = BDR_freq,
                   ABDR_freq = ABDR_freq)

    txb[i]<-txComp1(id = ids[i], sensib = sensib, grupos = grupos,
                    A_freq = A_freq, B_freq = B_freq, DR_freq = DR_freq,
                    AB_freq = AB_freq, ADR_freq = ADR_freq, BDR_freq = BDR_freq,
                    ABDR_freq = ABDR_freq,
                    ABO_freq = ABO_freq)
  }

  data.frame(ID = ids, cPRA = cpra,TxComp = txb) %>% left_join(grupos,by="ID")
}



#' Transplantability percentage for a single candidate AB0 **identical**
#'
#' This function computes an AB0 identical transplantability percentual score for a single cnadidate
#' given their own HLA sensitization, it's AB0 blood group, and HLA and AB0 genetic frequencies
#'
#' @param id a numeric value with patients identification
#' @param sensib a 2 columns data frame with patient's identification 'ID'
#' and HLA antibodies 'acs'
#' @param grupos a 2 columns data frame with patient's identification 'ID'
#' and AB0 blod type 'abo'
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
#' @param ABO_freq a 2 columns data frame with AB0 blood types frequencies;
#' A column (named 'grupo') with AB0 blood types,
#' another column (named 'freq') with respective frequencies.
#' @return a percentual value corresponding to compatible transplantability percentage
#' @author Bruno A Lima
#' @details This function returns a single percentage transplantability of a single candidate
#' considering AB0 identical, candidate's own HLA immunization and HLA and AB0 genetic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors;
#' and AB0 frequencies are from Portuguese blood donors.
#' @export
#' @import tidyverse
#'
txIdent1<-function(id = 1, sensib = acHLA, grupos = ABO,

                   A_freq = read.csv2("data/hla_A.csv"),
                   B_freq = read.csv2("data/hla_B.csv"),
                   DR_freq = read.csv2("data/hla_DRB1.csv"),
                   AB_freq = read.csv2("data/hla_AB.csv"),
                   ADR_freq = read.csv2("data/hla_ADR.csv"),
                   BDR_freq = read.csv2("data/hla_BDR.csv"),
                   ABDR_freq = read.csv2("data/hla_ABDR.csv"),

                   ABO_freq = read.csv2("data/AB0.csv")){


  if(filter(grupos,ID == id)$abo == "A")
  {probAB0 <- (1-filter(ABO_freq,grupo == "B")$freq)^2 - (filter(ABO_freq,grupo == "O")$freq)^2}
  else if(filter(grupos,ID == id)$abo == "B")
  {probAB0 <- (1-filter(ABO_freq,grupo == "A")$freq)^2 - (filter(ABO_freq,grupo == "O")$freq)^2}
  else if(filter(grupos,ID == id)$abo == "AB")
  {probAB0 <- 1-((1-filter(ABO_freq,grupo == "B")$freq)^2 + (1-filter(ABO_freq,grupo == "A")$freq)^2 - (filter(ABO_freq,grupo == "O")$freq)^2)}
  else probAB0 <- (filter(ABO_freq,grupo == "O")$freq)^2

  round(
    (1-cpra1(id = id ,sensib = sensib,
             A_freq = A_freq,
             B_freq = B_freq,
             DR_freq = DR_freq,
             AB_freq = AB_freq,
             ADR_freq = ADR_freq,
             BDR_freq = BDR_freq,
             ABDR_freq = ABDR_freq)/100) * probAB0 * 100,
    2)

}


#' Transplantability percentage for multiple candidates, AB0 **identical**
#'
#' This function computes with \code{txIdent1} AB0 identical transplantability percentual scores for multiple cnadidates
#' given their own HLA sensitization, it's AB0 blood group, and HLA and AB0 genetic frequencies
#'
#' @param sensib a 2 columns data frame with patients' identification 'ID'
#' and HLA antibodies 'acs'
#' @param grupos a 2 columns data frame with patients' identification 'ID'
#' and AB0 blod type 'abo'
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
#' @param ABO_freq a 2 columns data frame with AB0 blood types frequencies;
#' A column (named 'grupo') with AB0 blood types,
#' another column (named 'freq') with respective frequencies.
#' @return a 4 columns data frame with patients' identification ('ID'),
#' respective cPRA value computed withn \code{cpras} ('cPRA'), compatible transplantability percentage ('TxComp'),
#' and AB0 blood type ('abo')
#' @author Bruno A Lima
#' @details This function returns percentage transplantability of multiple candidates
#' considering AB0 compatible, candidates' own HLA immunization and HLA and AB0 genetic frequencies.
#' By default HLA frequencies used are those from 37.993 Portuguese bone marrow voluntary donors;
#' and AB0 frequencies are from Portuguese blood donors.
#' @export
#' @import tidyverse
#'
txIdents<-function(sensib = acHLA, grupos = ABO,

                   A_freq = read.csv2("data/hla_A.csv"),
                   B_freq = read.csv2("data/hla_B.csv"),
                   DR_freq = read.csv2("data/hla_DRB1.csv"),
                   AB_freq = read.csv2("data/hla_AB.csv"),
                   ADR_freq = read.csv2("data/hla_ADR.csv"),
                   BDR_freq = read.csv2("data/hla_BDR.csv"),
                   ABDR_freq = read.csv2("data/hla_ABDR.csv"),

                   ABO_freq = read.csv2("data/AB0.csv")){

  dados<-sensib %>% inner_join(grupos,by="ID")
  ids<-unique(dados$ID)
  cpra<-NA
  txb<-NA
  for(i in 1:length(ids)){

    cpra[i]<-cpra1(id = ids[i], sensib = sensib,
                   A_freq = A_freq, B_freq = B_freq, DR_freq = DR_freq,
                   AB_freq = AB_freq, ADR_freq = ADR_freq, BDR_freq = BDR_freq,
                   ABDR_freq = ABDR_freq)

    txb[i]<-txIdent1(id = ids[i], sensib = sensib, grupos = grupos,
                     A_freq = A_freq, B_freq = B_freq, DR_freq = DR_freq,
                     AB_freq = AB_freq, ADR_freq = ADR_freq, BDR_freq = BDR_freq,
                     ABDR_freq = ABDR_freq,
                     ABO_freq = ABO_freq)
  }

  data.frame(ID = ids, cPRA = cpra,TxIdent = txb) %>% left_join(grupos,by="ID")

}

