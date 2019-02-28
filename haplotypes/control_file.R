###Make Rqtl2 control file

library("qtl2")

for (i in 1:10){
    write_control_file(output_file =sprintf('Biogemma_c%d.json',i),
    crosstype = "riself16",geno_file=sprintf("../Biogemma_DHgenos/DH_geno_chr%d_121718.csv",i),
    founder_geno_file = sprintf("../Biogemma_foundergenos/Founder_genos_chr%d_121718.csv",i),
    gmap_file=sprintf("Biogemma_gmap_c%d.csv",i),pmap_file = sprintf("Biogemma_pmap_c%d.csv",i),covar_file="Biogemma_covar.csv",crossinfo_file = "Biogemma_cross_info.csv",
    geno_codes=c(A=1L,B=3L),
    alleles=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85"),
    sep=",",na.strings=c("NA"),comment.char="#",geno_transposed = FALSE,
    founder_geno_transposed = FALSE,description=sprintf("344 DH MAGIC Lines from 16 Biogemma Founders, Chromosome %d",i)
    )
}