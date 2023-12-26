library(immuneSIM)
number_of_sequences <- 69
hs_igh_sim <- immuneSIM(
                     number_of_seqs = number_of_sequences,
                     vdj_list = list_germline_genes_allele_01,
                     species = 'hs',
                     receptor = 'tr',
                     chain = 'b',
                     insertions_and_deletion_lengths = insertions_and_deletion_lengths_df,
                     name_repertoire = 'hs_trb_sim',
                     length_distribution_rand = length_dist_simulation,
                     random = FALSE,
                     vdj_noise = 0,
                     vdj_dropout = c(V=0, D=0, J=0),
                     ins_del_dropout = '',
                     equal_cc = FALSE,
                     freq_update_time = round(0.5 * number_of_sequences),
                     verbose = TRUE,
                     airr_compliant = TRUE,
                     shm.mode = 'none',
                     shm.prob = 15 / 350)
save(hs_igh_sim,file="hs_igh_sim")
write.csv(hs_igh_sim, "hs_igh_sim.csv")