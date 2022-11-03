##Parameters entry
haploXpop_filename <- "C:/Users/royau/Downloads/Galaxy97-[Table_Compute_on_data_96].tabular"
dist_filename <- "C:/Users/royau/Downloads/Galaxy99-[SNP_distance_matrix_on_data_98].tabular"
structure_filename <- "C:/Users/royau/Downloads/Galaxy114-[Regex_Find_And_Replace_on_data_113].tabular"
sep <- "\t"
dec <- "."

# on modifie ! 

# c'est bon ça !

##Steps
###Haplotypes x populations abundance matrix lol
tab_haploXpop <- read.table(file = haploXpop_filename, sep = sep, dec = dec, header = TRUE)
tab_haploXpop[is.na(tab_haploXpop)] <- 0
tab_haploXpop_c <- tab_haploXpop[,-1]

###Distance Matrix
# Open distance matrix file
tab_dist_haplo <- read.table(file = dist_filename, sep = sep, dec = dec)
dist_haplo <- as.dist(tab_dist_haplo[-1,-1])

# Create Euclidean distance matrix
euc_dist_haplo <- sqrt(dist_haplo)

# Verif if it is euclidean
ade4::is.euclid(euc_dist_haplo)
#TRUE

###Structure table (optional)
#tab_structure <- read.table(file = structure_filename, sep = sep, dec = dec)

###AMOVA TIME
# Perform the amova
Amv <- ade4::amova(tab_haploXpop_c, euc_dist_haplo)
Amv

# Test its significance through random test (permutation of matrix)
set.seed(1997)
Amvsignif <- ade4::randtest(Amv, nrepet = 999, alter = "two-sided")
plot(Amvsignif)
Amvsignif
