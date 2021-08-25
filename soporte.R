# Episodio 1 Conceptos básicos de R -----


## Objetos básicos ----

### Valores ----

a <- 1

v1 <- 2021

# Restar a v1 mi edad
v1 - 42

anio_nac <- v1 - 42

yo <- "Sergio"

#Create an object that has the value of number of pairs of human chromosomes
#Create an object that has a value of your favorite gene name
#Create an object that has this URL as its value: “ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/”
#Create an object that has the value of the number of chromosomes in a diploid human cell

human_chr_number <- 23
gene_name <- 'pten'
ensemble_url <- 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/'
human_diploid_chr_num <-  2 * human_chr_number


# Reasignar nombres de objetos
gene_name <- "tp53"

# Eliminar objeto
rm(gene_name)

gene_name


#### Ejercicio: Crear objetos y chequear su tipo con mode() ----
chromosome_name <- 'chr02'
od_600_value <- 0.47
chr_position <- '1001701'
spock <- TRUE
pilot <- Earhart  # No va a funcionar porque faltan las comillas

# Operaciones matemáticas 

# +	addition
# -	subtraction
# *	multiplication
# /	division
# ^ or **	exponentiation
# a%/%b	integer division (division where the remainder is discarded)
# a%%b	modulus (returns the remainder after division)

(1 + (5 ** 0.5))/2

human_chr_number * 2


### Vectores ----


# Create the SNP gene name vector

snp_genes <- c("OXTR", "ACTN3", "AR", "OPRM1")

# Check the mode, length, and structure of 'snp_genes'
mode(snp_genes)
length(snp_genes)
str(snp_genes)

#### Crear y subset de objetos

# Some interesting human SNPs
# while accuracy is important, typos in the data won't hurt you here

snps <- c('rs53576', 'rs1815739', 'rs6152', 'rs1799971')
snp_chromosomes <- c('3', '11', 'X', '6')
snp_positions <- c(8762685, 66560624, 67545785, 154039662)

# get the 3rd value in the snp_genes vector
snp_genes[3]

# get the 1st through 3rd value in the snp_genes vector
snp_genes[1:3]

# get the 1st, 3rd, and 4th value in the snp_genes vector
snp_genes[c(1, 3, 4)]

# get the 1st through the 3rd value, and 4th value in the snp_genes vector
# yes, this is a little silly in a vector of only 4 values.
snp_genes[c(1:3,4)]

# add the gene 'CYP1A1' and 'APOA5' to our list of snp genes
# this overwrites our existing vector
snp_genes <- c(snp_genes, "CYP1A1", "APOA5")


## Cargar datos ----

variants <- read.csv("data/SRR2584863_variants.csv", sep = ",")
vcf <- read.csv("data/combined_tidy_vcf.csv", sep = ",")

str(vcf)
