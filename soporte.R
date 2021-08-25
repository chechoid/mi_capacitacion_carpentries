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
snp_genes <- c(snp_genes, "CYP1A1","APOA5")

# Seleccionar todo menos un elemento
snp_genes[-6]

# Modificar un elemento específico
snp_genes[7]<- "APOA5"
snp_genes

# Seleccionar lógicamente
snp_positions[snp_positions > 100000000] # 100 Millones

# Operator	Description
# <	less than
# <=	less than or equal to
# >	greater than
# >=	greater than or equal to
# ==	exactly equal to
# !=	not equal to
# !x	not x
# a | b	a or b
# a & b	a and b

snp_positions > 100000000

snp_positions[c(FALSE, FALSE, FALSE, TRUE)]

which(snp_positions > 100000000)


# Crear un objeto de límite
snp_marker_cutoff <- 100000000
snp_positions[snp_positions > snp_marker_cutoff]

snp_genes

# test to see if "ACTN3" or "APO5A" is in the snp_genes vector
# if you are looking for more than one value, you must pass this as a vector

c("ACTN3","APOA5") %in% snp_genes


# Episodio 2 factors and data frames ----

?read.csv()

## Cargar datos ----

variants <- read.csv("data/combined_tidy_vcf.csv")


## get summary statistics on a data frame

summary(variants)


## put the first three columns of variants into a new data frame called subset

subset <- data.frame(variants[,c(1:3,6)])

# Ver la estructura del data frame
str(subset)


## Factores ----

alt_alleles <- subset$ALT

# Creamos un vector con los elementos que sólo contengan un dígito
snps <- c(alt_alleles[alt_alleles=="A"],
          alt_alleles[alt_alleles=="T"],
          alt_alleles[alt_alleles=="G"],
          alt_alleles[alt_alleles=="C"])

str(snps)

plot(snps)

# Los errores surgen porque los elementos de snps no son de tipo factor

factor_snps <- factor(snps)

str(factor_snps)

summary(factor_snps)

plot(factor_snps)


# Si quiero generar un gráfico ordenado
table(as.factor(snps))

ordered_factor_snps <- factor(factor_snps, 
                              levels = names(sort(table(factor_snps),decreasing = T)))

plot(ordered_factor_snps)


# create a new data frame containing only observations from SRR2584863

SRR2584863_variants <- variants[variants$sample_id == "SRR2584863",]

# check the dimension of the data frame

dim(SRR2584863_variants)

summary(SRR2584863_variants)



str(SRR2584863_variants)


# Forzar (coerce) tipo de dato

snp_positions_2 <- c("8762685", "66560624", "67545785", "154039662")
typeof(snp_positions_2)

snp_positions_2[2]

as.numeric(snp_positions_2[2])

snp_positions_2 <- as.numeric(snp_positions_2)
typeof(snp_positions_2)


snp_chromosomes_2 <- c(3, 11, 'X', 6)
typeof(snp_chromosomes_2)

snp_chromosomes_2 <- as.numeric(snp_chromosomes_2)
snp_chromosomes_2


### Ejercicio ordenar factores ----
sorted_by_DP <- variants[order(variants$DP, decreasing = TRUE), ]

head(sorted_by_DP$DP)

# Cambiar el nombre a una columna
colnames(variants)[colnames(variants) == "sample_id"] <- "strain"

names(variants)


# Guardar un archivo
write.csv(SRR2584863_variants, file = "data/SRR2584863_variants.csv")

### Abrir Excel ---
# Desde el menú File -> Import Dataset


# Episodio 3 dplyr ----

install.packages("dplyr") # Instala el paquete

library(dplyr) # Carga el paquete

# Cómo me doy cuenta que un paquete está cargado

# Cargar datos 
variants <- read.csv("data/combined_tidy_vcf.csv")

## glimpse ----
glimpse(variants)
str(variants)

## select ----

select(variants, sample_id, REF, ALT, DP)

# Seleccionar todas menos una
select(variants, -CHROM)

variants_numeric <- select_if(variants, is.numeric)

glimpse(variants_numeric)

## filter ----

filter(variants, sample_id == "SRR2584863")

variants[variants$sample_id == "SRR2584863",]

## pipe ----

variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP) %>%
  head()


SRR2584863_variants <- variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)

### head tail slice ----

SRR2584863_variants %>% 
  head() # Muestra los primeros 6 registros

SRR2584863_variants %>% 
  tail() # Muestra los últimos 6 registros

SRR2584863_variants %>% 
  slice(10:15)

### Ejericicio

variants %>%
  filter(sample_id == "SRR2584863" & DP >= 10) %>%
  select(REF, ALT, POS)

unique(variants$sample_id)

# Filtrar por dos condiciones de la misma variable
variants %>%
  filter(sample_id %in% c("SRR2584863", "SRR2584866")
         & DP >= 10) %>%
  select(REF, ALT, POS)

## mutate ----

variants %>%
  mutate(POLPROB = 1 - (10 ^ -(QUAL/10))) %>%
  head()

# Ejercicio
variants %>%
  mutate(POLPROB = 1 - 10 ^ -(QUAL/10)) %>%
  select(sample_id, POS, QUAL, POLPROB) %>%
  head

## group_by & summarise ----

variants %>%
  group_by(sample_id) %>%
  summarize(n())

variants %>%
  group_by(sample_id) %>% 
  tally()

variants %>%
  count(sample_id)

variants %>%
  group_by(sample_id) %>%
  summarize(max(DP))


# Episodio 4 ggplot2 ---- 

# Instalar el paquete
install.packages("ggplot2")

library(ggplot2)

# Cargamos los datos
variants <- read.csv("data/combined_tidy_vcf.csv")

ggplot(variants)

ggplot(variants, aes(x = POS, y = DP))

# * `geom_point()` for scatter plots, dot plots, etc.
# * `geom_boxplot()` for, well, boxplots!
# * `geom_line()` for trend lines, time series, etc.

ggplot(variants, aes(x = POS, y = DP)) +
  geom_point()


# Podemos asignar un ggplot a un objeto
coverage_plot <- ggplot(data = variants, aes(x = POS, y = DP))

# Draw the plot
coverage_plot +
  geom_point()

## Proceso iterativo de diseño ----
coverage_plot +
  geom_point()

coverage_plot +
  geom_point(alpha = 0.5)


ggplot(data = variants, aes(x = POS, y = DP)) +
  geom_point(alpha = 0.5, color = "blue")


ggplot(data = variants, aes(x = POS, y = DP, color = sample_id)) +
  geom_point(alpha = 0.5)

ggplot(data = variants, aes(x = POS, y = DP, color = sample_id)) +
  geom_jitter(alpha = 0.5)


ggplot(data = variants, aes(x = POS, y = DP, color = sample_id)) +
  geom_jitter(alpha = 0.5) +
  labs(x = "Base Pair Position",
       y = "Read Depth (DP)")

## Faceting ----
ggplot(data = variants, aes(x = POS, y = MQ, color = sample_id)) +
  geom_point() +
  labs(x = "Base Pair Position",
       y = "Mapping Quality (MQ)") +
  facet_grid(. ~ sample_id)


ggplot(data = variants, aes(x = POS, y = MQ, color = sample_id)) +
  geom_point() +
  labs(x = "Base Pair Position",
       y = "Mapping Quality (MQ)") +
  facet_grid(sample_id ~ .)

ggplot(data = variants, aes(x = POS, y = MQ, color = sample_id)) +
  geom_point() +
  labs(x = "Base Pair Position",
       y = "Mapping Quality (MQ)") +
  facet_grid(sample_id ~ .) +
  theme_bw() +
  theme(panel.grid = element_blank())

## bar plots ----
ggplot(data = variants, aes(x = INDEL, fill = sample_id)) +
  geom_bar() +
  facet_grid(sample_id ~ .)

ggplot(data = variants, aes(x = INDEL, fill = sample_id)) +
  geom_bar() +
  facet_grid(sample_id ~ .) +
  theme_light()


# https://ggplot2.tidyverse.org/reference/ggtheme.ht

