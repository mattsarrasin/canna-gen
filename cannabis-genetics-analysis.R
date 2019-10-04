library(rmarkdown)
library(adegenet)
library(ape)
library(vcfR)
library(poppr)
library(dplyr)
library(magrittr)
library(plotly)

#################################################
# setup
#################################################
vcfs <- c('Lynch.vcf.bgz', 'Sawler.vcf.bgz', 'Phylos.vcf.bgz', 'Phylos-pk.vcf.bgz')
vcfnames <- vcfs %>%
	strsplit(".vcf.bgz", fixed=TRUE) %>%
	unlist
srs2name <- read.table("sample_description_list-2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# make genlight objects from vcfs
gls <- lapply(vcfs, function(x) {read.vcfR(x) %>% vcfR2genlight()})
names(gls) <- vcfnames
# load in the phylos mds data separately
phylos.plink.mds <- read.table('phylos-plink.mds', header=TRUE, sep='\t', stringsAsFactors=FALSE)

# replace the srs IDs with variety names
for (gli in 1:length(gls)) {
	indNames(gls[[gli]]) <- srs2name$name[match(indNames(gls[[gli]]), srs2name$sample)]
}

# format the (pre-made) plink phylos data
phylos.plink.mds <- phylos.plink.mds %>%
	rename(sample=FID) %>%
	left_join(srs2name, on='sample')

# normalize function for setting an rgb color scheme in 3d plots
normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
}

#################################################
# run PCAs on allele frequencies
#################################################

gl.pcas <- lapply(gls, function(x) glPca(x, nf=3))

#################################################
# build 3d graph parameters for the PCAs
#################################################

# compute variance explained by axes and create colors for each pca dot
for (gli in 1:length(gl.pcas)) {
	gl.pcas[[gli]]$var.explain <- gl.pcas[[gli]]$eig[1:3]/sum(gl.pcas[[gli]]$eig) * 100 %>%
		round(3)
	gl.pcas[[gli]]$colors <- apply(gl.pcas[[gli]]$scores, 2, normalize) %>%
		apply(1, function(x) rgb(x[1], x[2], x[3])) %>%
		as.vector
}

# make the 3d plots for each dataset
pcas.3d <- lapply(gl.pcas, function(x) {
		plot_ly(
			x=x$scores[,1], y=x$scores[,2], z=x$scores[,3],
			type="scatter3d",
			mode="markers",
			text=~rownames(x$scores),
			marker=list(
				color=x$colors,
				size=4
			)
			) %>%
			layout(
				showlegend=FALSE,
				scene=list(
					xaxis=list(title=paste("PC1 (", x$var.explain[1], "%)")),
					yaxis=list(title=paste("PC2 (", x$var.explain[2], "%)")),
					zaxis=list(title=paste("PC3 (", x$var.explain[3], "%)"))
				)
			)
	}
)

phylos.plink.mds.3d <- plot_ly(x=phylos.plink.mds$C1, y=phylos.plink.mds$C2, z=phylos.plink.mds$C3,
	type="scatter3d",
	mode="markers",
	text=~phylos.plink.mds$name,
	marker=list(
		color=phylos.plink.mds.colors,
		size=5
	)
	) %>%
	layout(
		showlegend=FALSE,
		scene=list(
			xaxis=list(title="C1"),
			yaxis=list(title="C2"),
			zaxis=list(title="C3"))
	)

# visualize the 3d plots
lapply(pcas.3d, function(x) plotly_build)
phylos.plink.mds.3d 

#################################################
# genetic distance tree building
#################################################

# subset varieties of interest
sample.names <- c(
	'white_widdow',
	'chocolope',
	'cannatonic',
	'berry_haze',
	'rocky_mountain_blueberry',
	'og_kush_',
	'purple_kryptonite',
	'girl_scout_cookies',
	'hawaiian',
	'grape_kush',
	'grapefruit_diesel',
	'golden_goat',
	'blue_dream',
	'jack_herer',
	'jack_flash',
	'grape_ape',
	'grand_daddy_purp',
	'purple_urkle',
	'blue_cheese',
	'black_jack',
	'char_tango',
	'kool_aid_kush',
	'bubba_kush',
	'headband',
	'king_loule_cookies',
	'holy_grail',
	'skywalker_og',
	'kosher_kush',
	'sour_willie',
	'goast_train_haze',
	'the_sauce',
	'deadhead',
	'larry_og',
	'boss_hogg',
	'dog_walker',
	'old_school',
	'chinese_hemp',
	'carmagnola')
# find the indices of samples of interest
sample.indx <- sapply(sample.names, function(x) agrep(x, indNames(gls$Lynch))) %>% unlist %>% as.vector %>% unique
# subset the genlight object with those individuals
lynch.gli.select <- gls$Lynch[sample.indx]
# build genetic distance tree
lynch.tree <- aboot(lynch.gli.select, dist=bitwise.dist, tree="upgma", sample=100000, missing="ignore", showtree=FALSE, quiet=TRUE)
# plot the cladogram
plot.phylo(lynch.tree, cex=0.8, font=2, adj=0, no.margin=FALSE, show.tip.label=TRUE)
nodelabels(lynch.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8, font = 3, xpd = TRUE)
axisPhylo(3, cex=2, las=1, lty=3)
