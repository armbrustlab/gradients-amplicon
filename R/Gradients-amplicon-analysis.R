library("biomformat")
library("vegan")
library("dplyr")
library("ggplot2")

setwd('Documents/gradients-amplicon')

# Load in 16S ASVs and convert to matrix
proks = read_biom('data/feature-table-16S.biom')
proks = as(biom_data(proks), "matrix")

# Same as above for 18S
euks = read_biom('data/feature-table-18S.biom')
euks = as(biom_data(euks), "matrix")

# Load sample metadata
metadata = read.csv('data/Gradients-amplicon-metadata.csv')

# Trim excel relic extra 'X' column
metadata$X = NULL

#----------------------------------------------------------
#
#
# Analysis of the 'surface' prokaryotic and eukaryotic communities
# The provided metadata is a little wonky, but I define surface
# as samples taken at 15 m or with the depth label 'surface'.
#
#
#----------------------------------------------------------

# Deal with 16S first
proks = as.data.frame(proks)
euks = as.data.frame(euks)

surf_metadata = metadata[which(metadata$depth=='surface' | metadata$depth == 15), ]

samples = surf_metadata$sample.id

# Now pull out the corresponding samples
surf_proks = proks[ names(proks)[names(proks) %in% samples] ]

# Now separate the proks into particle and free-living
surf_metadata_pa = surf_metadata[which(surf_metadata$filter==3), ]
surf_metadata_fl = surf_metadata[which(surf_metadata$filter==0.2), ]

samples_pa = surf_metadata_pa$sample.id
samples_fl = surf_metadata_fl$sample.id

surf_proks_fl = surf_proks[ names(surf_proks)[names(surf_proks) %in% samples_fl] ]
surf_proks_pa = surf_proks[ names(surf_proks)[names(surf_proks) %in% samples_pa] ]

# Remove rows containing all zeros. The samples were processed together, so there
# is a chance that some ASVs were detected in the PA, but not the FL

surf_proks_fl = surf_proks_fl[apply(surf_proks_fl[,-1], 1, function(x) !all(x==0)),]
# This almost drops the number of ASVs in the FL to a 1/3

surf_proks_pa = surf_proks_pa[apply(surf_proks_pa[,-1], 1, function(x) !all(x==0)),]
# More than halves the ASV count in the PA.

rarecurve(t(surf_proks_fl), step = 50, main='Gradients Free-Living Prokaryotes')

rarecurve(t(surf_proks_pa), step = 50, main='Gradients Particle-attached Prokaryotes')
# Looks like there are a couple particle attached samples that might be problematic. 
# Check what happens to the agreement in species number if you rarefy.

# Free-living comparison
S <- specnumber(t(surf_proks_fl)) # observed number of species
Srare <- rarefy(t(surf_proks_fl), 18000)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients FL')
abline(0, 1,  col="red", lwd=3, lty=2)

# PA comparison
S <- specnumber(t(surf_proks_pa)) # observed number of species
Srare <- rarefy(t(surf_proks_pa), 9000)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients PA')
abline(0, 1,  col="red", lwd=3, lty=2)

# Both deviate slightly at high species richness, but not dramatically

fl = rrarefy(t(surf_proks_fl), 18000)
pa = rrarefy(t(surf_proks_pa), 9000)

fl_shannon = diversity(fl)
pa_shannon = diversity(pa)
  
surf_metadata_fl$S.index = fl_shannon  
surf_metadata_pa$S.index = pa_shannon
  
f = ggplot(surf_metadata_fl, aes(latitude, S.index))
f + geom_point(aes(shape=cruise, color=cruise), size=3) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Latitude") + ylab("Shannon Index")


f = ggplot(surf_metadata_pa, aes(latitude, S.index))
f + geom_point(aes(shape=cruise, color=cruise), size=3) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Latitude") + ylab("Shannon Index")



# 18S ASVs


# Now pull out the corresponding samples
surf_euks = euks[ names(euks)[names(euks) %in% samples] ]

# Now separate the euks into large and small
surf_metadata_l = surf_metadata[which(surf_metadata$filter==3), ]
surf_metadata_s = surf_metadata[which(surf_metadata$filter==0.2), ]

samples_l = surf_metadata_l$sample.id
samples_s = surf_metadata_s$sample.id

surf_euks_l = surf_euks[ names(surf_euks)[names(surf_euks) %in% samples_l] ]
surf_euks_s = surf_euks[ names(surf_euks)[names(surf_euks) %in% samples_s] ]

# Remove rows containing all zeros. The samples were processed together, so there
# is a chance that some ASVs were detected in the PA, but not the FL

surf_euks_l = surf_euks_l[apply(surf_euks_l[,-1], 1, function(x) !all(x==0)),]
# This almost drops the number of ASVs in the FL to a 1/3

surf_euks_s = surf_euks_s[apply(surf_euks_s[,-1], 1, function(x) !all(x==0)),]
# More than halves the ASV count in the PA.

rarecurve(t(surf_euks_l), step = 50, main='Gradients Large fraction eukaryotes')

rarecurve(t(surf_euks_s), step = 50, main='Gradients Small fraction eukaryotes')
# Looks like there are a couple particle attached samples that might be problematic. 
# Check what happens to the agreement in species number if you rarefy.

# Small comparison
S <- specnumber(t(surf_euks_s)) # observed number of species
Srare <- rarefy(t(surf_euks_s), 7900)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients Small-Euks')
abline(0, 1,  col="red", lwd=3, lty=2)

# Large comparison
S <- specnumber(t(surf_euks_l)) # observed number of species
Srare <- rarefy(t(surf_euks_l), 2600)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients Large-Euks')
abline(0, 1,  col="red", lwd=3, lty=2)

# So our first problem arises. We will have to filter out some eukaryotic samples to get better rarefaction depth.
# We need to sort out which samples have very low numbers of reads. 
# This can be done from the Qiime2 run metadata.

# First list of sample ids are those that have less than 5,000 reads as input
# The extended list has <10,000

flagged_samples = c(63693, 63837, 64041, 68584, 68585, 68728)

# Not all Qiime2 samples are in this list. Find the intersect!
intersect(samples, flagged_samples)

# Only two: 63837 64041
# Check the extended list
flagged_samples_ext = c(63693, 63837, 64041, 68584, 68585, 68728, 68486,
                        68529, 68534, 68620, 68887)

intersect(samples, flagged_samples_ext)
# 63837 64041 68486 68529 68887
samples_to_remove = c(63837, 64041, 68486, 68529, 68887)

# Let's remove the extended set and check the rarefaction information out again.

samples_filt = samples[!samples %in% samples_to_remove]

# Oh yeahhhhh, let's do it agaiinnnnn. Dooo youuu, youuuuuu feeeel like I doooo.
surf_euks_f = euks[ names(euks)[names(euks) %in% samples_filt] ]

surf_metadata_f = surf_metadata[which(surf_metadata$sample.id %in% samples_filt),]

# Now separate the euks into large and small
surf_metadata_l = surf_metadata_f[which(surf_metadata_f$filter==3), ]
surf_metadata_s = surf_metadata_f[which(surf_metadata_f$filter==0.2), ]

samples_l = surf_metadata_l$sample.id
samples_s = surf_metadata_s$sample.id

surf_euks_l = surf_euks_f[ names(surf_euks_f)[names(surf_euks_f) %in% samples_l] ]
surf_euks_s = surf_euks_f[ names(surf_euks_f)[names(surf_euks_f) %in% samples_s] ]

surf_euks_l = surf_euks_l[apply(surf_euks_l[,-1], 1, function(x) !all(x==0)),]
surf_euks_s = surf_euks_s[apply(surf_euks_s[,-1], 1, function(x) !all(x==0)),]



rarecurve(t(surf_euks_l), step = 50, main='Gradients Large fraction eukaryotes - post filter')

S <- specnumber(t(surf_euks_l)) # observed number of species
Srare <- rarefy(t(surf_euks_l), 25000)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients Large-Euks',
     xlim = c(0,1500), ylim = c(0,1500))
abline(0, 1,  col="red", lwd=3, lty=2)


rarecurve(t(surf_euks_s), step = 50, main='Gradients Small fraction eukaryotes - post filter')

S <- specnumber(t(surf_euks_s)) # observed number of species
Srare <- rarefy(t(surf_euks_s), 20000)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", pch=16, main = 'Gradients Small-Euks',
     xlim = c(0,1500), ylim = c(0,1500))
abline(0, 1,  col="red", lwd=3, lty=2)

# Unfortunately, it looks like our 18S sequencing efforts came up a bit short.
# Except for G3. We can try to use iNEXT to determine what sample depth we need for
# xx% completeness

################################################################################
library(iNEXT)

# Don't run. This takes forever.
out = iNEXT(t(surf_euks_l), q = c(0), datatype = "abundance")
ggiNEXT(out, type = 2)
################################################################################

# Alternative is to estimate the number of 'species' present and use the bootstrap
# index to exclude samples that don't cut it.





