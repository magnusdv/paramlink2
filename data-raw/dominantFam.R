## code to prepare `DATASET` dataset goes here

library(ibdsim2)
x = nuclearPed(4, sex = c(2,1,2,1))
x = addChildren(x, mother=3, nch = 3, sex = c(1,2,1))
x = addDaughter(x, 5)
x = addChildren(x, 6, nch = 1, sex = 1)
x = addChildren(x, 6, nch = 2, sex = c(2,2,2))
x = addSon(x, 10)
plot(x)

# Affection status
pattern = list(carriers = c(2,3,6,9,10,16,17,19),
               noncarriers = c(4,5,8,14))
plot(x, shaded = pattern$carriers)

# Decode map of Chr 1
chrom1 = loadMap(chrom = 1)

for(i in 1:1000) {
  s = ibdsim(x, map = chrom1, verbose = F, seed = i)[[1]]
  good = findPattern(s, pattern)
  if(nrow(good) == 1 && good[, 'length'] > 12) {
    message(i); print(good)
  }
}

i = 161
s = ibdsim(x, map = chrom1, verbose = F, seed = i)[[1]]
findPattern(s, pattern)
haploDraw(x, s, col = c(0,0,2,4), dist = 0.1, height = 2.5,
          pos = ifelse(1:pedsize(x) %in% unlist(pattern), 1, 0),
          margin = c(5,3,3,3), shaded = pattern$carriers)


# Attach ca 500 linked SNP markers on chromosome 1
pos = seq(1, 248, by = 0.5)   # marker positions in megabases
mlist = lapply(pos, function(i)
  marker(x, alleles = 1:2, chrom = 1, posMb = i))
x = setMarkers(x, mlist)

# Simulate onto the pattern
dominantFam = profileSimIBD(x, s, seed = 1234)

# Save
usethis::use_data(dominantFam, overwrite = TRUE)


#############
# LOD
aff = ifelse(1:pedsize(x) %in% pattern$carriers, 2, 1)
mod = diseaseModel("AD")
mld = merlinLod(dominantFam, aff, mod)
sld = lod(dominantFam, aff, mod)
plot(sld,ylim = c(-1,3.5), col = 4)
par(new = T)
plot(mld, ylim = c(-1,3.5))

