library(comparison)
data(glass)
measurements = subset(glass, item %in% c("s1", "s2"))
control = makeCompItem(item ~ logKO + logCaO+ logFeO, data = measurements[1:6, ])

# calculate a compitem object representing the recovered item
# known to be from the same item (item 1)
recovered.1 = makeCompItem(item ~ logKO + logCaO + logFeO, data = measurements[7:12, ])

# calculate a compitem object representing the recovered item
# known to be from a different item (item 2)
recovered.2 = makeCompItem(item ~ logKO + logCaO + logFeO, data = measurements[19:24, ])


background = subset(glass, select = c(item, logKO, logCaO, logFeO))
Z = makeCompVar(background, 1)


# calculate the likelihood ratio for a known
# same source comparison 
# This value is 51.1424311 (9dp) in this version and the last version David wrote (1.0-4)
lr.1 = two.level.normal.LR(control, recovered.1, Z)

stopifnot(all.equal(lr.1, 51.1424311))

# calculate the likelihood ratio for a known
# different source comparison 
# This value is 0.0289990757 (9dp) in this version and the last version David wrote (1.0-4)
lr.2 = two.level.normal.LR(control, recovered.2, Z)

stopifnot(all.equal(lr.2, 0.0289990757))


# test the new calcLR routine - should be identical
lr.3 = calcLR(control, recovered.1, Z, "mvn")
stopifnot(all.equal(lr.3, 51.1424311))

lr.4 = calcLR(control, recovered.2, Z, "mvn")
stopifnot(all.equal(lr.4, 0.0289990757))


# calculate the likelihood ratio for a known
# same source comparison 
# 2020-08-01 Both this version and the previous version return 20.5896672
lr.5 = two.level.density.LR(control, recovered.1, Z)

stopifnot(all.equal(lr.5, 20.5896672))

# calculate the likelihood ratio for a known
# different source comparison
# 2020-08-01 Both this version and the previous version return 0.0116139204
lr.6 = two.level.density.LR(control, recovered.2, Z)

stopifnot(all.equal(lr.6, 0.0116139204))

# Test new calcLR - should be identical
lr.7 = calcLR(control, recovered.1, Z, "kde")

stopifnot(all.equal(lr.7, 20.5896672))

lr.8 = calcLR(control, recovered.2, Z, "kde")

stopifnot(all.equal(lr.8, 0.0116139204))
