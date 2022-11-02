# =========================================================================== #
# Step 5 Normalization and Round 2 Filter
# --------------------------------------------------------------------------- #
# input
    # df.input.data.SP.norm
    # cutoff.sp

# Output:
    # df.input.data.Filt2           -> Save

# =========================================================================== #

                    # --------------------------------- #
                    #        Second Round Filter        #
                    # --------------------------------- #
# input: df.input.data.SP.norm
# output: df.input.data.Filter2
df.input.data.Filt2 = df.input.data.SP.norm
df.input.data.Filt2[, -c(1,2)] = apply(df.input.data.Filt2[, -c(1,2)],
                                       2,
                                       function(x) ifelse(x > cutoff.sp,
                                                          NA,
                                                          x))

# Check Point
# NA1 = apply(df.input.data.Filt1[,-c(1,2)], 2, function(x) sum(is.na(x)))
# NA2 = apply(df.input.data.Filt2[,-c(1,2)], 2, function(x) sum(is.na(x)))
# test = apply(df.input.data.SP.norm[,-c(1,2)], 2, function(x) length(which(x>cutoff.sp)))
# NA1 + test == NA2

write.csv(df.input.data.Filt2, file.path(dir.out.tbl, 'Data SP-inNorm.csv'))





