pdfConvolve2 <- function(p1, p2, n1, n2, k) {
    # Compute region of overlap
    beg <- max(0, k - n2)
    end <- min(k, n1)
    xs <- beg:end
    ys <- k - xs
    sum(dbinom(xs, n1, p1) * dbinom(ys, n2, p2))
}

count <- 10

p1s <- runif(count)
p2s <- runif(count)
n1s <- as.integer(runif(count, max=20) + 2)
n2s <- as.integer(runif(count, max=20) + 2)
ks <- as.integer((n1s + n2s) / 2)

results <- NULL

for (i in 1:count) {
    results <- rbind(results, pdfConvolve2(p1s[i], p2s[i], n1s[i], n2s[i], ks[i]))
}

df <- cbind(p1s, p2s, n1s, n2s, ks, results)

print(df)
