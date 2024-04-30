library(terra)
filename <- system.file("ex/elev.tif", package="terra")
filename

r <- rast(filename)

plot(r, main = "DEM of some watershed I don't know")

hasValues(r)

r1 <- r2 
r2 <- r3 
r3 <- rast(nrow=10, ncol=10)

values(r1) <- runif(ncell(r1))
values(r2) <- runif(ncell(r2))
values(r3) <- runif(ncell(r3))

hasValues(r3)


stack <- c(r1,r2,r3)
#double brackets show layers individually

filename <- system.file("ex/logo.tif", package="terra")

logo_rast <- rast(filename)
plot(logo_rast)


plot(logo_rast[[1]])

r <- rast(ncol=10, nrow=10)
values(r) <- 1:ncell(r)

plot(r)

plot(r*100)

s <- r + 10
plot(s)
s <- sqrt(s)
s <- s * r + 5


values(r) <- runif(ncell(r))
r <- round(r)
plot(r)

x <- (r == 1)
plot(x)
plot(x*r)

s[r] <- -0.5
plot(s)
s[!r] <- 5
s[s == 5] <- 15
plot(s)


