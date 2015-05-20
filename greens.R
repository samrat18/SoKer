suppressMessages(library(taRifx))
suppressMessages(library(bvpSolve))
suppressMessages(require(compiler))

path <- "/scratch/jishnu/kernel/greens/"
setwd(path)

#f = Sys.getenv('PBS_NODEFILE')
#nproc = strtoi(length(readLines(f)))
nproc=240
procid = strtoi(Sys.getenv('PBS_VNODENUM'))

derivfd <- function(f){
dim_f <- length(f)
f_fwd_shifted <- shift(f,1)
f_bwd_shifted <- shift(f,-1)
df <- 0.5*(f_fwd_shifted - f_bwd_shifted)
df[1] <- f[2] - f[1]
df[dim_f] <- f[dim_f] - f[dim_f-1]
return(df)
}

ellmax <- 450
numin <- 2.0
numax <- 4.5
ndiv <- 960
sampling <- 6
rsun <- 6.9598*10^10
diml <- rsun
dimc <- 10^6
dimrho <- 10^(-3)
dimp <- 10^9

g <- read.table("/home/jishnu/kernel/R_codes/m585q.4816")
dmpnu <- g[,"V3"]*10^(-6)
fwhm <- g[,"V5"]*10^(-6)
damping <- approxfun(dmpnu,fwhm)

model_s = read.table("/home/jishnu/kernel/R_codes/JCD_modelfull")
r <- model_s[, "V1"]
r <- rev(r)
dimr <- length(r)

cp <- model_s[, "V2"]
cp <- rev(cp)
c2 <- cp^2/dimc^2

rho <- model_s[, "V3"]
rho <- rev(rho)/dimrho

p <- model_s[, "V4"]
p <- rev(p)/dimp

#sampling
r <- r[seq(1, dimr, sampling)]
c2 <- c2[seq(1, dimr, sampling)]
rho <- rho[seq(1, dimr, sampling)]
p <- p[seq(1, dimr, sampling)]

#constants 
width <- (50. * 10^5)/rsun;
center_src <- 1. -(200. * 10^5)/rsun;
center_rcv <- 1. +(200. * 10^5)/rsun
nu <- seq (numin, numax, (numax-numin)/ndiv)*.001
#adding small damping to the frequency grid
for (count in 1:length(nu)){
nu[count] <- complex(real=nu[count],imaginary=(damping(nu[count])))*(diml/dimc)
}
#----------------------------------------
nu0 <- 3*0.001*(diml/dimc)
sigma <- 0.001*(diml/dimc)
dnu <- derivfd(nu)
fnu <- (-1)*(nu^2)*exp(-(nu-nu0)^2 / (2*sigma^2))

nr <- length(r)
scaling <- derivfd(r)
gravity <- (-1)*derivfd(p)/scaling/rho
goverc2 <- gravity/c2
hrho <- derivfd(log(rho))/scaling
N2 <- (-1)* gravity*( goverc2 + hrho)

goverc2 <- goverc2 + hrho*0.5
rho0 <- rho
rho <- 1.0

#rhoomega2 <- rho*(omega^2 - N2)
s <- exp((-1)*(r-center_src)^2 / (2.0*width^2))
rcv <- exp((-1)*(r-center_rcv)^2 / (2.0*width^2))

source <- function(x ,ell , omega ){
if (x < r[nr-5]) s <- approxfun (r ,s)(x)
else s <- 0.0
rhoomega2 <- rho*(omega^2 - N2)
sl2 <- ell*(ell+1)*c2/r^2
oneoverrhoc2 <- (sl2 /omega^2 -1)/ (rho*c2)
k1 <- approxfun(r, goverc2)(x)
k2 <- approxfun(r, oneoverrhoc2)(x)
k3 <- approxfun(r, rhoomega2)(x)
list(s, k1, k2, k3)
}

receiver <- function(x ,ell , omega ){
rcv <- approxfun (r ,rcv)(x)
rhoomega2 <- rho*(omega^2 - N2)
sl2 <- ell*(ell+1)*c2/r^2
oneoverrhoc2 <- (sl2 /omega^2 -1)/ (rho*c2)
k1 <- approxfun(r, goverc2)(x)
k2 <- approxfun(r, oneoverrhoc2)(x)
k3 <- approxfun(r, rhoomega2)(x)
list(rcv, k1, k2, k3)
}

diffeq_src <- function (x, y, pars){
params <- source(x,ell,omega)
s <- params[[1]]
k1 <- params[[2]]
k2 <- params[[3]]
k3 <- params[[4]]
list(c((k1-2/x)*y[1] +k2*y[2], k3*y[1] -k1*y[2] + s))
}

diffeq_rcv <- function (x, y, pars){
params <- receiver(x,ell,omega)
rcv <- params[[1]]
k1 <- params[[2]]
k2 <- params[[3]]
k3 <- params[[4]]
list(c((k1-2/x)*y[1] +k2*y[2], k3*y[1] -k1*y[2] + rcv))
}

#t <- proc.time()

for ( iomega in (procid*(ndiv/nproc)+1):(procid*(ndiv/nproc)+4)) {

#print (c('id-', procid,'iOmega',iomega,'nu[iOmega]',nu[iomega],'frequency-', nu[iomega]*dimc*1000/diml))
omega <- 2*pi*nu[iomega]
fname <- paste("omega-",iomega,".RData",sep="")

xisrc <- list()
psrc <- list()
xisrcatrcv <- list()
xisrcatsrc <- list()
psrcatrcv <- list()
xircv <- list()
prcv <- list()
xisrc_denkernel <- list()
xisrc_sskernel <- list()
psrc_denkernel <- list()
psrc_sskernel <- list()
xircv_sskernel <- list()
prcv_sskernel <- list()
tell <- proc.time()

for (ell in 0:ellmax) {

sol_src <- bvptwp( yini=c(0,NA), yend=c(NA,0),x = r, func=diffeq_src)
sol_rcv <- bvptwp( yini=c(0,NA), yend=c(NA,0),x = r, func=diffeq_rcv)

xisrc[[ell+1]] <- sol_src[,"1"]/sqrt(rho0)
psrc[[ell+1]] <- sol_src[,"2"]*sqrt(rho0)

xisrcatrcv[ell+1] <- approxfun(r,xisrc[[ell+1]])(center_rcv)
xisrcatsrc[ell+1] <- approxfun(r,xisrc[[ell+1]])(center_src)
psrcatrcv[ell+1] <- approxfun(r,psrc[[ell+1]])(center_rcv)
xircv [[ell+1]] <- sol_rcv[,"1"]/sqrt(rho0)
prcv[[ell+1]] <- sol_rcv[,"2"]*sqrt(rho0)

#terms required for kernels
xisrckernel <- ((2*pi)^3 *fnu[iomega]*dnu[iomega])*xisrc[[ell+1]]
xisrc_denkernel[[ell+1]] <-  xisrckernel * rho0
dxisrc <- derivfd( xisrc[[ell+1]] )
dxisrc_dr <- ((2*pi)^3 *fnu[iomega]* dnu[iomega] / scaling)*dxisrc
xisrc_sskernel [[ell+1]]<- dxisrc_dr * rho0 * c2
psrckernel <- ((2*pi)^3 *fnu[iomega]* dnu[iomega] )* psrc[[ell+1]]
psrc_denkernel[[ell+1]] <- (psrckernel / (r ^2 * rho0 ))/ omega^4
psrc_sskernel[[ell+1]] <- (psrckernel * c2 / r )/ omega^2
xircv_sskernel[[ell+1]] <- derivfd ( xircv[[ell+1]] ) / scaling
prcv_sskernel[[ell+1]] <- (prcv[[ell+1]] / (r * rho0) )/ omega^ 2
}
save(r, xisrc, xisrcatrcv, xisrcatsrc, psrc, psrcatrcv, xircv, prcv, xisrc_denkernel,xisrc_sskernel, psrc_denkernel, xircv_sskernel,psrc_sskernel, prcv_sskernel, file = fname)
#print (proc.time()-tell)
#plot(sol_src[,1], sol[,2], type = "l",lwd = 0.5, col = "red")
#plot(sol_rcv[,1], sol[,2], type = "l",lwd = 0.5, col = "green")
}
#print(proc.time()-t)

#stopCluster(cl)
