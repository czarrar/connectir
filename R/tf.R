tf.pspec <- function(vec) {
    library(fftw)
    N <- length(vec)
    x <- FFT(vec)/N
    x <- Mod(x[1:ceiling(N/2)])
    return(x)
}

tf.freqs <- function(N, Fs) {
    k <- 0:(N-1)     # Create vector from 0 to N-1
    T <- N/Fs        # Get the frequency interval
    freqs <- k/T     # Create the frequency range
    
    # only want 1st half of the FFT (positive portion)
    cutOff <- ceiling(N/2)
    freqs <- freqs[1:cutOff]
    
    return(freqs)
}
