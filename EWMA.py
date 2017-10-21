def smoothed_z (z, L):
    """
    Computes the exponentially weighted moving average (with memory L) of input data z.
    Ported from original MATLAB function to Python by Z. Royston (Coleman lab):
    
    z is an input vector/1D array
    L is a time window in s
    """

    lam = 1-2/(L+1)

    smoothed = z[:] #'slice' the input array (i.e. copy)

    for j in range(1, len(z)):
        smoothed[j] = lam * smoothed[j-1] + (1-lam) * z[j]

    return smoothed

# Example
#j = [[5,5,6,7,8,1,2,9],[2,3,4,1,5,2,4,8]]
#l = 16
#
#s = smoothed_z(j, l)
#
#print(s)