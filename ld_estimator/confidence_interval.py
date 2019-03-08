
def ld_confidence_interval(lsurface, loglike1, alpha=0.05):
    ''' estimate confidence interval for D'
    '''

    total = 0.0
    for i in range(len(lsurface)):
        lsurface[i] -= loglike1
        lsurface[i] = 10 ** lsurface[i]
        total += lsurface[i]

    low_i = 0
    summed = 0.0
    for i in range(100):
        summed += lsurface[i]
        if summed > alpha * total and (summed - lsurface[i]) < alpha * total:
            low_i = i - 1
            break

    high_i = 0
    summed = 0.0
    for i in range(99, 0, -1):
        summed += lsurface[i]
        if summed > alpha * total and (summed - lsurface[i]) < alpha * total:
            high_i = i + 1
            break

    if high_i > 100:
        high_i = 100

    return low_i / 100, high_i / 100
