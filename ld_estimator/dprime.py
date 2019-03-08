
def get_d(freqs):
    return freqs.aa * freqs.bb - freqs.ab * freqs.ba

def get_denominator(freqs):
    a = (freqs.aa + freqs.ba) * (freqs.ba + freqs.bb)
    b = (freqs.aa + freqs.ab) * (freqs.ab + freqs.bb)
    return min(a, b)
