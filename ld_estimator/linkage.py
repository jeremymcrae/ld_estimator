''' class to store LD results
'''

class LD(object):
    ''' simple class to store linkage disequilibrium results
    '''
    def __init__(self, dprime, loglikelihood, r_squared, freqs, phase):
        self.dprime = dprime
        self.loglikelihood = loglikelihood
        self.r_squared = r_squared
        self.freqs = freqs
        self.phase = phase

    def __repr__(self):
        return f'LD({self.dprime}, {self.loglikelihood}, {self.r_squared}, ' \
            f'{self.freqs}, {self.phase})'
