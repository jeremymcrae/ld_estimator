''' package to calculate linkage disequilibrium
'''
from pkg_resources import get_distribution

from ld_estimator.pairwise import pairwise_ld
from ld_estimator.region import region_ld
from ld_estimator.site import site_ld

__version__ = get_distribution('ld_estimator').version
