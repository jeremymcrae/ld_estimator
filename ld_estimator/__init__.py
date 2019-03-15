from pkg_resources import get_distribution

name = 'ld_estimator'
__version__ = get_distribution(name).version

from ld_estimator.pairwise import pairwise_ld
from ld_estimator.region import region_ld
from ld_estimator.site import site_ld
