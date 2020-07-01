import os
import platform
if platform.system() == "Linux" and os.environ.get('DISPLAY', '') == '':
    import matplotlib
    matplotlib.use('Agg')

 
from .fc_gaia import *
from .fc_skyquery import *
from .fcutils import *
