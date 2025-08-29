"""
TuMagPipe package: tools for image destretching and reduction pipeline.
"""

__version__ = "0.1.0"

from . import (image_handler,master_dark,master_flatfield,image_filtering,
               destretch,demodulation,fits_files_handling,prefilter_removal,config,
               zernike,math_func2,shift_func,get_rotation,
               alignment)
# from . import (alignment,Calibration_finder,Check_ocs,config,demodulation,destretch,image_filtering,image_handler,master_dark,master_flat)

import logging
from .logutils import log_memory


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)
logger.info("mypackage initialized.")

