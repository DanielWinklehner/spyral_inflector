from .settings import SettingsHandler
from scipy import constants as const
import os
import sys
import shutil
from PyPATools.colors import MyColors

# --- Set global variables from settings.txt file--- #
settings = SettingsHandler()
DEBUG = (settings["DEBUG"] == "True")
DECIMALS = int(settings["DECIMALS"])
TEMP_DIR = settings["TEMP_DIR"]
RELATIVISTIC = (settings["RELATIVISTIC"] == "True")
USE_MULTIPROC = True  # In case we are not using mpi or only using 1 processor, fall back on multiprocessing
ON_WINDOWS = "win" in sys.platform.lower()
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

if ON_WINDOWS:
    LOG_FONT = "Monospac821 BT"
    LOG_FONT_SIZE = 8
else:
    LOG_FONT = "Monospace"
    LOG_FONT_SIZE = 10

# Temporary directory for saving intermittent files
if os.path.exists(TEMP_DIR):
    shutil.rmtree(TEMP_DIR)
os.mkdir(TEMP_DIR)

# Other variables
COLORS = MyColors()

# --- Set global constants from scipy --- #
CLIGHT = const.speed_of_light
ECHARGE = const.elementary_charge
EPS0 = const.epsilon_0
AMU_MEV = const.value("atomic mass constant energy equivalent in MeV")
EMASS_MEV = const.value("electron mass energy equivalent in MeV")
PMASS_MEV = const.value('proton mass energy equivalent in MeV')

# --- Define factors to go to SI units --- #
kV = 1e3
MV = 1e6
GV = 1e9
kHz = 1e3
MHz = 1e6
km = 1e3
dm = 1e-1
cm = 1e-2
mm = 1e-3
pA = 1e-12
nA = 1e-9
uA = 1e-6
mA = 1e-3
