"""HCHC-60 spiral inflector, pg5L geometry: build it, solve it, send the RFQ core beam through.

Besides spyral_inflector, PyPATools and py_electrodes this needs only the cyclotron
B-field map and the RFQ particle file (paths below). Outputs (STEP files, fields, bunch
results, exit metrics, geometry figures, openPMD hand-off files, report.md) go to OUT.
"""
import os

from spyral_inflector.tracking import build_geometry, run_final, mean_energy_mev

DECK = r"D:\MIT Dropbox\Daniel Winklehner\Projects\IsoDAR\60 MeV Cyclotron\Spiral_inflector"
BFIELD = os.path.join(DECK, "Fields", "HCHC-60_CentralBField_z-40to5cm_1mm.pickle")
PARTICLES = os.path.join(DECK, "Particles", "MIT RFQ Beamdynamics", "ext3_exit_core_as_txt.txt")
OUT = os.path.join(DECK, "Results", "example_pg5L")

# pg5L: 18 mm quad bore, 55/60 mm quads with one shared plate and 5 mm plate gaps, 19 mm entrance slot,
# 23 x 40 mm housing exit opening; the whole system shifted by dz = +3.50 mm (the chain's fringe result)
KNOBS = dict(gap=0.019, quad_bore=0.018, aper_hole=0.0175, slot_width=0.019, exit_opening=(0.023, 0.040),
             plate_gap=0.005, quad_z1=-0.264, quad_len=0.055, quad_len2=0.060, shared_plates=True)

if __name__ == "__main__":
    geo = build_geometry(os.path.join(OUT, "geometry"), os.path.join(OUT, "steps"), BFIELD, mean_energy_mev(PARTICLES),
                         knobs=KNOBS, fix_dz=0.0035, res=0.0025)          # optimizes the spiral voltage (~11.6 kV)
    run_final(geo["steps_dir"], geo["voltages_csv"], geo["state_pickle"], os.path.join(OUT, "run"),
              q1=6575, q2=-8450, rotate_quads=(16, 24), particles=PARTICLES, bfield=BFIELD, phi=90, sc=True)
