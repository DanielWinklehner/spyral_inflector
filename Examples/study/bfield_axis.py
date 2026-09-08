"""B-field on the axis and the Larmor rotation across the quadrupole region."""
import os
import numpy as np
from PyPATools.field import Field

DECK = os.environ.get("SI_DECK") or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
bf = Field.from_file(os.path.join(DECK, "Fields", "HCHC-60_CentralBField_zflipped.pickle"))
z = np.linspace(-0.29, 0.02, 63)
b = bf(np.column_stack([np.zeros_like(z), np.zeros_like(z), z]))
print("   z [mm]    Bz [T]    Br(x=10mm) [T]")
bx = bf(np.column_stack([np.full_like(z, 0.01), np.zeros_like(z), z]))[:, 0]
for zi, bi, bxi in zip(z, b[:, 2], bx):
    print("  {:+7.1f}   {:+7.4f}   {:+8.5f}".format(1e3 * zi, bi, bxi))
# Larmor rotation angle of a 69.34 keV H2+ over the quadrupole region
brho = 0.0538  # T m at 69.3 keV H2+ (r_cyc 59.7 mm at 0.9 T)
sel = (z >= -0.275) & (z <= -0.11)
theta = np.trapezoid(b[sel, 2], z[sel]) / (2.0 * brho)
print("Larmor rotation from z = -275 to -110 mm: {:.1f} deg (int Bz dz = {:.4f} T m)".format(
    np.degrees(theta), np.trapezoid(b[sel, 2], z[sel])))
sel2 = (z >= -0.19) & (z <= -0.145)
print("across quad 2 (-190..-145): {:.1f} deg; quad 1 (-270..-225): {:.1f} deg".format(
    np.degrees(np.trapezoid(b[sel2, 2], z[sel2]) / (2 * brho)),
    np.degrees(np.trapezoid(b[(z >= -0.27) & (z <= -0.225), 2], z[(z >= -0.27) & (z <= -0.225)]) / (2 * brho))))
