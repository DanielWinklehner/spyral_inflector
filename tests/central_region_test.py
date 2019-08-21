import pytest
from spyral_inflector import *


class TestCentralRegion:

    def test_central_region_setup(self):

        h2p = IonSpecies("H2_1+", 0.02)
        h2p.calculate_from_energy_mev(0.4 / h2p.a())

        cr = CentralRegion(r_cr=[0.1, 0.3],
                           dee_voltage=70e3,
                           dee_opening_angle=35,
                           rf_phase=0.0,
                           ion=h2p)

        cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
                      vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

        cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

        gap = 0.06
        thickness = 0.025
        cl = 0.05

        my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                             r_init=0.05,
                             char_len=0.05,
                             gap=gap,
                             thickness=thickness)
        my_dee.initialize()
        my_dee.rotate(45, angle_unit="deg")
        my_second_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_second_dee.initialize()
        my_second_dee.rotate(90 + 45, angle_unit="deg")

        my_third_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                   r_init=0.05,
                                   char_len=0.05,
                                   gap=gap,
                                   thickness=thickness)
        my_third_dee.initialize(bottom_angle_offset=0, angle_unit="deg")
        my_third_dee.rotate(180 + 45, angle_unit="deg")

        my_fourth_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_fourth_dee.initialize()
        my_fourth_dee.rotate(270 + 45, angle_unit="deg")

        dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]
        cr._abstract_dees = dees

        for dee in dees:
            dee.make_transforms()
            for k in range(10):
                if k == 0 and dee is my_third_dee:
                    dee.next_bottom_segment(angle_offset=0, angle_unit='deg')
                else:
                    dee.next_bottom_segment()
                dee.next_top_segment()

    def test_central_region_scipy_orbits(self):

        h2p = IonSpecies("H2_1+", 0.02)
        h2p.calculate_from_energy_mev(0.4 / h2p.a())

        cr = CentralRegion(r_cr=[0.1, 0.3],
                           dee_voltage=70e3,
                           dee_opening_angle=35,
                           rf_phase=0.0,
                           ion=h2p)

        cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
                      vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

        cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

        res = orbit_finder(cr, energy_mev=0.3, verbose=False)

    def test_central_region_sector_field_calcs(self):

        df = TwoDeeField(left_voltage=0.0, right_voltage=70.0e3)

    def test_central_region_sector_tracking(self):

        h2p = IonSpecies("H2_1+", 0.02)
        h2p.calculate_from_energy_mev(0.4 / h2p.a())

        cr = CentralRegion(r_cr=[0.1, 0.3],
                           dee_voltage=70e3,
                           dee_opening_angle=35,
                           rf_phase=0.0,
                           ion=h2p)

        cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
                      vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

        cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

        gap = 0.06
        thickness = 0.025
        cl = 0.05

        my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                             r_init=0.05,
                             char_len=0.05,
                             gap=gap,
                             thickness=thickness)
        my_dee.initialize()
        my_dee.rotate(45, angle_unit="deg")
        my_second_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_second_dee.initialize()
        my_second_dee.rotate(90 + 45, angle_unit="deg")

        my_third_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                   r_init=0.05,
                                   char_len=0.05,
                                   gap=gap,
                                   thickness=thickness)
        my_third_dee.initialize(bottom_angle_offset=0, angle_unit="deg")
        my_third_dee.rotate(180 + 45, angle_unit="deg")

        my_fourth_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_fourth_dee.initialize()
        my_fourth_dee.rotate(270 + 45, angle_unit="deg")

        dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]
        cr._abstract_dees = dees

        for dee in dees:
            dee.make_transforms()
            for k in range(10):
                if k == 0 and dee is my_third_dee:
                    dee.next_bottom_segment(angle_offset=0, angle_unit='deg')
                else:
                    dee.next_bottom_segment()
                dee.next_top_segment()

        rstart = h2p.v_m_per_s() / (1.04 * h2p.q_over_m())
        r_init = np.array([rstart, 0.0, 0.0])
        v_init = np.array([0.0, h2p.v_m_per_s(), 0.0])

        simple_tracker(cr, r_start=r_init, v_start=v_init, dt=1e-11, nturns=1, phase=0.0)

    def test_central_region_full_field_calcs(self):

        h2p = IonSpecies("H2_1+", 0.02)
        h2p.calculate_from_energy_mev(0.4 / h2p.a())

        cr = CentralRegion(r_cr=[0.1, 0.3],
                           dee_voltage=70e3,
                           dee_opening_angle=35,
                           rf_phase=0.0,
                           ion=h2p)

        cr.initialize(xi=np.array([0.05, 0.0, 0.0]),
                      vi=h2p.v_m_per_s() * np.array([0.0, 1.0, 0.0]))

        cr.load_bfield(bfield=Field(dim=0, field={"x": 0.0, "y": 0.0, "z": -1.04}))

        gap = 0.06
        thickness = 0.025
        cl = 0.05

        my_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                             r_init=0.05,
                             char_len=0.05,
                             gap=gap,
                             thickness=thickness)
        my_dee.initialize()
        my_dee.rotate(45, angle_unit="deg")
        my_second_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_second_dee.initialize()
        my_second_dee.rotate(90 + 45, angle_unit="deg")

        my_third_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                   r_init=0.05,
                                   char_len=0.05,
                                   gap=gap,
                                   thickness=thickness)
        my_third_dee.initialize(bottom_angle_offset=0, angle_unit="deg")
        my_third_dee.rotate(180 + 45, angle_unit="deg")

        my_fourth_dee = AbstractDee(opening_angle=cr.dee_opening_angle,
                                    r_init=0.05,
                                    char_len=0.05,
                                    gap=gap,
                                    thickness=thickness)
        my_fourth_dee.initialize()
        my_fourth_dee.rotate(270 + 45, angle_unit="deg")

        dees = [my_dee, my_second_dee, my_third_dee, my_fourth_dee]
        cr._abstract_dees = dees

        for dee in dees:
            dee.make_transforms()
            for k in range(10):
                if k == 0 and dee is my_third_dee:
                    dee.next_bottom_segment(angle_offset=0, angle_unit='deg')
                else:
                    dee.next_bottom_segment()
                dee.next_top_segment()

        cr.make_dees(dees, n=4, voltage=cr.dee_voltage, gap=gap, thickness=thickness)
        cr.make_dummy_dees(gap=gap, thickness=thickness)
        cr.solve_bempp()
        calculate_potential(cr,
                            limits=((-0.25, 0.25), (-0.25, 0.25), (-0.01, 0.01)),
                            res=0.01,
                            domain_decomp=(4, 4, 4),
                            overlap=0)

        calculate_efield_bempp(cr)

    def test_central_region_modulated_track(self):
        pass
