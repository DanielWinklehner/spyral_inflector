from spyral_inflector.geometry import SIAperture

if __name__ == "__main__":
    sia1 = SIAperture("Aperture 1", 0)
    sia1.create_geo_str(0.1, 0.003, 0.015, 0.006, hole_type="rectangle", h=0.005, load=True)
    sia1.show()

    sia2 = SIAperture("Aperture 2", 0)
    sia2.create_geo_str(0.1, 0.003, 0.015, 0.006, hole_type="ellipse", h=0.005, load=True)
    sia2.show()
