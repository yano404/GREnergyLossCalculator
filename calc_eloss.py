#!/bin/env python3

"""
GR Energy Loss Calculator
"""

import argparse
import numpy as np
import pycatima as atima

# Constants
um2cm = 1e-4
mm2cm = 1e-1
sqrt2 = np.sqrt(2)
rho = 3.0 # meters

# Layers
exit_window = atima.Layers()
vdc = atima.Layers()
sci_pla1 = atima.Layers()
al_plate = atima.Layers()
sci_pla2 = atima.Layers()
sci_gagg = atima.Layers()
mat = atima.Layers()
idx_pla1 = 0
idx_pla2 = 0
idx_gagg = 0

# Materials
kapton = atima.get_material(atima.material.Kapton)
air = atima.get_material(atima.material.Air)
vdc_film = atima.get_material(atima.material.Kapton)
vdc_cathode = atima.get_material(atima.material.Kapton)
# Ar/Iso-Buthane = 9/1
vdc_gas = atima.Material(
        [[0,1,30],
         [0,6,12],
         [0,18,7]],
        density=0.001865 #g/cm^3
        )
black_sheet = atima.get_material(atima.material.CH2)
pla = atima.get_material(atima.material.BC_400)
al = atima.Material([[0,13,1]], density=2.7)
gagg = atima.Material(
        [[0,64,3],
         [0,13,2],
         [0,31,3],
         [0,8,12],
         [0,58,1]],
        density=6.63 #g/cm^3
        )

# Construct Layers

# GR exit window
# Kapton 50um
kapton.thickness_cm(50.0*um2cm*sqrt2)
exit_window.add(kapton)

# VDC
vdc_film.thickness_cm(25.0*um2cm*sqrt2)
vdc_cathode.thickness_cm(6.0*um2cm*sqrt2)
vdc_gas.thickness_cm(10.0*mm2cm*sqrt2)
vdc.add(vdc_film) # Entrance window film
vdc.add(vdc_gas)
vdc.add(vdc_cathode) # 1st cathode film
vdc.add(vdc_gas) # Gas
# x anode wire
vdc.add(vdc_gas)
vdc.add(vdc_cathode) # 2nd cathode film
vdc.add(vdc_gas)
# u anode film
vdc.add(vdc_gas)
vdc.add(vdc_cathode) # 3rd cathode film
vdc.add(vdc_gas)
vdc.add(vdc_film) # Exit window film

# Plastic1 3mmt
black_sheet.thickness_cm(100.0*um2cm*sqrt2)
pla.thickness_cm(3.0*mm2cm*sqrt2)
sci_pla1.add(black_sheet)
sci_pla1.add(pla)
sci_pla1.add(black_sheet)
# Al Plate 1mmt
al.thickness_cm(1.0*mm2cm*sqrt2)
al_plate.add(al)
# Plastic2 10mmt
pla.thickness_cm(10.0*mm2cm*sqrt2)
sci_pla2.add(black_sheet)
sci_pla2.add(pla)
sci_pla2.add(black_sheet)
# GAGG
black_sheet.thickness_cm(200.0*um2cm)
gagg.thickness_cm(35.0*mm2cm)
sci_gagg.add(black_sheet)
sci_gagg.add(gagg)

# Construct layers
mat.add_layers(exit_window) # GR exit window
air.thickness_cm((10.0+12.0)*mm2cm*sqrt2) # VDC window frame 12mmt
mat.add(air)
mat.add_layers(vdc) # VDC1
air.thickness_cm((166.0+12.0+12.0)*mm2cm*sqrt2) # VDC window frame 12mmt
mat.add(air)
mat.add_layers(vdc) # VDC2
air.thickness_cm((12.0+107.25)*mm2cm*sqrt2) # VDC window frame 12mmt
mat.add(air)
mat.add_layers(sci_pla1) # Pla1 (3mmt)
idx_pla1 = mat.num() - 2 # Index for pla1
air.thickness_cm(35.25*mm2cm*sqrt2)
mat.add(air)
mat.add_layers(al_plate) # Al plate 1mmt
air.thickness_cm(45.25*mm2cm*sqrt2)
mat.add(air)
mat.add_layers(sci_pla2) # Pla2 (10mmt)
idx_pla2 = mat.num() - 2 # Index for pla1
air.thickness_cm(294.0*mm2cm)
mat.add(air)
mat.add_layers(sci_gagg) # GAGG
idx_gagg = mat.num() - 1 # Index for pla1

mp = 1.007_276_467 # u
md = 2.013_553_213 # u
mt = 3.015_500_716 # u
qp = 1.0 # e
qd = 1.0 # e
qt = 1.0 # e
Ap = 1.0
Ad = 2.0
At = 3.0


def Brho2p(Brho):
    # Brho to momentum (q=1)
    # Brho in mT*m units
    # q in e units
    # return in MeV/c
    p = 0.299792458 * Brho
    return p


def p2T(p, mu):
    # Momentum to Kinetic Energy
    # p in MeV/c
    # mu in amu
    amu = 931.49410242 # MeV/c^2
    m = mu*amu # MeV/c^2
    E = np.sqrt(np.power(p,2)+np.power(m,2))
    T = E - m
    return T


def main():
    # Parser
    parser = argparse.ArgumentParser(
            description='Energy loss calculator for GR')
    parser.add_argument(
            'b',
            type=float,
            help='Magnetic field strength(mT)')
    args = parser.parse_args()
    # Rigidity
    b = args.b
    # Momentum
    k = Brho2p(b*rho)
    kp = k * qp
    kd = k * qd
    kt = k * qt
    # Energy
    Tp = p2T(kp, mp)
    Td = p2T(kd, md)
    Tt = p2T(kt, mt)

    # Projectile
    p = atima.Projectile(mp, qp, qp, Tp/Ap)
    d = atima.Projectile(md, qd, qd, Td/Ad)
    t = atima.Projectile(mt, qt, qt, Tt/At)

    # Calculate energy loss
    resp = atima.calculate_layers(p, mat)
    resd = atima.calculate_layers(d, mat)
    rest = atima.calculate_layers(t, mat)

    resp_pla1 = resp.results[idx_pla1]
    resd_pla1 = resd.results[idx_pla1]
    rest_pla1 = rest.results[idx_pla1]
    resp_pla2 = resp.results[idx_pla2]
    resd_pla2 = resd.results[idx_pla2]
    rest_pla2 = rest.results[idx_pla2]
    resp_gagg = resp.results[idx_gagg]
    resd_gagg = resd.results[idx_gagg]
    rest_gagg = rest.results[idx_gagg]
    
    # Print results
    print("# Proton")
    print( "| Material |     Ein |    Eout |   Eloss |")
    print( "|:--------:|--------:|--------:|--------:|")
    print(f"|   Pla1   | {resp_pla1.Ein*Ap:7.3f} | {resp_pla1.Eout*Ap:7.3f} | {resp_pla1.Eloss:7.3f} |")
    print(f"|   Pla2   | {resp_pla2.Ein*Ap:7.3f} | {resp_pla2.Eout*Ap:7.3f} | {resp_pla2.Eloss:7.3f} |")
    print(f"|   GAGG   | {resp_gagg.Ein*Ap:7.3f} | {resp_gagg.Eout*Ap:7.3f} | {resp_gagg.Eloss:7.3f} |")
    print("")

    print("# Deuteron")
    print( "| Material |     Ein |    Eout |   Eloss |")
    print( "|:--------:|--------:|--------:|--------:|")
    print(f"|   Pla1   | {resd_pla1.Ein*Ad:7.3f} | {resd_pla1.Eout*Ad:7.3f} | {resd_pla1.Eloss:7.3f} |")
    print(f"|   Pla2   | {resd_pla2.Ein*Ad:7.3f} | {resd_pla2.Eout*Ad:7.3f} | {resd_pla2.Eloss:7.3f} |")
    print(f"|   GAGG   | {resd_gagg.Ein*Ad:7.3f} | {resd_gagg.Eout*Ad:7.3f} | {resd_gagg.Eloss:7.3f} |")
    print("")

    print("# Triton")
    print( "| Material |     Ein |    Eout |   Eloss |")
    print( "|:--------:|--------:|--------:|--------:|")
    print(f"|   Pla1   | {rest_pla1.Ein*At:7.3f} | {rest_pla1.Eout*At:7.3f} | {rest_pla1.Eloss:7.3f} |")
    print(f"|   Pla2   | {rest_pla2.Ein*At:7.3f} | {rest_pla2.Eout*At:7.3f} | {rest_pla2.Eloss:7.3f} |")
    print(f"|   GAGG   | {rest_gagg.Ein*At:7.3f} | {rest_gagg.Eout*At:7.3f} | {rest_gagg.Eloss:7.3f} |")


if __name__ == '__main__':
    main()

