GR Energy Loss Calculator
=========================

## Requirements

- numpy
- pandas
- [pycatima](https://github.com/hrosiak/pycatima)

## Usage

```
./calc_eloss.py -b <magnetic field(mT)>
./calc_eloss.py --mag <magnetic field(mT)>
```

or

```
./calc_eloss.py -n <run number>
./calc_eloss.py --run <run number>
```

or

```
./calc_eloss.py -p <momentum(MeV/c)>
./calc_eloss.py --momentum <momentum(MeV/c)>
```

Example:

```
./calc_eloss.py -b 600
# B=600.00mT Brho=1800.00mTm p=539.63MeV/c
# Proton   T = 144.110MeV
| Material |     Ein |    Eout |   Eloss |  sigmaE |
|:--------:|--------:|--------:|--------:|--------:|
|   Pla1   | 143.413 | 140.989 |   2.442 |   0.237 |
|   Pla2   | 139.074 | 130.684 |   8.452 |   0.492 |
|   GAGG   | 130.296 |   0.000 | 131.244 |   0.498 |

# Deuteron T =  76.084MeV
| Material |     Ein |    Eout |   Eloss |  sigmaE |
|:--------:|--------:|--------:|--------:|--------:|
|   Pla1   |  74.144 |  67.054 |   7.138 |   0.230 |
|   Pla2   |  61.219 |  26.032 |  35.425 |   0.648 |
|   GAGG   |  23.499 |   0.000 |  23.658 |   0.653 |

# Triton   T =  51.365MeV
| Material |     Ein |    Eout |   Eloss |  sigmaE |
|:--------:|--------:|--------:|--------:|--------:|
|   Pla1   |  47.565 |  31.436 |  16.212 |   0.261 |
|   Pla2   |  12.863 |   0.000 |  12.929 |   0.362 |
|   GAGG   |   0.000 |   0.000 |   0.000 |   0.362 |
```
