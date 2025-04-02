import numpy as np
import sympy as sp
import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI
# from Numerical_simulation import *

class const_var:
    def __init__(self):
        # Constants
        self.D_w_f = 9 * 10 ** -7  # Effective moisture diffusivity in membrane [m2/s]
        self.k_fiber = 0.17 / 1000  # Thermal conductivity of fiber [kW/mK]
        self.k_air = 0.0263  # Placeholder for thermal conductivity function [W/mK]
        self.k_water = 0.612  # Placeholder for thermal conductivity function [W/mK]
        self.nu_a = 1.52e-5  # Placeholder for kinematic viscosity function [m2/s]
        self.D_va = 2.2e-5  # Moisture diffusivity [m2/s]
        self.Pr_a = 0.71  # Placeholder for Prandtl number function [-]
        self.Pr_w = 7.008 # Prandtl number of water [-]
        self.nu_w_in = 1.007e-6  # Placeholder for kinematic viscosity function [m2/s]

class water:
    def __init__(self, temp, flowrate, rho, cp):
        self.T = temp  # water temperature [C]
        self.V_dot = flowrate  # volume flow rate of water [lpm]
        self.cp = cp
        self.rho = rho
        self.Molar_mass = 0.018  # molecule weigh of vapor [kg/mol]
        self.V_dot_cms = self.V_dot / (1000 * 60)  # [m3/s]
        self.m_dot = self.rho * self.V_dot_cms  # [kg/s]
        self.c = self.m_dot * self.cp  # Heat capacity of water [kW/K]
        self.enthalpy = CP.PropsSI('H', 'T', self.T + 273.15, 'Q', 0, 'Water') / 1000
        self.velocity = None


class air:
    def __init__(self, temp, rh, flowrate, rho, cp):
        self.T = temp  # Air dry-bulb temperature [C]
        self.RH = rh  # Air relative humidity [%]
        self.V_dot = flowrate  # volume flow rate of air [m3/h]
        self.rho = rho  # Placeholder for density function [kg/m3]
        self.cp = cp  # Placeholder for specific heat capacity function [J/kgK]
        self.Molar_mass = 0.029  # molecule weigh of dry air [kg/mol]
        self.w = HAPropsSI('W', 'T', self.T + 273.15, 'R', self.RH / 100, 'P', 101325)
        self.m_dot = self.rho * self.V_dot / 3600  # mass flow rate of air [kg/s]
        self.C_dot = self.m_dot * self.cp  # Heat capacity rate of air [kW/k]
        self.enthalpy = HAPropsSI('H', 'T', self.T + 273.15, 'R', self.RH / 100, 'P', 101325) / 1000
        self.velocity = None
