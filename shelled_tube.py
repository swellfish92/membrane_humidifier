from material_properties import *
import numpy as np
import sympy as sp
from CoolProp.HumidAirProp import HAPropsSI
import matplotlib.pyplot as plt


class module_class:
    def __init__(self, var_dict):
        try:
            self.module_characteristic = var_dict[
                'tube_characteristic']  # Tube characteristic of module: aligned or staggered
            self.L = var_dict['L']  # Length of hollow fiber [m]
            self.D_module = var_dict['D_module']  # Diameter of membrane module [m]
            self.d_o = var_dict['d_o']  # Fiber outer diameter [m]
            self.d_i = var_dict['d_i']  # Fiber inner diameter [m]
            self.D_w = var_dict['D_w']  # Effective moisture diffusivity in membrane [m2/s]
            self.k_f = var_dict['k_f']  # Thermal conductivity of fiber [kW/mK]
            self.n_fiber = var_dict['n_fiber']  # Thermal conductivity of fiber [kW/mK]
        except:
            print('필요한 인자 값이 부족합니다.')
            raise IOError

        self.d_avg = (self.d_i + self.d_o) / 2  # Arithmetic mean diameter of fiber [m]
        self.thick_f = (self.d_o - self.d_i) / 2  # Thickness of fiber [m]
        self.A_tot_m = self.n_fiber * self.L * np.pi * self.d_o  # Total outer surface area of fibers [m2]

        # thickness of membrane layer
        self.thickness_membrane_porous = self.d_o - self.d_i  # [m]

        # Packing fraction = Hollow fiber cross section area / (H X D)
        self.packing_fraction = self.n_fiber * self.d_o ** 2 / (self.D_module ** 2)  # [m^2 / m^2]
        # print(self.packing_fraction)

        # Air side cross area
        self.air_side_cross_area = 3.14 / 4 * self.D_module ** 2 * (1 - self.packing_fraction)

        # Packing density = Hollow fiber cylinder surface area / module volume
        self.packing_density = 4 * self.n_fiber * self.d_o / (self.D_module ** 2)  # [m^2 / m^3]

        # section area of shell and tube module
        self.A_module = 3.14 * self.D_module ** 2 / 4
        # shell air side section area
        self.A_sec_shell = self.A_module * (1 - self.packing_fraction)

        # Hydraulic diameter
        self.d_h_shell = (1 - self.packing_fraction) * self.D_module ** 2 / (self.n_fiber * self.d_o + self.D_module)

        self.d_h_tube = self.d_i

    def simulate(self, water_in, air_in, constants):

        self.Nu_lim = 3.658  # Reference : Coupled heat and mass transfer in a counter flow hollow fiber membrane module for air humidification

        # calculate inlet air velocity
        air_in.velocity = air_in.V_dot / (3600 * self.A_sec_shell)  # [m/s]

        water_in.velocity = (water_in.m_dot / water_in.rho) / (
                    self.n_fiber * 3.14 / 4 * self.d_i ** 2)  # Inlet water velocity [m/s]

        self.Re_a = (air_in.velocity * self.d_h_shell) / constants.nu_a  # Reynolds number of inlet air [-]

        self.Re_w = (water_in.velocity * self.d_i) / constants.nu_w_in  # Reynolds number of water [-]

        # friction factor of laminar flow in water side
        self.friction_factor_w = 64 / self.Re_w

        # friction factor of laminar flow in air side
        self.friction_factor_a = 41.3 / self.Re_a  # 41.3 from Lz.Zhang

        # pressure drop of laminar flow
        self.pressure_drop_a = self.friction_factor_a * self.L * air_in.rho * air_in.velocity ** 2 / (
                    2 * self.d_h_shell)

        self.pressure_drop_w = self.friction_factor_w * self.L * water_in.rho * water_in.velocity ** 2 / (
                    2 * self.d_h_tube)

        self.Sc_a = constants.nu_a / constants.D_va  # Schmidt number of inlet air [-]

        self.Sh_a = (0.53 - 0.58 * self.packing_fraction) * self.Re_a ** 0.53 * self.Sc_a ** 0.33  # Sherwood number for air [-]

        self.alpha_a = constants.k_air / (air_in.rho * air_in.cp)  # thermal diffusivity of air

        # Lewis number
        self.Le_a = constants.Pr_a / self.Sc_a

        # Nusselt number of air side
        self.Nus_a = self.Sh_a * self.Le_a ** 0.33

        # Nusselt number of water side
        self.Nus_w = self.Nu_lim + (0.085 * self.Re_w * constants.Pr_w * self.d_h_tube / self.L) / (
                    1 + 0.047 * self.Re_w * constants.Pr_w * self.d_h_tube / self.L) ** 0.67

        # Heat and mass transfer coefficients at air-water interface in membrane surface
        self.h_a = (self.Nus_a * constants.k_air / self.d_h_shell) / 1000  # Convective heat transfer coefficient of air [kW/m2K]

        self.h_w = (self.Nus_w * constants.k_water / self.d_h_tube) / 1000  # convective heat transfer coefficient of water [kW/m2K]

        self.h_m_a = constants.D_va * self.Sh_a / self.d_h_shell  # Convective mass transfer coefficient of air [m/s]

        # heat and mass transfer coefficients of membrane
        self.effective_diffusivity_membrane = self.D_w

        # mass transfer resistance of porous membrane
        self.Resistance_mass_transfer_membrane = self.thickness_membrane_porous / self.effective_diffusivity_membrane

        self.effective_mass_diffusivity = self.thickness_membrane_porous / self.Resistance_mass_transfer_membrane

        # thermal conductance of porous membrane layer
        self.conductance_porous_membrane = self.k_f

        # Overall heat and mass transfer coefficients
        # mass transfer resistance of water phase is neglected because it is so small

        # Overall mass transfer coefficient
        self.h_m_tot = 1 / (self.thickness_membrane_porous / self.effective_mass_diffusivity * (
                    self.d_o / self.d_avg) + 1 / self.h_m_a)  # [m/s]

        # Overall heat transfer coefficient
        self.h_tot = 1 / (1 / self.h_w * (self.d_o / self.d_i) + self.thickness_membrane_porous / (
            self.conductance_porous_membrane) * (self.d_o / self.d_avg) + 1 / self.h_a)  # [kW/m2K]

        self.NTU = (self.h_tot * self.A_tot_m) / (air_in.m_dot * air_in.cp)  # Number of heat transfer units [-]

        # self.NTU = 1

        self.NTU_m = (self.h_m_tot * self.A_tot_m) / (air_in.m_dot / air_in.rho)  # Number of mass transfer units [-]

        print(self.NTU, self.NTU_m)

        # self.NTU_m = 2
        # print(self.A_tot_m, air_in.m_dot)
        # print(self.h_tot, self.h_m_tot)
        # print(self.NTU, self.NTU_m)
        # raise IOError

        self.Lewis_effective = self.NTU / self.NTU_m  # Lewis number

        self.R_thermal_membrane_porous = self.thickness_membrane_porous / self.conductance_porous_membrane

        self.R_thermal_water = 1 / self.h_w

        self.R_thermal_air = 1 / self.h_a

        self.R_mass_air = 1 / self.h_m_a

    def numerical(self, water_in, air_in, constants, grids=[25]):
        # EES 코드 포팅 과정에서 grid 셋업이 생긴 것으로 추정. append를 사용하도록 추후 수정할 것

        # inlet & constants setting
        self.simulate(water_in, air_in, constants)

        # 1-D counter flow grid setting
        n = grids[0]

        # Initialize arrays
        Ta = np.zeros(n + 1)
        wa = np.zeros(n + 1)
        Tw = np.zeros(n + 1)
        ww = np.zeros(n + 1)

        # Define initial conditions
        Ta[1] = air_in.T
        wa[1] = HAPropsSI('W', 'T', air_in.T + 273.15, 'R', air_in.RH / 100, 'P', 101325)

        # initial water side outlet assumption
        Tw[1] = air_in.T
        ww[1] = HAPropsSI('W', 'T', water_in.T + 273.15, 'R', 1, 'P', 101325)

        # Maximum number of iterations
        iteration = 0
        max_iterations = 1000
        tolerance = 0.05

        while iteration < max_iterations:
            # Perform calculations
            for i in range(1, n):
                Ta[i + 1] = (1 / n) * (self.NTU * (Tw[i] - Ta[i])) + Ta[i]

                wa[i + 1] = (1 / n) * (self.NTU_m * (ww[i] - wa[i])) + wa[i]

                Tw[i + 1] = Tw[i] + (1 / n) * ((self.NTU * (air_in.C_dot / water_in.c) * (Tw[i] - Ta[i])) + (
                            self.NTU_m * (air_in.m_dot / water_in.c) * 2501 * (ww[i] - wa[i])))

                ww[i + 1] = HAPropsSI('W', 'T', Tw[i] + 273.15, 'R', 1, 'P', 101325)
            dT_w = water_in.T - Tw[n]

            if abs(dT_w) <= tolerance:
                break
            else:
                # initial water side outlet assumption
                Tw[1] += dT_w
                ww[1] = HAPropsSI('W', 'T', Tw[1] + 273.15, 'R', 1, 'P', 101325)

        # Plot wa[i] versus i (or n)
        # plt.figure(figsize=(8, 6))
        # plt.plot(range(n + 1), wa, label='wa[i]')
        # plt.plot(range(n + 1), ww, label='ww[i]')
        # plt.xlabel('i (discretized points)')
        # plt.ylabel('wa[i] (Specific Humidity)')
        # plt.title('Specific Humidity Distribution across Discretized Points')
        # plt.grid(True)
        # plt.legend()
        # plt.show()
        #
        # plt.figure(figsize=(8, 6))
        # plt.plot(range(n + 1), Ta, label='Ta[i]')
        # plt.plot(range(n + 1), Tw, label='Tw[i]')
        # plt.xlabel('i (discretized points)')
        # plt.ylabel('Ta[i] (Temperature)')
        # plt.title('Temperature Distribution across Discretized Points')
        # plt.grid(True)
        # plt.legend()
        # plt.show()

        # Calculate output values
        T_a_out = Ta[n]
        w_a_out = wa[n]
        T_w_out = Tw[1]
        w_w_out = ww[1]

        try:
            air_out = air(T_a_out, HAPropsSI('R', 'T', T_a_out + 273.15, 'W', w_a_out, 'P', 101325) * 100,
                          air_in.V_dot, air_in.rho, air_in.cp)
        except Exception as e:
            air_out = air(T_a_out, HAPropsSI('R', 'T', T_a_out + 273.15, 'R', 1, 'P', 101325) * 100,
                          air_in.V_dot, air_in.rho, air_in.cp)

        error_humrat = (air_out.w - w_a_out) / w_a_out * 100

        water_out = water(T_w_out, water_in.V_dot, water_in.rho, water_in.cp)

        # error of heat and moisture transfer
        dQ_a = air_in.m_dot * (air_out.enthalpy - air_in.enthalpy)
        dQ_w = water_in.m_dot * (water_in.enthalpy - water_out.enthalpy)
        error = (dQ_a - dQ_w) / dQ_a * 100

        # humidification rate [kg/s]
        m_humidified = air_in.m_dot * (air_out.w - air_in.w)

        # Effectiveness
        self.Effectiveness_mass = (air_out.w - air_in.w) / (ww[1] - air_in.w)

        numerical_calc_dict = {'error_humrat': error_humrat, 'dQ_a': dQ_a, 'dQ_w': dQ_w, 'error': error,
                               'm_humidified': m_humidified}

        return water_out, air_out, numerical_calc_dict

