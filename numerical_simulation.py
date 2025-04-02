from material_properties import *
import numpy as np
import sympy as sp
import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI
# from Material_properties import air
import matplotlib.pyplot as plt

# 아 집에가고싶다

class module_class:
    def __init__(self, var_dict):
        try:
            self.align = var_dict['tube_characteristic']
            self.L = var_dict['L']  # Length of hollow fiber [m]
            self.H = var_dict['H']  # Height of membrane module[m]
            self.D = var_dict['D']  # Depth of membrane module [m]
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
        self.A_tot_m_in = self.n_fiber * self.L * np.pi * self.d_i
        self.n_ch = self.A_tot_m / (self.L * self.D)  # Equilibrium channel number [-]

        # Packing fraction = Hollow fiber cross section area / (H X D)

        self.packing_fraction = self.n_fiber * np.pi * self.d_o ** 2 / (4 * self.H * self.D)  # [m^2 / m^2]
        print('packing fraction: {}'.format(self.packing_fraction))

        # Packing density = Hollow fiber cylinder surface area / module volume
        self.packing_density = self.n_fiber * np.pi * self.d_o / (self.H * self.D)  # [m^2 / m^3]

        # Tube pitch calculation
        p_t = sp.symbols('p_t')
        eqn_p_t = p_t ** 2 * self.n_fiber - (p_t + self.D - self.d_o) * (p_t + self.H - self.d_o)
        solve_p_t = sp.solve(eqn_p_t, p_t)

        self.pitch_tranverse = float(max(solve_p_t))  # Calculated p_t [m]
        self.pitch_longitudinal = self.pitch_tranverse  # Tranverse pitch [m]
        self.pitch_diagonal = np.sqrt(self.pitch_longitudinal ** 2 + (self.pitch_tranverse / 2) ** 2)

        # self.pitch_tranverse = 1.66 / 1000  # Calculated p_t [m]
        # self.pitch_longitudinal = 1.66 / 1000  # Tranverse pitch [m]
        # self.pitch_diagonal = np.sqrt(self.pitch_longitudinal**2 + (self.pitch_tranverse / 2)**2)

    def simulate(self, water_in, air_in, constants):
        # calculate inlet air velocity
        air_in.velocity = air_in.V_dot / (3600 * self.H * self.L)  # [m/s]
        # Max velocity in aligned tube
        self.align = False
        if self.align == True:
            max_air_velocity = air_in.velocity * self.pitch_tranverse / (self.pitch_tranverse - self.d_o)
        elif self.align == False:
            max_air_velocity = air_in.velocity * self.pitch_tranverse / (2 * (self.pitch_diagonal - self.d_o))

        self.Re_a = (max_air_velocity * self.d_o) / constants.nu_a  # Reynolds number of inlet air [-]

        print(max_air_velocity / air_in.velocity)

        print(self.Re_a / max_air_velocity)

        self.Sc_a = constants.nu_a / constants.D_va  # Schmidt number of inlet air [-]

        print('Re: {}'.format(self.Re_a))

        water_in.velocity = (water_in.m_dot / water_in.rho) / (self.n_fiber * np.pi * 0.25 * self.d_o ** 2)  # Inlet water velocity [m/s]

        self.Re_w = (water_in.velocity * self.d_o) / constants.nu_w_in  # Reynolds number of inlet water [-]

        # Heat and mass transfer coefficients
        self.Nus_a = 1.04 * self.Re_a ** 0.4 * constants.Pr_a ** 0.33  # Nusselt number for air [-]
        self.h_a = (self.Nus_a * constants.k_air / self.d_o) / 1000  # Convective heat transfer coefficient of air [kW/m2K]
        self.Sh_a = (0.53 - 0.58 * self.packing_fraction) * self.Re_a ** 0.53 * self.Sc_a ** 0.33  # Sherwood number for air [-]
        self.h_m_a = constants.D_va * self.Sh_a / self.d_o  # Convective mass transfer coefficient of air [m/s]

        # Total heat and mass transfer coefficients

        R_h_membrane = (self.thick_f / self.k_f)  # K/kW
        R_h_air = 1 / self.h_a  # K/kW
        R_h_tot = (self.thick_f / self.k_f) * (self.d_o / self.d_avg) + 1 / self.h_a  # Total heat transfer resistance [m2K/kW]

        R_m_membrane = (self.thick_f / self.D_w) * (self.d_o / self.d_avg)  # s/m
        R_m_air = 1 / self.h_m_a  # s/m

        R_m_tot = (self.thick_f / self.D_w) * (self.d_o / self.d_avg) + 1 / self.h_m_a  # Total mass transfer resistance [s/m]
        self.h_tot = 1 / R_h_tot  # Total heat transfer coefficient [kW/m2K]
        self.h_m_tot = 1 / R_m_tot  # Total mass transfer coefficient [m/s]

        self.NTU = (self.h_tot * self.A_tot_m) / (air_in.m_dot * air_in.cp)  # Number of heat transfer units [-]

        self.NTU_m = (self.h_m_tot * self.A_tot_m) / (air_in.m_dot / air_in.rho) # Number of mass transfer units [-]
        # self.NTU_m = 0.48
        print('R_h_a: {}'.format(R_h_air), 'R_h_membrane: {}'.format(R_h_membrane))
        print('R_m_a: {}'.format(R_m_air), 'R_m_membrane: {}'.format(R_m_membrane))
        print('NTU: {}'.format(self.NTU), 'NTU_m: {}'.format(self.NTU_m))


    def numerical(self, water_in, air_in, constants, grids=[25, 25]):

        # inlet & constants setting
        self.simulate(water_in, air_in, constants)

        # grid setting
        nx = grids[0]
        ny = grids[1]

        # Initialize arrays
        Ta = np.zeros((nx + 2, ny + 2))
        wa = np.zeros((nx + 2, ny + 2))
        Tw = np.zeros((nx + 2, ny + 2))
        ww = np.zeros((nx + 2, ny + 2))

        # Define initial conditions
        Ta[1, :] = air_in.T
        wa[1, :] = HAPropsSI('W', 'T', air_in.T + 273.15, 'R', air_in.RH/100, 'P', 101325)

        Tw[:, 1] = water_in.T
        ww[:, 1] = HAPropsSI('W', 'T', water_in.T + 273.15, 'R', 1, 'P', 101325)

        # Perform calculations
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                Ta[i + 1, j] = (1 / ny) * (self.NTU * (Tw[i, j] - Ta[i, j])) + Ta[i, j]

                wa[i + 1, j] = (1 / ny) * (self.NTU_m * (ww[i, j] - wa[i, j])) + wa[i, j]

                # Tw[i, j + 1] = -(1 / nx) * ((self.NTU * (air_in.C_dot / water_in.c) * (Tw[i, j] - Ta[i, j])) +
                #                             (self.NTU_m * (air_in.m_dot / water_in.c) * 2501 * (ww[i, j] - wa[i, j]))) + Tw[i, j]

                Tw[i, j + 1] = -(1 / nx) * ((self.NTU * (air_in.C_dot / water_in.c) * (Tw[i, j] - Ta[i, j])) + (
                        self.NTU_m * (air_in.m_dot / water_in.c) * 2501 * (ww[i, j] - wa[i, j]))) + Tw[i, j]


                ww[i, j + 1] = HAPropsSI('W', 'T', Tw[i, j + 1] + 273.15, 'R', 1, 'P', 101325)

        # Ta[i + 1, j] = (1 / ny) * ((NTU * (Tw[i, j] - Ta[i, j]))) + Ta[i, j]
        #
        # wa[i + 1, j] = (1 / ny) * ((NTU_m * (ww[i, j] - wa[i, j]))) + wa[i, j]
        #
        # Tw[i, j + 1] = -(1 / nx) * ((NTU * (C_a / C_w) * (Tw[i, j] - Ta[i, j])) + (
        #             NTU_m * (m_a_in / C_w) * 2501 * (ww[i, j] - wa[i, j]))) + Tw[i, j]
        #
        # ww[i, j + 1] = w_eq_w0(Tw[i, j + 1])

        # Calculate output values
        T_a_out = np.mean(Ta[nx + 1, 1:ny + 1])
        w_a_out = np.mean(wa[nx + 1, 1:ny + 1])
        T_w_out = np.mean(Tw[1:nx + 1, ny + 1])
        w_w_out = np.mean(ww[1:nx + 1, ny + 1])

        T_a_internal = Ta[1:nx + 1, 1:ny + 1]
        w_a_internal = wa[1:nx + 1, 1:ny + 1]


        # # 물리적 길이 설정
        # x_length = 240  # x 방향 길이 (mm)
        # y_length = 54
        # x = np.linspace(0, x_length, nx)  # 0에서 240mm까지 nx개의 점 생성
        # y = np.linspace(0, y_length, ny)  # y 방향은 그대로 그리드 크기로 유지
        # X, Y = np.meshgrid(x, y)
        #
        # # 컨투어 플롯 그리기
        # plt.figure(figsize=(8, 6))
        # contour = plt.contour(X, Y, T_a_internal, levels=10, cmap='viridis')  # levels 조정 가능
        # plt.clabel(contour, inline=True, fontsize=12)  # 레벨 값 표시
        # # plt.title('Air temperature (C)')
        # plt.xlabel('Length (mm)', fontsize=12)
        # plt.ylabel('Depth (mm)', fontsize=12)
        # plt.show()
        #
        # # 컨투어 플롯 그리기
        # plt.figure(figsize=(8, 6))
        # contour = plt.contour(X, Y, w_a_internal, levels=10, cmap='viridis')  # levels 조정 가능
        # plt.clabel(contour, inline=True, fontsize=12)  # 레벨 값 표시
        # # plt.title('Air temperature (C)')
        # plt.xlabel('Length (mm)', fontsize=12)
        # plt.ylabel('Depth (mm)', fontsize=12)
        # plt.show()
        #
        #
        # Tw_internal = Tw[1:nx + 1, 1:ny + 1]
        #
        # x = np.linspace(0, nx, nx)
        # y = np.linspace(0, ny, ny)
        # X, Y = np.meshgrid(x, y)
        #
        # print(np.mean(Tw[1:nx + 1, ny + 1]))
        # print(np.mean(Tw[nx + 1, 1]))
        #
        # # 컨투어 플롯 그리기
        # plt.figure(figsize=(8, 6))
        # contour = plt.contourf(X, Y, Tw_internal, levels=50, cmap='viridis')
        # plt.colorbar(contour)
        # plt.title('Temperature (C)')
        # plt.xlabel('Grid X')
        # plt.ylabel('Grid Y')
        # plt.show()

        try:
            air_out = air(T_a_out, HAPropsSI('R', 'T', T_a_out + 273.15, 'W', w_a_out, 'P', 101325) * 100,
                          air_in.V_dot, air_in.rho, air_in.cp)
        except Exception as e:
            air_out = air(T_a_out, HAPropsSI('R', 'T', T_a_out + 273.15, 'R', 1, 'P', 101325) * 100,
                          air_in.V_dot, air_in.rho, air_in.cp)
        error_humrat = (air_out.w - w_a_out) / w_a_out * 100

        water_out = water(T_w_out, water_in.V_dot, water_in.rho, water_in.cp)

        # error of heat and moisture transfer
        dQ_a = air_in.m_dot * (air_out.enthalpy - air_in.enthalpy)  # kW
        dQ_w = water_in.m_dot * (water_in.enthalpy - water_out.enthalpy)  # kW
        error = (dQ_a - dQ_w) / dQ_a * 100

        # humidification rate
        m_humidified = air_in.m_dot * (air_out.w - air_in.w)  # kg/s

        # Mass transfer effectiveness
        self.Effectiveness_mass = (w_a_out - wa[1,1]) / (ww[1,1] - wa[1,1])

        numerical_calc_dict = {'error_humrat': error_humrat, 'dQ_a': dQ_a, 'dQ_w': dQ_w, 'error': error,
                               'm_humidified': m_humidified}
        print(T_a_out, w_a_out, T_w_out)
        return water_out, air_out, numerical_calc_dict
#
#
