from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_plotting
from material_properties import *
from CoolProp.HumidAirProp import HAPropsSI
import shelled_tube as Shellandtube
import numerical_simulation as rectangular_module
import traceback
import os
from pathlib import Path


def get_onedrive_path():
    # On Windows, the OneDrive path is stored in an environment variable
    onedrive_path = os.getenv('OneDrive')

    if not onedrive_path:
        raise EnvironmentError("OneDrive path not found. Ensure OneDrive is installed correctly.")

    return onedrive_path


# Function to calculate missing humidity parameters
def calculate_humidity_params(T, RH=None, abs_hum=None, T_wb=None):
    if RH is not None:
        abs_hum = HAPropsSI('W', 'T', T + 273.15, 'P', 101325, 'RH', RH / 100.0)
        T_wb = HAPropsSI('Twb', 'T', T + 273.15, 'P', 101325, 'RH', RH / 100.0) - 273.15
    elif abs_hum is not None:
        RH = HAPropsSI('RH', 'T', T + 273.15, 'P', 101325, 'W', abs_hum) * 100
        T_wb = HAPropsSI('Twb', 'T', T + 273.15, 'P', 101325, 'W', abs_hum) - 273.15
    elif T_wb is not None:
        RH = HAPropsSI('RH', 'T', T + 273.15, 'P', 101325, 'Twb', T_wb + 273.15) * 100
        abs_hum = HAPropsSI('W', 'T', T + 273.15, 'P', 101325, 'Twb', T_wb + 273.15)
    else:
        raise ValueError("One of RH, abs_hum, or T_wb must be provided")

    return RH, abs_hum, T_wb


def membrane_simulation(input_dict, module_var_dict):
    final_result_dict = {
        # input
        'n_fiber': [],
        'T_a_in': [],
        'w_a_in': [],
        'rh_a_in': [],
        'h_a_in': [],
        'v_a_in': [],
        'T_w_in': [],
        'h_w_in': [],
        'v_w_in': [],
        # output
        'T_a_out': [],
        'w_a_out': [],
        'rh_a_out': [],
        'h_a_out': [],
        'T_w_out': [],
        'h_w_out': [],
        'error_humrat': [],
        'dQ_a': [],
        'dQ_w': [],
        'error': [],
        'm_humidified': [],
        'Reynolds_a': [],
        'Nusselt_a': [],
        'sherwood_a': [],
        'h_total': [],
        'h_m_total': [],
        'NTU': [],
        'NTU_m': [],
        'Effectiveness_mass': []
    }

    n_fiber_iter = module_var_dict['n_fiber'] # Module의 Hollow fiber 총 개수

    temp_air_in_iter = input_dict['t_air_in']
    temp_water_in_iter = input_dict['t_water_in']
    v_air_in_iter = input_dict['v_air_in']
    v_water_in_iter = input_dict['v_water_in']

    # air_in의 속성값 3개에 대한 로딩을 시도한다. 3개 전부 None일 경우는 진행이 불가능하므로 중단한다.
    try:
        RH_air_in = input_dict['rh_air_in']
    except:
        RH_air_in = None
    try:
        T_wb_air_in = input_dict['wbt_air_in']
    except:
        T_wb_air_in = None
    try:
        abs_hum_air_in = input_dict['abs_hum_air_in']
    except:
        abs_hum_air_in = None

    if RH_air_in == None and T_wb_air_in == None and abs_hum_air_in == None:
        raise IOError('air_in에 대한 습도설정값이 정의되지 않았습니다. rh_air_in, wbt_air_in, abs_hum_air_in 중 최소 하나를 정의하십시오.')

        # RH, HR, WBT 중 아무거나 하나만 넣어도 RH 값 계산
    try:
        RH_air_in, abs_hum_air_in, T_wb_air_in = calculate_humidity_params(temp_air_in_iter,
                                                                           RH=RH_air_in,
                                                                           abs_hum=abs_hum_air_in,
                                                                           T_wb=T_wb_air_in)

        final_result_dict['n_fiber'].append(n_fiber_iter)
        final_result_dict['T_a_in'].append(temp_air_in_iter)
        final_result_dict['rh_a_in'].append(RH_air_in)
        final_result_dict['v_a_in'].append(v_air_in_iter)
        final_result_dict['T_w_in'].append(temp_water_in_iter)
        final_result_dict['v_w_in'].append(v_water_in_iter)


        if 'H' in module_var_dict and 'D' in module_var_dict:
            module = rectangular_module.module_class(module_var_dict)
            const_class = const_var()

        # Checking if 'D_module' exists in the dictionary
        elif 'D_module' in module_var_dict:
            module = Shellandtube.module_class(module_var_dict)
            const_class = const_var()
        else:
            raise ValueError(
                "Neither 'H' and 'D' nor 'D_module' is present in the module_var_dict.")

        water_inlet = water(temp_water_in_iter, v_water_in_iter, 1000, 4.179)
        air_inlet = air(temp_air_in_iter, RH_air_in, v_air_in_iter, 1.2, 1.02)

        # 1-D numerical
        # water_outlet, air_outlet, result = module.numerical(water_inlet, air_inlet, const_class,
        #                                                     [60])

        # 2-D numerical
        water_outlet, air_outlet, result = module.numerical(water_inlet, air_inlet, const_class,
                                                            [40, 40])

        # Output parameters
        final_result_dict['w_a_in'].append(air_inlet.w)
        final_result_dict['T_a_out'].append(air_outlet.T)
        final_result_dict['w_a_out'].append(air_outlet.w)
        final_result_dict['rh_a_out'].append(air_outlet.RH)
        final_result_dict['h_a_in'].append(air_inlet.enthalpy)
        final_result_dict['h_a_out'].append(air_outlet.enthalpy)
        final_result_dict['T_w_out'].append(water_outlet.T)
        final_result_dict['h_w_in'].append(water_inlet.enthalpy)
        final_result_dict['h_w_out'].append(water_outlet.enthalpy)
        final_result_dict['error_humrat'].append(result['error_humrat'])
        final_result_dict['dQ_a'].append(result['dQ_a'])
        final_result_dict['dQ_w'].append(result['dQ_w'])
        final_result_dict['error'].append(result['error'])
        final_result_dict['m_humidified'].append(result['m_humidified'])
        final_result_dict['Reynolds_a'].append(module.Re_a)
        final_result_dict['Nusselt_a'].append(module.Nus_a)
        final_result_dict['sherwood_a'].append(module.Sh_a)
        final_result_dict['h_total'].append(module.h_tot)
        final_result_dict['h_m_total'].append(module.h_m_tot)
        final_result_dict['NTU'].append(module.NTU)
        final_result_dict['NTU_m'].append(module.NTU_m)
        final_result_dict['Effectiveness_mass'].append(module.Effectiveness_mass)
    except:
        final_result_dict['n_fiber'].append(n_fiber_iter)
        final_result_dict['T_a_in'].append(temp_air_in_iter)
        final_result_dict['rh_a_in'].append(RH_air_in)
        final_result_dict['v_a_in'].append(v_air_in_iter)
        final_result_dict['T_w_in'].append(temp_water_in_iter)
        final_result_dict['v_w_in'].append(v_water_in_iter)
        final_result_dict['w_a_in'].append('error')
        final_result_dict['T_a_out'].append('error')
        final_result_dict['w_a_out'].append('error')
        final_result_dict['rh_a_out'].append('error')
        final_result_dict['h_a_in'].append('error')
        final_result_dict['h_a_out'].append('error')
        final_result_dict['T_w_out'].append('error')
        final_result_dict['h_w_in'].append('error')
        final_result_dict['h_w_out'].append('error')
        final_result_dict['error_humrat'].append('error')
        final_result_dict['dQ_a'].append('error')
        final_result_dict['dQ_w'].append('error')
        final_result_dict['error'].append('error')
        final_result_dict['m_humidified'].append('error')
        final_result_dict['Reynolds_a'].append('error')
        final_result_dict['Nusselt_a'].append('error')
        final_result_dict['sherwood_a'].append('error')
        final_result_dict['h_total'].append('error')
        final_result_dict['h_m_total'].append('error')
        final_result_dict['NTU'].append('error')
        final_result_dict['NTU_m'].append('error')
        final_result_dict['Effectiveness_mass'].append('error')
        # print(f"Error occurred: {e}")
        print(traceback.format_exc())

    return final_result_dict

    result_df = pd.DataFrame(final_result_dict)

    return result_df


# calculate simulation result
if __name__ == '__main__':


    df_dict_arr = []

    for depth in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
        for v_water_in in [1, 2, 3]:
            # 실 체적
            room_volume = 8 * 5 * 2.4

            # 외기조건
            t_ext_air = 0
            rh_ext_air = 55

            # 환기조건
            ventilation_ratio = 0.5

            # 멤브레인 특성값
            D_w = 1.85E-06
            k_f = 4.63E-03 / 1000  # kW/mK
            length = 0.4   # length of fiber [m]
            height = 0.5    # height of fiber [m]
            # depth = 0.04   # depth of fiber [m] (Thickness)

            # 중공사 밀도가 일정하게 유지되도록 중공사 수를 조정. 중공사당 유량도 동일하게 수정
            #n_fiber = 1800
            # 기존 값 = 0.054 depth, 0.155 height에 1800개이므로 이에 맞추어 계산
            density_fiber = 1800 / (0.054 * 0.155)
            n_fiber = height * depth * density_fiber
            v_water_in = v_water_in * (n_fiber / 1800)

            # 모듈에 들어가는 물 및 공기의 조건 정의 (공기 - 실내공기, 물 - 급수온도 (해당 절기 시수온도))
            t_air_in = 20
            rh_air_in = 58.92
            t_water_in = 10
            v_air_in = 1300
            # v_water_in = 1

            # Module information
            module_var_dict = {
                'tube_characteristic': False, # Tube characteristic: set True to aligned and false to staggered and shell&tube

                # Cross flow module specification
                'L': length,  # Length of hollow fiber [m]
                'H': height,  # Height of membrane module[m]
                'D': depth,  # Depth of membrane module [m]

                # shell and tube module specification
                # 'L': 0.29,  # Length of hollow fiber [m]
                # 'D_module': 33.27 / 1000,  # module diameter of shell and tube module [m]

                'd_o': 1.1 / 1000,  # Fiber outer diameter [m]
                'd_i': 0.85 / 1000,  # Fiber inner diameter [m]

                # Assumed values
                'D_w': D_w,  # Effective moisture diffusivity in membrane [m2/s]
                'k_f': k_f,  # Thermal conductivity of fiber [kW/mK]
                'n_fiber': n_fiber  # Number of hollow fiber
            }

            input_dict = {
                't_air_in' : t_air_in,
                'rh_air_in' : rh_air_in,
                't_water_in' : t_water_in,
                'v_air_in' : v_air_in,
                'v_water_in' : v_water_in
            }

            result_simulation = membrane_simulation(input_dict, module_var_dict)
            result_df = pd.DataFrame(result_simulation)
            result_df['thickness'] = depth
            result_df['v_w_in'] = result_df['v_w_in'] / n_fiber

            # 외기조건과 내기조건, 가열후 조건을 비교하여 가습량 및 수분손실량 정의
            # 시간당 수분손실량 정의
            outdoor_air_absh = CP.HAPropsSI('W', 'T', t_ext_air + 273.15, 'R', rh_ext_air/100, 'P', 101325)
            indoor_air_absh = CP.HAPropsSI('W', 'T', t_air_in + 273.15, 'R', rh_air_in/100, 'P', 101325)
            absh_loss = indoor_air_absh - outdoor_air_absh
            # 가습량 정의
            humidified_air_absh = CP.HAPropsSI('W', 'T', float(result_simulation['T_a_out'][0]) + 273.15, 'R', float(result_simulation['rh_a_out'][0]) / 100, 'P', 101325)
            absh_gain = humidified_air_absh - indoor_air_absh
            # (시간당 수분손실량 * 실체적) 과 (가습량 * 히트펌프 실내기 처리유량) 값을 비교
            if absh_loss * room_volume * ventilation_ratio <= absh_gain * v_air_in:
                result_df['design_cond'] = 'PASS'
            else:
                result_df['design_cond'] = 'FAIL'

            # 나머지 값들도 확인을 위해 기록
            result_df['absh_loss'] = absh_loss
            result_df['absh_gain'] = absh_gain

            df_dict_arr.append(result_df)
            print(result_simulation)



    result = pd.concat(df_dict_arr)
    result.to_csv('./export_result.csv', encoding='ANSI')
    # result_simulation.to_csv('./export.csv', encoding='ANSI')

# Finding D_w and k_f from experimental result
# if __name__ == '__main__':
#
#     D_w = []
#     k_f = []
#
#     D_w += [1 * 10 ** (-7), 1.1 * 10 ** (-7)]  # Effective moisture diffusivity in membrane [m2/s]
#     k_f += [1 * 10 ** (-7), 1.1 * 10 ** (-7)]  # Thermal conductivity of fiber [kW/mK]
#
#     T_a_out_exp = 24.1206896551724
#     w_a_out_exp = 0.0130606169585817
#     T_w_out_exp = 25.3865517241379
#
#     tolerance = 0.03  # 허용 오차
#
#     max_iterations = 100  # 최대 실행 횟수
#
#     # D_w_factor = D_w / 5  # D_w 변화량 조절 계수
#     # k_f_factor = k_f / 5
#
#     T_a_out = []
#     w_a_out = []
#     T_w_out = []
#
#     T_a_out_error = []
#     w_a_out_error = []
#     T_w_out_error = []
#
#     w_a_out_error_change = []
#     T_a_out_error_change = []
#     T_w_out_error_change = []
#
#     result_simulation = membrane_simulation(D_w[0], k_f[0])
#
#     T_a_out.append(result_simulation['T_a_out'][0])
#     w_a_out.append(result_simulation['w_a_out'][0])
#     T_w_out.append(result_simulation['T_w_out'][0])
#
#     T_a_out_error.append((T_a_out[0] - T_a_out_exp) / T_a_out_exp)
#     w_a_out_error.append((w_a_out[0] - w_a_out_exp) / w_a_out_exp)
#     T_w_out_error.append((T_w_out[0] - T_w_out_exp) / T_w_out_exp)
#
#     for i in range(1, max_iterations):
#         print(i)
#
#         result_simulation = membrane_simulation(D_w[i], k_f[i])
#
#         T_a_out.append(result_simulation['T_a_out'][0])
#         w_a_out.append(result_simulation['w_a_out'][0])
#         T_w_out.append(result_simulation['T_w_out'][0])
#
#         T_a_out_error.append((T_a_out[i] - T_a_out_exp) / T_a_out_exp)
#         w_a_out_error.append((w_a_out[i] - w_a_out_exp) / w_a_out_exp)
#         T_w_out_error.append((T_w_out[i] - T_w_out_exp) / T_w_out_exp)
#         print(T_a_out_error[i], w_a_out_error[i], T_w_out_error[i])
#
#         total_error = max(abs(T_a_out_error[i]), abs(w_a_out_error[i]), abs(T_w_out_error[i]))  # 최대 오차 사용
#         print(f"total_error {total_error}")
#
#         if total_error <= tolerance:
#             print("Converged within tolerance!")
#             print(D_w[i], k_f[i])
#             break
#         print(T_a_out, w_a_out, T_w_out)
#
#         D_w += [D_w[i] - ((D_w[i] - D_w[i - 1]) / (w_a_out[i] - w_a_out[i - 1]) * (w_a_out[i] - w_a_out_exp))]
#         k_f += [k_f[i] - ((k_f[i] - k_f[i - 1]) / (T_w_out[i] - T_w_out[i - 1]) * (T_w_out[i] - T_w_out_exp))]
#
#         print(D_w[i + 1], k_f[i + 1])
#
#     if i == max_iterations - 1:
#         print("Did not converge within max iterations.")
#
#     print(T_a_out, w_a_out, T_w_out)
