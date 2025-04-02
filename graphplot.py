import matplotlib.pyplot as plt
import pandas as pd
from CoolProp.CoolProp import HAPropsSI

import matplotlib
matplotlib.use('TkAgg')

data = pd.read_csv('./export_result.csv', encoding='ANSI')
print(data)

data['m_humidified'] = data['m_humidified'] * 3600 # kg/s -> g/h
# 외기/보충 공기/배기 조건 설정
m_oa = 0.01654 # 외기 유량
omega_oa = 0.0021 # 외기 절대 습도(0℃, RH 55%)
m_makeup = 0.4168 # 보충 공기 유량 (전체 유량 - 외기 유량)
omega_mu = 0.0086 # 보충 공기 절대 습도 (실내 조건: 20℃, RH 58.92%)
m_ea = 0.4333 # 배기 유량 = 전체 유입량 (steady-state)
omega_ea = 0.008 # 배기 습도 = 실내 조건 기준

left = m_oa * omega_oa + m_makeup * omega_mu
# 외기 도입량 * (실내 절대 습도 - 실외 절대 습도) = 가습량 = 실내기 유량 * (실내기 출구 절대 습도 - 실내기 입구 절대 습도)
right = m_ea * omega_ea
delta = right - left

print("[Steady-State 검토]")
print("유입 수분량 (left):", left)
print("유출 수분량 (right):", right)
print("차이 (right - left):", delta)

if abs(delta) < 1e-5:
    print("Steady-state")
else:
    print("Steady-state 아님")

# 수분 손실량 계산 (g/h)
data['m_loss'] = (m_ea * omega_ea - (m_oa * omega_oa + m_makeup * omega_mu)) * 3600 * 1000

# 조건 판단
data['design_cond'] = data.apply(
    lambda row: 'PASS' if row['m_humidified'] >= row['m_loss'] else 'FAIL',
    axis=1
)

fig, ax = plt.subplots(1, 3)
for v_water in data['v_w_in'].unique():
    sliced_data = data[(data['v_w_in'] == v_water)]
    valid_data = sliced_data[sliced_data['design_cond'] == 'PASS']

    # v_lpm = v_water * 1000 * 60  # m³/s → L/min 변환
    v_lpm = v_water * 1800
    label_T_a_out = f"T_a_out: v = x {v_lpm:.0f}" # 소수점 1자리 정수로 출력
    label_rh_a_out = f"RH_out: v = x {v_lpm:.0f}"
    label_m_humidified = f"Humidified : v = x {v_lpm:.0f}"

    label_v=f"{v_water:.4f}"

    print(sliced_data)

    ax[0].plot(sliced_data['thickness'], sliced_data['T_a_out'], label=label_T_a_out)
    ax[0].scatter(valid_data['thickness'], valid_data['T_a_out'], marker="*")
    # ax[0].legend()
    ax[1].plot(sliced_data['thickness'], sliced_data['rh_a_out'], label=label_rh_a_out)
    ax[1].scatter(valid_data['thickness'], valid_data['rh_a_out'], marker="*")
    # ax[1].legend()
    ax[2].plot(sliced_data['thickness'], sliced_data['m_humidified'], label=label_m_humidified)
    ax[2].scatter(valid_data['thickness'], valid_data['m_humidified'], marker="*")
    # ax[2].legend()

    ax[0].set_xlabel('Thickness [m]')
    ax[0].set_ylabel('Outlet Air Temperature [℃]')
    ax[0].set_title('Outlet Air Temperature vs Membrane Thickness')
    ax[0].legend()

    ax[1].set_xlabel('Thickness [m]')
    ax[1].set_ylabel('Outlet Relative Humidity [%]')
    ax[1].set_title('Outlet Relative Humidity vs Membrane Thickness')
    ax[1].legend()

    ax[2].set_xlabel('Thickness [m]')
    ax[2].set_ylabel('Humidification Amount [g/h]')
    ax[2].set_title('Humidification amount vs Membrane Thickness')
    ax[2].legend()

fig.text(0.9, 0.05, '* V is the ratio of the multiplier of the flowrate of each membrane tube (1/1800 lpm)', fontsize = 12, horizontalalignment='right')


ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()
    # print(data)

