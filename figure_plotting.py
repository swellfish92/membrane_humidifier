import matplotlib.pyplot as plt
import pandas as pd
# import reading_csv as rc
from matplotlib import colormaps

color_arr = ['orange', 'red', 'blue', 'purple', 'green', 'navajowhite', 'lightcoral', 'plum', 'slategray', 'navy']

    # Plot function
def plotting(data, x_axis_label, y_axis_label_dict):
    # 딕셔너리는 무조건 최대 2개까지의 축을 가지며, 구조는 아래와 같다.
    # sample_dict = {
    #     'left_axis_name':[value_label_1, value_label_2, ...],
    #     'right_axis_name': [value_label_1, value_label_2, ...],
    #     'axis_labels':[left_axis_label, right_axis_label]
    # }
    # Plot the 'V_OC_TEG' column
    fig, ax1 = plt.subplots(figsize=(10, 6))
    # print(list(y_axis_label_dict.keys()))

    if len(list(y_axis_label_dict.keys())) == 2:
        minmax_value_arr = [[], []]
        for left_axis_label in y_axis_label_dict[list(y_axis_label_dict.keys())[0]]:
            print(left_axis_label)
            ax1.plot(data.index, data[left_axis_label], label=left_axis_label, color=color_arr[len(minmax_value_arr[0])])  # Using black line for V_OC_TEG data
            ax1.set_ylabel(y_axis_label_dict['axis_labels'][0], color='black')
            ax1.tick_params(axis='y', labelcolor='black')
            # select adequate y scale range
            minmax_value_arr[0].append(data[left_axis_label].min())
            minmax_value_arr[1].append(data[left_axis_label].max())

            # # select adequate major y scale ticker
            # ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
            #
            # # Enable minor ticks for ax1
            # ax1.minorticks_on()
            # ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
            # ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

            # Grid only for major ticks for ax1
            ax1.grid(True, which='major', linestyle='--', linewidth=0.5)
        delta = max(minmax_value_arr[1]) - min(minmax_value_arr[0])
        ax1.set_ylim(min(minmax_value_arr[0]) - delta*0.1, max(minmax_value_arr[1]) + delta*0.1)
        ax1.set_xlabel(x_axis_label)

    elif len(list(y_axis_label_dict.keys())) == 3:
        minmax_value_arr = [[], []]
        for left_axis_label in y_axis_label_dict[list(y_axis_label_dict.keys())[0]]:
            print(left_axis_label)
            ax1.plot(data.index, data[left_axis_label], label=left_axis_label, color=color_arr[len(minmax_value_arr[0])])  # Using black line for V_OC_TEG data
            ax1.set_ylabel(y_axis_label_dict['axis_labels'][0], color='black')
            ax1.tick_params(axis='y', labelcolor='black')
            # select adequate y scale range
            minmax_value_arr[0].append(data[left_axis_label].min())
            minmax_value_arr[1].append(data[left_axis_label].max())
        delta = max(minmax_value_arr[1]) - min(minmax_value_arr[0])
        ax1.set_ylim(min(minmax_value_arr[0]) - delta * 0.1, max(minmax_value_arr[1]) + delta * 0.1)
        ax1.set_xlabel(x_axis_label)
        ax1.grid(True, which='major', linestyle='--', linewidth=0.5)
        color_step = len(minmax_value_arr[0])
        ax2 = ax1.twinx()
        minmax_value_arr = [[], []]
        for right_axis_label in y_axis_label_dict[list(y_axis_label_dict.keys())[1]]:
            print(right_axis_label)
            ax2.plot(data.index, data[right_axis_label], label=right_axis_label, color=color_arr[color_step + len(minmax_value_arr[0])])  # Using black line for V_OC_TEG data
            ax2.set_ylabel(y_axis_label_dict['axis_labels'][1], color='black')
            ax2.tick_params(axis='y', labelcolor='black')
            # select adequate y scale range
            minmax_value_arr[0].append(data[right_axis_label].min())
            minmax_value_arr[1].append(data[right_axis_label].max())

        delta = max(minmax_value_arr[1]) - min(minmax_value_arr[0])
        ax2.set_ylim(min(minmax_value_arr[0]) - delta*0.1, max(minmax_value_arr[1]) + delta*0.1)
        # ax2.set_xlabel(x_axis_label)

    else:
        raise IOError('부적절한 입력 딕셔너리 구조')
    fig.legend(draggable=True)
    plt.show()
