import parser

import numpy as np
import matplotlib.pyplot as plt
import sympy
from math import factorial


def calc_lagrange(function, x):
    result = 0
    for i in range(len(function)):
        numerator = 1
        denominator = 1
        for j in range(len(function)):
            if i != j:
                numerator *= x - function[j][0]
                denominator *= function[i][0] - function[j][0]
        result += function[i][1] * numerator / denominator
    return result


def calc_newton(function, x):
    n = len(function)
    h = function[1][0] - function[0][0]
    finite_diff = [[0] * n for _ in range(n)]
    for i in range(n):
        finite_diff[i][0] = function[i][1]
    for i in range(1, n):
        for j in range(n - i):
            finite_diff[j][i] = finite_diff[j + 1][i - 1] - finite_diff[j][i - 1]
    if x <= function[int(n / 2)][0]:
        x0_ind = n - 1
        for i in range(n):
            if x <= function[i][0]:
                x0_ind = i - 1
                break
        x0_ind = max(x0_ind, 0)
        t = (x - function[x0_ind][0]) / h
        result = finite_diff[x0_ind][0]
        for i in range(1, n):
            t_mult = t
            for j in range(1, i):
                t_mult *= t - j
            result += (t_mult / factorial(i)) * finite_diff[x0_ind][i]
    else:
        t = (x - function[n - 1][0]) / h
        result = finite_diff[n - 1][0]
        for i in range(1, n):
            t_mult = t
            for j in range(1, i):
                t_mult *= t + j
            result += (t_mult / factorial(i)) * finite_diff[n - i - 1][i]
    return result


if __name__ == '__main__':
    print("Выберите тип ввода дланных\n Функция - 1\n Набор точек - 2:")
    while True:
        try:
            input_type = int(input())
            if input_type not in [1, 2]:
                raise Exception
            break
        except Exception:
            print("Введенное число должно быть 1 или 2. Повторите ввод:")
    if input_type == 1:
        print("Введите выражение:")
        while True:
            try:
                func = input()
                x = sympy.symbols('x')
                sympy.sympify(func).evalf(subs={x: 1})
                break
            except Exception:
                print("Некорректное выражение. Повторите ввод:")
        print("Введите количество узлов:")
        while True:
            try:
                num = int(input())
                if num < 2:
                    raise Exception
                break
            except Exception:
                print("Некорректное число узлов. Число узлов должно быть целым >= 2. Повторите ввод:")
        print("Введите границы отрезка интерполяции:")
        while True:
            try:
                a, b = map(float, input().split(" "))
                if a > b:
                    c = a
                    a = b
                    b = c
                break
            except Exception:
                print("Некорректные границы. Повторите ввод:")
        function = []
        h = (b - a) / (num - 1)
        for i in range(num):
            function.append((a, sympy.sympify(func).evalf(subs={x: a})))
            a += h
    else:
        print("Введите количество точек:")
        while True:
            try:
                num = int(input())
                if num < 2:
                    raise Exception
                break
            except Exception:
                print("Некорректное число точек. Число узлов должно быть целым >= 2. Повторите ввод:")
        cnt = 0
        function = []
        print("Введите точки в формате x y :")
        while True:
            try:
                a, b = map(float, input().split(" "))
                function.append((a, b))
                cnt += 1
                if cnt == num:
                    break
            except Exception:
                print("Некорректная точка. Повторите ввод:")
    print("Введите координату x точки для рассчета значения:")
    while True:
        try:
            x_coord = float(input())
            break
        except Exception:
            print("Некорректная точка. Повторите ввод:")

    x = np.array([point[0] for point in function])
    y = np.array([point[1] for point in function])
    plot_x = np.linspace(np.min(x), np.max(x))
    print("Значение функции в точке по Лагранжу: " + str(calc_lagrange(function, x_coord)))
    plot_y_lagrange = [calc_lagrange(function, x) for x in plot_x]
    print("Значение функции в точке по Ньютону: " + str(calc_newton(function, x_coord)))
    plot_y_newton = [calc_newton(function, x) for x in plot_x]
    ax = plt.gca()
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_title('Интерполяция по формуле Лагранжа')
    plt.plot(x, y, 'o', plot_x, plot_y_lagrange)
    plt.show()
    ax = plt.gca()
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_title('Интерполяция по формуле Ньютона')
    plt.plot(x, y, 'o', plot_x, plot_y_newton)
    plt.show()
