import numpy as np
import sympy as sp
import sympy.abc as abc

# Текстовый документ
Res = open("Res_part_1.txt", "w")

# Исходные данные
L = 130         # Длина, м
B = 16.25       # Ширина, м
T = 6.450       # Осадка, м
H = 9.6105       # Высота борта, м
h_0 = 0.402907  # Метацентрическая высота, м
r = 3.298555  # Поперечный метацентрический радиус, м
alpha = 0.704  # Коэффициент полноты ватерлинии, -
delta = 0.567  # Коэффициент общей полноты, -
chi = delta / alpha  # Коэффициент вертикальной полноты, -
rho = 1.025  # Плотность воды, т/м^3
g = 9.81  # Ускорение свободного паденияя. м/с^2
V = 8049.09  # Объёмное водоизмещение, м^3
D = rho * g * V  # Массовое водоизмещение, т
J_xx = D / g * (B ** 2 * alpha ** 2 / (11.4 * delta) + H ** 2 / 12)  # Момент инерции массы корпуса, т*м^2
rho_x = (J_xx / (rho * V)) ** (1 / 2)  # Приведенный радиус инерции массы корпуса, м
lambda_44 = J_xx / (0.28 + 1.8 * rho_x / (alpha * B * (1 + 1 / 6 * B / T)))  # Момент инерции присоединенных масс, т*м^2
rho_x_ = ((J_xx + lambda_44) / (rho * V)) ** (
            1 / 2)  # Приведенный радиус инерции массы корпуса с учётом момента инерции присоединённых масс, м
q = J_xx / (J_xx + lambda_44)
omega_lin = (D * h_0 / (J_xx + lambda_44)) ** (1 / 2)  # Собственная частота бортовой качки по линейной теории, 1/с
w = 0.24 + 1.42 * h_0 / B * (
            alpha * B * (1 + 1 / 6 * B / T) / (4 * rho_x_)) ** 2  # Безразмерный коэффициент квадратичного сопротивления
_2_nu_theta = 0.3 * omega_lin * w ** (1 / 2)  # Коэффициент демпфирования, 1/с
nu_theta = 0.15 * omega_lin * w ** (1 / 2)
R = chi * ((B * T * chi * r / h_0) ** (1 / 2) / (2 * np.pi * g)) ** (1 / 2)

# Коэффициенты аппроксимации ДСО и ДДО
A_l = 0.0000002
B_l = -0.00003
C_l = 0.0007
D_l = 0.0125

A_d = -0.00000002
B_d = -0.0000009
C_d = 0.0002
D_d = -0.0004

#
theta_max_ = sp.solve(A_l * abc.x ** 4 + B_l * abc.x ** 3 + C_l * abc.x ** 2 + D_l * abc.x, dict=True)
theta_max_ = sp.re(theta_max_[2][abc.x])

# Способ В.Г. Власова
def Vlasov():
    with open("Res_part_1.txt", "a", encoding = "utf-8") as Res:
        print("СПОСОБ В.Г. ВЛАСОВА", file=Res)
        print("theta\tT\tomega", file=Res)
        print(0, 2 * np.pi / omega_lin, omega_lin, sep="\t", file=Res)
        Res.close()

    theta_delta = 10
    theta = theta_delta
    theta_max = 60

    while theta <= theta_max:
        l = A_l * theta ** 4 + B_l * theta ** 3 + C_l * theta ** 2 + D_l * theta
        d = A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta

        theta_rad = np.radians(theta)

        T = 2 * np.pi * ((J_xx + lambda_44) * theta_rad / (D * (l / 2 + d / theta_rad))) ** (1 / 2)
        omega = 2 * np.pi / T

        with open("Res_part_1.txt", "a") as Res:
            print(theta, T, omega, sep="\t", file=Res)
            Res.close()

        theta = theta + theta_delta

    with open("Res_part_1.txt", "a") as Res:
        print(theta_max_, "-", 0, sep="\t", file=Res)
        Res.close()

# Способ А.Б. Карпова
def Karpov():
    with open("Res_part_1.txt", "a", encoding = "utf-8") as Res:
        print("\nСПОСОБ А.Б. КАРПОВА", file=Res)
        print("theta\tT\tomega", file=Res)
        print(0, 2 * np.pi / omega_lin, omega_lin, sep="\t", file=Res)
        Res.close()

    theta_delta = 10
    theta = theta_delta
    theta_max = 60

    while theta <= theta_max:
        theta_i = 0
        theta_i_delta = 2.5

        d = A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta

        sum_1 = 0

        while theta_i < theta:
            d_i = A_d * theta_i ** 4 + B_d * theta_i ** 3 + C_d * theta_i ** 2 + D_d * theta_i

            sum_1 = sum_1 + 1 / (d - d_i) ** (1 / 2)

            theta_i = theta_i + theta_i_delta

        d_1 = A_d * (theta - theta_i_delta) ** 4 + B_d * (theta - theta_i_delta) ** 3 + C_d * (
                    theta - theta_i_delta) ** 2 + D_d * (theta - theta_i_delta)

        P = 1 / 2 * (1 / (d - d_1) ** (1 / 2) + 1 / (d - d_i) ** (1 / 2))

        theta_i_delta_rad = np.radians(theta_i_delta)

        p = 2 / (d - d_1) ** (1 / 2)

        T = 4 * ((J_xx + lambda_44) / (2 * D)) ** (1 / 2) * theta_i_delta_rad * (sum_1 - P + p)

        omega = 2 * np.pi /T
        with open("Res_part_1.txt", "a") as Res:
            print(theta, T, omega, sep="\t", file=Res)
            Res.close()

        theta = theta + theta_delta

    with open("Res_part_1.txt", "a") as Res:
        print(theta_max_, "-", 0, sep="\t", file=Res)
        Res.close()


# Способ В.Г. Сизова
def Sizov():
    with open("Res_part_1.txt", "a", encoding = "utf-8") as Res:
        print("\nСПОСОБ В.Г. СИЗОВА", file=Res)
        print("theta\tT\tomega", file=Res)
        print(0, 2 * np.pi / omega_lin, omega_lin, sep="\t", file=Res)
        Res.close()

    theta_delta = 10
    theta = theta_delta
    theta_max = 60

    while theta <= theta_max:
        phi_delta = 10
        phi = phi_delta
        phi_max = 90

        sum_1 = 0

        d = A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta

        while phi <= phi_max:
            phi_rad = np.radians(phi)

            theta_i = theta * np.cos(phi_rad)

            d_i = A_d * theta_i ** 4 + B_d * theta_i ** 3 + C_d * theta_i ** 2 + D_d * theta_i

            sum_1 = sum_1 + np.sin(phi_rad) / (d - d_i) ** (1 / 2)

            phi = phi + phi_delta

        l = A_l * theta ** 4 + B_l * theta ** 3 + C_l * theta ** 2 + D_l * theta

        theta_rad = np.radians(theta)

        P = 1 / 2 * ((2 / (theta_rad * l)) ** (1 / 2) + np.sin(phi_rad) / (d - d_i) ** (1 / 2))

        phi_delta_rad = np.radians(phi_delta)

        sum_2 = phi_delta_rad * (((2 / (theta_rad * l)) ** (1 / 2) + sum_1) - P)

        T = (8 * (J_xx + lambda_44) / D) ** (1 / 2) * theta_rad * sum_2

        omega = 2 * np.pi / T

        with open("Res_part_1.txt", "a") as Res:
            print(theta, T, omega, sep="\t", file=Res)
            Res.close()

        theta = theta + theta_delta

    with open("Res_part_1.txt", "a") as Res:
        print(theta_max_, "-", 0, sep="\t", file=Res)
        Res.close()

        # Способ Г.Е. Павленко
def Pavlenko():
    with open("Res_part_1.txt", "a", encoding = "utf-8") as Res:
        print("\nСПОСОБ Г.Е. ПАВЛЕНКО", file=Res)
        print("theta\tT\tomega", file=Res)
        print(0, 2 * np.pi / omega_lin, omega_lin, sep="\t", file=Res)
        Res.close()

    theta_delta = 10
    theta = theta_delta
    theta_max = 60

    while theta <= theta_max:
        xi_delta = 10
        xi = xi_delta
        xi_max = 90

        sum_1 = 0

        while xi <= xi_max:
            xi_rad = np.radians(xi)

            d = (A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta) * (np.sin(xi_rad)) ** 2

            theta_x = sp.solve(A_d * abc.x ** 4 + B_d * abc.x ** 3 + C_d * abc.x ** 2 + D_d * abc.x - d, dict=True)

            l = A_l * sp.re(theta_x[2][abc.x]) ** 4 + B_l * sp.re(theta_x[2][abc.x]) ** 3 + C_l * sp.re(
                theta_x[2][abc.x]) ** 2 + D_l * sp.re(theta_x[2][abc.x])

            sum_1 = sum_1 + np.sin(xi_rad) / l

            xi = xi + xi_delta

        d = A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta

        xi_rad = np.radians(xi)

        P = 1 / 2 * (1 / (2 * d * h_0) ** (1 / 2) + np.sin(xi_rad) / l)

        xi_delta_rad = np.radians(xi_delta)

        sum_2 = xi_delta_rad * ((1 / (2 * d * h_0) ** (1 / 2) + sum_1) - P)

        T = 8 * ((J_xx + lambda_44) / (2 * D) * d) ** (1 / 2) * sum_2

        omega = 2 * np.pi / T

        with open("Res_part_1.txt", "a") as Res:
            print(theta, T, omega, sep="\t", file=Res)
            Res.close()

        theta = theta + theta_delta

    with open("Res_part_1.txt", "a") as Res:
        print(theta_max_, "-", 0, sep="\t", file=Res)
        Res.close()

Vlasov()
Karpov()
Sizov()
Pavlenko()

print(D, J_xx, rho_x, lambda_44, rho_x_, w, _2_nu_theta, q, R, omega_lin, rho, g)