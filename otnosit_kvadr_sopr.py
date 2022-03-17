# Библиотеки
import numpy as np
import sympy as sp
import sympy.abc as abc

# Текстовый документ
Res = open("Res_part_2.2.txt", "w")

# Исходные данные
L = 90         # Длина, м
B = 12.35       # Ширина, м
T = 3.782       # Осадка, м
H = 6.4       # Высота борта, м
h_0 = 0.49  # Метацентрическая высота, м
r = 3.516  # Поперечный метацентрический радиус, м
alpha = 0.75  # Коэффициент полноты ватерлинии, -
delta = 0.65  # Коэффициент общей полноты, -
chi = delta / alpha  # Коэффициент вертикальной полноты, -
rho = 1.025  # Плотность воды, т/м^3
g = 9.81  # Ускорение свободного паденияя. м/с^2
V = 2732.40045  # Объёмное водоизмещение, м^3
D = rho*g*V         # Массовое водоизмещение, т
J_xx = D/g*(B**2*alpha**2/(11.4*delta) + H**2/12)               # Момент инерции массы корпуса, т*м^2
rho_x = (J_xx/(rho*V))**(1/2)                                   # Приведенный радиус инерции массы корпуса, м
lambda_44 = J_xx/(0.28 + 1.8*rho_x/(alpha*B*(1 + 1/6*B/T)))     # Момент инерции присоединенных масс, т*м^2
rho_x_ = ((J_xx + lambda_44)/(rho*V))**(1/2)                    # Приведенный радиус инерции массы корпуса с учётом момента инерции присоединённых масс, м
q = J_xx/(J_xx + lambda_44)
omega_lin = (D*h_0/(J_xx + lambda_44))**(1/2)                   # Собственная частота бортовой качки по линейной теории, 1/с
w = 0.24 + 1.42*h_0/B*(alpha*B*(1 + 1/6*B/T)/(4*rho_x_))**2     # Безразмерный коэффициент квадратичного сопротивления
_2_nu_theta = 0.3*omega_lin*w**(1/2)                            # Коэффициент демпфирования, 1/с
nu_theta = 0.15*omega_lin*w**(1/2)
R = chi*((B*T*chi*r/h_0)**(1/2)/(2*np.pi*g))**(1/2)

# Коэффициенты аппроксимации ДСО и ДДО
A_l = 0.0000000527
B_l = -0.000013
C_l = 0.000437
D_l = 0.00919

A_d = -0.0000000321
B_d = 0.00000156
C_d = 0.0000828
D_d = 0.00025

#
theta_max_ = sp.solve(A_l * abc.x ** 4 + B_l * abc.x ** 3 + C_l * abc.x ** 2 + D_l * abc.x, dict=True)
theta_max_ = sp.re(theta_max_[2][abc.x])

# Способ В.Г. Власова
def Vlasov():
    with open("Res_part_2.2.txt", "a") as Res:
        print("ОТНОСИТЕЛЬНЫЕ КООРДИНАТЫ. КВАДРАТИЧНОЕ СОПРОТИВЛЕНИЕ", file=Res)
        print("theta\tomega_0\t\tomega_1\t\tomega_2", file=Res)

    theta_delta = 2.5
    theta = theta_delta
    theta_max = 50

    while theta <= theta_max:
        l = A_l * theta ** 4 + B_l * theta ** 3 + C_l * theta ** 2 + D_l * theta
        d = A_d * theta ** 4 + B_d * theta ** 3 + C_d * theta ** 2 + D_d * theta

        theta_rad = np.radians(theta)

        T = 2 * np.pi * ((J_xx + lambda_44) * theta_rad / (D * (l / 2 + d / theta_rad))) ** (1 / 2)
        omega = 2 * np.pi / T

        # Относительные координаты. Квадратичное сопротивление
        _lambda = 2 * np.pi * g / omega ** 2
        h = 0.17 * _lambda ** 0.75
        a_w = h / 2
        alpha_0 = omega ** 2 / g * a_w

        aleph = np.exp(-4.2 * (R * omega) ** 2)
        alpha_m = aleph * alpha_0

        _A_ = q ** 2 * alpha_m ** 2 / theta_rad ** 2 - 64 / (9 * np.pi ** 2) * w ** 2 * theta_rad ** 2
        _B_ = 1 + 64 / (9 * np.pi ** 2) * w ** 2 * theta_rad ** 2 - q ** 2 * alpha_m ** 2 / theta_rad ** 2

        if _A_ > 0:
            if omega ** 2 * (1 + _A_ ** (1 / 2)) / _B_ > 0:
                omega_1 = (omega ** 2 * (1 + _A_ ** (1 / 2)) / _B_) ** (1 / 2)
            else:
                omega_1 = "-"

            if omega ** 2 * (1 - _A_ ** (1 / 2)) / _B_ > 0:
                omega_2 = (omega ** 2 * (1 - _A_ ** (1 / 2)) / _B_) ** (1 / 2)
            else:
                omega_2 = "-"

            theta_m = np.degrees((3 * np.pi * q * alpha_m / (8 * w)) ** (1 / 2))

            with open("Res_part_2.2.txt", "a") as Res:
                print(theta, theta_m, omega_1, omega_2, sep="\t", file=Res)

        else:
            omega_1 = "-"
            omega_2 = "-"

            theta_m = np.degrees((3 * np.pi * q * alpha_m / (8 * w)) ** (1 / 2))

            with open("Res_part_2.2.txt", "a") as Res:
                print(theta, theta_m, omega_1, omega_2, sep="\t", file=Res)

        theta = theta + theta_delta


Vlasov()
