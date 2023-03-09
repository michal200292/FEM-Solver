from math import pi
import matplotlib.pyplot as plt

# Parametry, które mogą być zmieniane przez użytkownika
# Stała G, można zmienić na inne wartości, żeby zobaczyć wykres inny niż linia prosta
G = 6.67 * pow(10, -11)

# Ilość punktów podziałowych
n = 1000


# Funkcja oblicza wartości, które znajdują się w macierzy kolumnowej
# Pod całką znajduje się funkcja bazowa
# h jest to szerokość podprzedziału, na jakie jest dzielony przedział (0, 3)
# x to współrzędna, w której znajduje się szpic funkcji, jako
# że funkcja p(x) nie zeruje się jedynie na przedziale (1, 2), wystarczy rozważać oblizanie całek.
# Jeśli funkcja bazowa należy do tego przedziału (1, 2) to obliczam jej wartość,
# w każdym przypadku jest to 4*pi*G * h^2
# Mogą zdarzyć się przypadki, że niezerowa część funkcji bazowej będzie tylko częściowo w przedziale (1, 2)
# Jednak niedokładność z tego wynikająca jest pomijalnie mała, tym bardziej dla dużych n

def calc_integral(x, h):
    return (1 <= x <= 2) * 4 * pi * G * pow(h, 2)


# Oblicza ostateczne wartości szukanej funkcji,
# a więc sumę wartości dla funkcji liniowej wprowadzonej ze względu na niezerowe warunki Dirichleta
# i wartości wynikającej z układu równań
# Wartości obliczam w takich punktach x, w których dokładnie jedna funkcja jest niezerowa
# Jej wartość wynosi tam h
def f(x, scalar, h):
    return (5 - x / 3) + h * scalar


# Funkcja rozwiązująca układ równań
'''
Macierz za każdym razem będzie wyglądała tak
[-2h h 0 0 0 .. 0]
[h -2h h 0 0 .. 0]
[0 h -2h h 0 .. 0]
[0 0 h -2h h .. 0]
...
[0 0 0 .. 0 h -2h]
Na początku można w jednej pętli sprowadzić ją do postaci macierzy trójkątnej
stosując technikę znaną z metody Gaussa
Po sprowadzeniu do macierzy trójkątnej, można już łatwo obliczać wartości kolejnych zmiennych, zaczynając
od zmiennej o numerze n
'''


def solve_linear_system(A, column, h):
    values = [0 for _ in range(n)]
    for j in range(n - 1):
        q = h / A[j]
        A[j + 1] -= q * h
        column[j + 1] -= q * column[j]

    values[-1] = column[-1] / A[-1]
    for j in range(n - 2, -1, -1):
        values[j] = (column[j] - h * values[j + 1]) / A[j]
    return values


def solve():
    # Wyznaczam szerokość podziału
    h = 3 / (n + 1)

    # Główna macierz układu
    A = [-2 * h for _ in range(n)]

    # Macierz kolumnowa
    column = [calc_integral((i + 1) * h, h) for i in range(n)]

    # Wyznaczenie wartości skalarów
    values = solve_linear_system(A, column, h)

    # Rysowanie wykresu
    mesh = [i * h for i in range(n + 2)]
    function_values = [0.0 for _ in range(n + 2)]
    function_values[0] = 5
    function_values[-1] = 4
    for i in range(n):
        function_values[i + 1] = f(mesh[i + 1], values[i], h)

    plt.plot(mesh, function_values)
    plt.show()


# Uruchomienie funkcji aproksymującej równanie
solve()
