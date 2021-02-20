import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import operator, math
import numpy
from numpy import sqrt, dot, cross
from numpy.linalg import norm

SOUND_SPEED_IN_SEA_WATER = 1.5216
POLUS_EARTH_RADIUS = 6357
ECVATOR_EARTH_RADIUS = 6378
F_EARTH = (ECVATOR_EARTH_RADIUS - POLUS_EARTH_RADIUS) / ECVATOR_EARTH_RADIUS


class CanCopyLabel(ttk.Entry):
    def __init__(self, master=None, text=''):
        super().__init__(master)
        self.insert(tk.END, text)
        self.config(state='readonly')
        self.bind('<Control-c>', self.clip_copy)

    def clip_copy(self, event=None):
        tk._default_root.clipboard_clear()
        tk._default_root.clipboard_append(self.get())


import fractions


def calc_slu2(f, s):
    cf = f[:]
    zero = s[0]
    zerof = f[0]
    for i in range(len(f)):
        f[i] = f[i] * zero
        s[i] = s[i] * zerof
    try:
        Y = (f[2] - s[2]) / (f[1] - s[1])
        X = (f[2] - f[1] * Y) / f[0]
    except Exception as E:
        raise ValueError('error in calc system of linear alg - ' + str(E))

    return X, Y


def calc_slu3(f, s, t):
    zerof = f[0]
    zeros = s[0]
    zerot = t[0]
    fc = f[:]
    for i in range(4):
        fc[i] = fc[i] * zeros
        s[i] = s[i] * zerof
    row1 = []
    for i in range(1, 4):
        row1.append(s[i] - fc[i])
    fc = f[:]
    for i in range(4):
        fc[i] = fc[i] * zerot
        t[i] = t[i] * zerof
    row2 = []
    for i in range(1, 4):
        row2.append(t[i] - fc[i])
    Y, Z = calc_slu2(row1, row2)
    X = (f[-1] - Z * f[2] - Y * f[1]) / f[0]
    return X, Y, Z


def initSpher(a, f):
    b = a * (1. - f)
    c = a / (1. - f)
    e2 = f * (2. - f)
    e12 = e2 / (1. - e2)
    return (b, c, e2, e12)


def N(B):
    return ECVATOR_EARTH_RADIUS ** 2 / math.sqrt(
        ECVATOR_EARTH_RADIUS ** 2 * math.cos(math.radians(B)) ** 2 + POLUS_EARTH_RADIUS ** 2 * math.sin(
            math.radians(B)))


III = 0

main_win = tk.Tk()
main_win.title('Определение локации донной станции')

root = ttk.Frame(main_win)  # Для того чтобы фон окошка не был убогим, и соответсвовал стилям ttk
root.pack(fill=tk.BOTH, expand=1)

edit_frm = ttk.Frame(root)  # Тулбар с кнопками
edit_frm.pack(side=tk.TOP, fill='x')


#  Google Координаты

def cos(alpha): return math.cos(math.radians(alpha))


def sin(alpha): return math.sin(math.radians(alpha))


def to_decart_system(lat, lon):
    h = 0
    a = ECVATOR_EARTH_RADIUS
    f = F_EARTH
    b, c, e2, e12 = initSpher(a, f)
    cos_lat = math.cos(math.radians(lat))
    n = c / math.sqrt(1. + e12 * cos_lat ** 2)
    p = (n + h) * cos_lat
    x = p * math.cos(math.radians(lon))
    y = p * math.sin(math.radians(lon))
    z = (n + h - e2 * n) * math.sin(math.radians(lat))
    return x, y, z


def to_geodezic_system(x, y, z):
    a = ECVATOR_EARTH_RADIUS
    f = F_EARTH

    b, c, e2, e12 = initSpher(a, f)
    p = math.hypot(x, y)
    t = z / p * (1. + e12 * b / math.sqrt(x ** 2 + y ** 2 + z ** 2))
    fi = math.degrees(math.atan((1 - f) * t))
    print('t:', math.degrees(math.atan(t)), 'fi:', fi)
    fi = math.degrees(math.atan(t))
    tg_B_1 = (z + e12 * POLUS_EARTH_RADIUS * math.sin(fi) ** 3) / (p - e2 * ECVATOR_EARTH_RADIUS * math.cos(fi) ** 3)
    print('tg_B_1 - B:', math.degrees(math.atan(tg_B_1)))
    for i in range(3):
        fi = tg_B_1
        tg_B_1 = (z + e12 * POLUS_EARTH_RADIUS * math.sin(fi) ** 3) / (
                    p - e2 * ECVATOR_EARTH_RADIUS * math.cos(fi) ** 3)
        print('tg_B_1 - B:', math.degrees(math.atan(tg_B_1)))
    return math.degrees(math.atan(tg_B_1)), math.degrees(math.atan2(math.radians(y), math.radians(x)))


# Find the intersection of three spheres                  
# P1,P2,P3 are the centers, r1,r2,r3 are the radii       
# Implementaton based on Wikipedia Trilateration article. 
def trilaterate(P1, P2, P3, r1, r2, r3):
    P1 = numpy.array(P1)
    P2 = numpy.array(P2)
    P3 = numpy.array(P3)
    temp1 = P2 - P1
    e_x = temp1 / norm(temp1)
    temp2 = P3 - P1
    i = dot(e_x, temp2)
    temp3 = temp2 - i * e_x
    e_y = temp3 / norm(temp3)
    e_z = cross(e_x, e_y)
    d = norm(P2 - P1)
    j = dot(e_y, temp2)
    x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j)
    temp4 = r1 * r1 - x * x - y * y
    z = sqrt(temp4)
    p_12_a = P1 + x * e_x + y * e_y + z * e_z
    p_12_b = P1 + x * e_x + y * e_y - z * e_z
    return p_12_a, p_12_b


def trilaterate2(P1, P2, P3, r1, r2, r3):  # Решение системы уравнений
    x1, y1, z1 = P1
    x2, y2, z2 = P2
    x3, y3, z3 = P3
    first_row = [
        2 * (x2 - x1), 2 * (y2 - y1), 2 * (z2 - z1),
        r1 ** 2 - r2 ** 2 - x1 ** 2 - y1 ** 2 - z1 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    ]
    second_row = [
        2 * (x3 - x1), 2 * (y3 - y1), 2 * (z3 - z1),
        r1 ** 2 - r3 ** 2 - x1 ** 2 - y1 ** 2 - z1 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    ]
    third_row = [
        2 * (x2 - x3), 2 * (y2 - y3), 2 * (z2 - z3),
        r3 ** 2 - r2 ** 2 - x3 ** 2 - y3 ** 2 - z3 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    ]
    print(first_row, second_row, third_row, sep='\n')
    X, Y, Z = calc_slu3(first_row, second_row, third_row)
    return X, Y, Z


def from_gcoods(coods):
    if '"' in coods:  # Format: 41°24'12.2"N 2°10'26.5"E
        SHIROTA = coods.split(' ')[0]
        DEEPOTA = coods.split(' ')[-1]
        s_SPIN = SHIROTA.split('"')[-1]
        d_SPIN = DEEPOTA.split('"')[-1]
        s_GRAD = SHIROTA.split('°')[0]
        s_MINUT = SHIROTA.split('°')[1].split("'")[0]
        s_SECOND = SHIROTA.split('°')[1].split("'")[1].split('"')[0]
        d_GRAD = DEEPOTA.split('°')[0]
        d_MINUT = DEEPOTA.split('°')[1].split("'")[0]
        d_SECOND = DEEPOTA.split('°')[1].split("'")[1].split('"')[0]
        try:
            s_GRAD = float(s_GRAD)
            d_GRAD = float(d_GRAD)
            s_MINUT = float(s_MINUT)
            d_MINUT = float(d_MINUT)
            s_SECOND = float(s_SECOND)
            d_SECOND = float(d_SECOND)
        except:
            messagebox.showerror('Определение локации донной станции', 'Вы указали невалидные координаты!')
            return
        sch = s_SECOND + s_MINUT * 60 + s_GRAD * 3600 if s_SPIN == 'N' else - s_SECOND + s_MINUT * 60 + s_GRAD * 3600
        dep = d_SECOND + d_MINUT * 60 + d_GRAD * 3600 if d_SPIN == 'E' else - d_SECOND + d_MINUT * 60 + d_GRAD * 3600
        return sch / 3600, dep / 3600
    if ', ' in coods:  # Format: 41.40338, 2.17403
        SHIROTA = coods.split(', ')[0]
        DEEPOTA = coods.split(', ')[-1]
        try:
            sch = float(SHIROTA)
            dep = float(DEEPOTA)
        except:
            messagebox.showerror('Определение локации донной станции', 'Вы указали невалидные координаты!')
            return
        return sch, dep


def to_gcoods(x, y):
    return str(x) + ', ' + str(y)


###->


def ask_point(event=None):
    def yes():
        global III
        treev.insert('', tk.END, values=(str(III), dte.get(), to_gcoods(*from_gcoods(cgg.get()))))
        III += 1
        ask_win.destroy()

    ask_win = tk.Toplevel()
    ask_win.title('Создать измерение...')
    ask = ttk.Frame(ask_win)
    ask.pack(fill=tk.BOTH, expand=1)
    ttk.Label(ask, text='Google Координаты: ').grid(row=0, column=0, padx=5, pady=5)
    ttk.Label(ask, text='Δt - Дельта t - Задержка сигнала: ').grid(row=1, column=0, padx=5, pady=5)
    cgg = ttk.Entry(ask, width=30)
    cgg.grid(row=0, column=1, padx=(5, 0), pady=5)

    def clipp(event=None): cgg.insert(tk.INSERT, main_win.clipboard_get())

    cgg.bind('<Double-Button-1>', clipp)
    dte = ttk.Entry(ask, width=15)
    dte.grid(row=1, column=1, columnspan=1, padx=5, pady=5)
    ttk.Button(ask, text='Ок', width=50, command=yes).grid(row=2, padx=15, pady=15, columnspan=3)
    ask_win.transient(tk._default_root)
    ask_win.resizable(0, 0)
    ask_win.tkraise()
    ask_win.grab_set()

    ask_win.wait_window()  # ждем окошечко


def del_point(event=None):
    if treev.selection():
        if messagebox.askyesno('Удаление измерения', 'Вы точно хотите удалить измерение?'):
            treev.delete(treev.selection())


def calc_point(event=None):
    points = [treev.item(x)['values'][1:] for x in treev.get_children()]
    points = [[float(s) / 2 * SOUND_SPEED_IN_SEA_WATER, list(from_gcoods(coods))] for s, coods in points]
    if len(points) < 3:
        messagebox.showerror('Ошибка', 'Измерений должно быть не меньше 3!')
        return
    print(points)
    for i in range(len(points)):
        points[i] = [points[i][0], list(to_decart_system(*points[i][1]))]
    points = sorted(points, key=operator.itemgetter(0))
    points = points[:3]
    result = trilaterate(points[0][1], points[1][1], points[2][1], points[0][0], points[1][0], points[2][0])
    it = 0
    while numpy.isnan(numpy.sum(result[0])):
        points[0][0] *= 1.05
        points[1][0] *= 1.05
        points[2][0] *= 1.05
        result = trilaterate(points[0][1], points[1][1], points[2][1], points[0][0], points[1][0], points[2][0])
        print(result)
        it += 1
        if it > 179:
            messagebox.showerror('Ошибка!', 'Ваши данные не верны! Система не находит ответ!')
            return
    a, b = result
    a = list(a)
    b = list(b)
    if a[2] > b[2]:
        result = b
    else:
        result = a

    result = to_geodezic_system(*result)
    print(result)
    topka = tk.Toplevel()

    win = ttk.Frame(topka)
    win.pack(fill='both', expand=1)
    ttk.Label(win, text='Координаты точки ~ ').grid(row=0, column=0, padx=(10, 5), pady=15)
    CanCopyLabel(win, text=to_gcoods(round(result[0], 6), round(result[1], 6))).grid(row=0, column=1)
    topka.resizable(0, 0)
    topka.transient(main_win)


add_butt = ttk.Button(edit_frm, text='Добавить измерение', width=30, command=ask_point)
add_butt.grid(row=0, column=0, pady=2, padx=2)

del_butt = ttk.Button(edit_frm, text='Удалить измерение', width=30, command=del_point)
del_butt.grid(row=0, column=2, pady=2, padx=2)

del_butt = ttk.Button(edit_frm, text='Подсчитать результаты', width=30, command=calc_point)
del_butt.grid(row=0, column=3, pady=2, padx=2)

columns = ("#1", "#2", "#3")
treev = ttk.Treeview(root, show="headings", columns=columns)
treev.heading("#1", text="Номер")
treev.heading("#2", text="Δt")
treev.heading("#3", text="Google Коорд.")
treev.pack(expand=1, side=tk.BOTTOM, fill=tk.BOTH)

lat, lon = 35, 45
to = to_decart_system(lat, lon)
print(to)
from_ = to_geodezic_system(*to)
print(from_)

main_win.mainloop()
