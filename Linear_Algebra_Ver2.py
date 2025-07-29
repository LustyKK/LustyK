import math as m
import numpy as np
import sys 
import os 
from tqdm import tqdm
import time
import tkinter as tk
from tkinter import filedialog
from colorama import init,Fore

init(autoreset=True)

root = tk.Tk()
root.withdraw()

#===================================================================================
def print_menu():
    print("+------------------------------------------------------------------------------+")
    print("|                          Solve Linear Algebra 2.0                            |")
    print("|                          +/Phạm Tuấn Kiệt (LustyK)/+                         |")
    print("+------------------------------------------------------------------------------+")
    print("| 1. Đưa ma trận về bậc thang (chỉ áp dụng với 3x3)                            |")
    print("| 2. Tìm phương trình Lambda (Ax = Lambda*x*I)                                 |")
    print("| 3. Tính số chiều (Dim) của ma trận                                           |")
    print("| 4. Tính hạng (Rank) của ma trận                                              |")
    print("| 5. Kiểm tra độc lập tuyến tính (Independent)                                 |")
    print("| 6. Giải hệ phương trình tuyến tính                                           |")
    print("| 7. Tính độ dài,khoảng cách,góc,tích có hướng giữa 2 vector                   |")
    print("| 8. Xóa màn hình                                                              |")
    print("| 9. Exit                                                                      |")
    print("+------------------------------------------------------------------------------+")
#===================================================================================
def choose_menu():
    print("Lựa chọn cách nhập ma trận: ")
    print("+---------------------+")
    print("|     Press 1         |")
    print("| Nhập tay ma trận    |")
    print("+---------------------+")
    print("=============================================================")
    print("+----------------------+")
    print("|     Press 2          |")
    print("| Chọn tệp ma trận txt |")
    print("+----------------------+")    
#====================================================================================
def choose_menu_vector():
    print("Lựa chọn kích thước của vector: ")
    print("+------------------------------+")
    print("|            Press 2           |")
    print("|         Vector 2 chiều       |")
    print("+------------------------------+")
    print("=============================================================")
    print("+------------------------------+")
    print("|            Press 3           |")
    print("|         Vector 3 chiều       |")
    print("+------------------------------+")  
#====================================================================================
def read_matrix_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            matrix = [list(map(float, line.split())) for line in lines]
            return np.array(matrix)
    except FileNotFoundError:
        print("Không tìm thấy file.")
#===================================================================================
def input_matrix(so_hang, so_cot):
    ma_tran = []
    print(f"Nhập các phần tử cho ma trận {so_hang}x{so_cot}:")
    for i in range(so_hang):
        hang = []
        for j in range(so_cot):
            phan_tu = float(input(f"Nhập phần tử ở vị trí ({i+1},{j+1}): "))
            hang.append(phan_tu)
        ma_tran.append(hang)
    return np.array(ma_tran)
#===================================================================================
def output_matrix(ma_tran):
    for hang in ma_tran:
        print(hang)
#===================================================================================
def Echolon(x):
    try:
        b = x.reshape(3,3)
        if b.shape[0]==b.shape[1]:
            for i in range(b.shape[0]):
                for j in range(b.shape[1]):
                    if i==1:
                        b[i]=b[i-i][j]*b[i]-[b[i][j]*b[i-i]]
                        i=i+1
                        b[i]=b[i-i][j]*b[i]-[b[i][j]*b[i-i]]
                        j=j+1
                        b[i]=b[i]-[b[i-1]*(b[i][j]/b[i-1][j])]

        print(b)
    except ValueError:
        print("Đã xảy ra lỗi,tự ngắt kết nối")
        sys.exit()
#====================================================================================
def det(A):
    det = 0
    for i in range(3):
        det += A[0][i] * (A[1][(i + 1) % 3] * A[2][(i + 2) % 3] - A[1][(i + 2) % 3] * A[2][(i + 1) % 3])
    return det
#===================================================================================
def Find_Lambda_Equ(A):
    DetA= det(A)
    Det1=(A[0,0]*A[1,1])-(A[0,1]*A[1,0])
    Det2=(A[1,1]*A[2,2])-(A[1,2]*A[2,1])
    Det3=(A[0,0]*A[2,2])-(A[0,2]*A[2,0])
    solve_eqs(A,0)
    print("-x^3 +",A[0,0]+A[1,1]+A[2,2],"x^2 -",Det1+Det2+Det3,"x +",DetA)
#===================================================================================
def Rank_matrix(ma_tran):
  return np.linalg.matrix_rank(ma_tran)
#===================================================================================
def is_linearly_independent(A):
  matrix_np = np.array(A)
  rank = np.linalg.matrix_rank(matrix_np)
  return rank == matrix_np.shape[0]
#==================================================================================
def Dim_matrix(matrix):
    num_variables = int(input("Nhập số lượng biến: "))
    mat_np = np.array(matrix)
    dimension = mat_np.ndim
    return print("Số chiều của ma trận là: DIM = ",dimension) 
#===================================================================================
def solve_eqs(a,b):
  # Giải hệ phương trình
  try:
    x = np.linalg.solve(a,b)
    for i in range(n):
        print(f"x{i+1} = {x[i]}")
  except np.linalg.LinAlgError:
      print("Hệ phương trình vô nghiệm")
#===================================================================================
def input_2d_vector(so_cot):
    ma_tran = []
    for i in range(so_cot):
        phan_tu = float(input(f"Nhập tọa độ vector ở vị trí thứ ({i+1}): "))
        ma_tran.append(phan_tu)
    return ma_tran
#===================================================================================
def input_3d_vector(so_cot):
    ma_tran = []
    for i in range(so_cot):
        phan_tu = float(input(f"Nhập tọa độ vector ở vị trí thứ ({i+1}): "))
        ma_tran.append(phan_tu)
    return ma_tran
#===================================================================================
def vector_2D(point_A,point_B):
    vector_AB = (vector2[0] - vector1[0], vector2[1] - vector1[1])    
    dAB = m.sqrt(vector_AB[0] ** 2 + vector_AB[1] ** 2)
    Etymology_AB=np.dot(point_A,point_B)
    length_A = np.linalg.norm(point_A)
    length_B = np.linalg.norm(point_B)
    cos_theta = Etymology_AB / (length_A * length_B)
    theta = np.arccos(cos_theta)
    print("Độ dài của vector 1 là: ",length_A)
    print("+--------------------------------------------------------+")
    print("Độ dài của vector 2 là: ",length_B)
    print("+--------------------------------------------------------+")
    print("Tích có hướng của vector 1 và vector 2 là:", Etymology_AB)
    print("+--------------------------------------------------------+")
    print(f"Khoảng cách giữa 1 và 2 là: {dAB}")
    print("+--------------------------------------------------------+")
    print("Góc giữa A và B (rad):", theta)
    print("+--------------------------------------------------------+")
    print("Góc giữa A và B (độ):", np.degrees(theta))
    print("+--------------------------------------------------------+")
    if Etymology_AB == 0:
        print("2 vector này vuông góc (Orthogonal)")
    elif theta < 90:
        print("Góc giữa 2 vector là góc nhọn (Acute)")
    elif theta > 90:
        print("Góc giữa 2 vector là góc tù ( Obtuse)")
#===================================================================================
def vector_3D(point_C,point_D):
    vector_CD = (point_D[0] - point_C[0], point_D[1] - point_C[1], point_D[2] - point_C[2])
    dCD = m.sqrt(vector_CD[0] ** 2 + vector_CD[1] ** 2 + vector_CD[2] ** 2)
    Etymology_CD=np.dot(point_C,point_D)
    length_C = np.linalg.norm(point_C)
    length_D = np.linalg.norm(point_D)
    cos_theta2 = Etymology_CD / (length_C * length_D)
    theta2 = np.arccos(cos_theta2)   
    print("Tích có hướng của vector 1 và vector 2 là:", Etymology_CD)
    print("+--------------------------------------------------------+")
    print(f"Khoảng cách giữa 1 và 2 là: {dCD}")
    print("+--------------------------------------------------------+")
    print("Độ dài của vector 1 là: ",length_C)
    print("+--------------------------------------------------------+")
    print("Độ dài của vector 2 là: ",length_D)
    print("+--------------------------------------------------------+")
    print("Góc giữa C và D (rad):", theta2)
    print("+--------------------------------------------------------+")
    print("Góc giữa C và S (độ):", np.degrees(theta2))
    print("+--------------------------------------------------------+")
    if Etymology_CD == 0:
        print("2 vector này vuông góc (Orthogonal)")
    elif theta2 < 90:
        print("Góc giữa 2 vector là góc nhọn (Acute)")
    elif theta2 > 90:
        print("Góc giữa 2 vector là góc tù ( Obtuse)")
#===================================================================================
def result():
        print("=============================================================")
        for _ in tqdm(range(100), desc="Đang tính toán...", unit="iter"):
            time.sleep(0.009)
        print("=============================================================")
        print("+--------------------------------------------------------+")
        print("|                    KẾT QUẢ                             |")
        print("+--------------------------------------------------------+")
#==================================================================================
def start():
    while True:
        try:
            choose_menu()
            a=int(input())
            break
        except ValueError:  
            print("Lỗi cú pháp !!!")
#===================================================================================
    if a==1:
        so_hang = int(input("Nhập số hàng: "))
        so_cot = int(input("Nhập số cột: "))
        ma_tran=input_matrix(so_hang,so_cot)
    if a==2:
        file_path = filedialog.askopenfilename()  
        ma_tran = read_matrix_from_file(file_path)
    if a!=1 and a!=2:
        print("Đã xảy ra lỗi,chương trình tự ngắt....")
        sys.exit()

    return ma_tran
#==================================================================================
def loading_screen():
    # Biểu tượng loading
    loading_symbols = ['⣾', '⣽', '⣻', '⢿', '⡿', '⣟', '⣯', '⣷']

    # Số lượng biểu tượng loading
    num_loading_symbols = len(loading_symbols)

    # Thời gian giữa các biểu tượng loading (giây)
    delay_time = 0.05

    for _ in range(50):
        symbol = loading_symbols[_ % num_loading_symbols]
        sys.stdout.write(f"\r{Fore.MAGENTA}                         Đang khởi động chương trình... {Fore.LIGHTCYAN_EX}{symbol}                      ")
        sys.stdout.flush()
        time.sleep(delay_time)

    sys.stdout.write("\r                                                                                \n")
    sys.stdout.flush()
#==================================================================================
def Disconect_screen():
    # Biểu tượng loading
    loading_symbols = ['⣾', '⣽', '⣻', '⢿', '⡿', '⣟', '⣯', '⣷']

    # Số lượng biểu tượng loading
    num_loading_symbols = len(loading_symbols)

    # Thời gian giữa các biểu tượng loading (giây)
    delay_time = 0.1

    for _ in range(50):
        symbol = loading_symbols[_ % num_loading_symbols]
        sys.stdout.write(f"\r{Fore.MAGENTA}                                 Disconecting... {Fore.LIGHTCYAN_EX}{symbol}   ")
        sys.stdout.flush()
        time.sleep(delay_time)

    sys.stdout.write(f"\r{Fore.MAGENTA}                                     See you again !!!                             \n")
    sys.stdout.flush()
#====================================================================================
loading_screen()
while True:
    print_menu()
    try:
        cn = int(input())
    except ValueError:
        print("Lỗi !!! Vui lòng nhập lại")
#===================================================================================
    if cn==8:
        os.system('cls')
#===================================================================================
    if cn==9:
        os.system('cls')
        Disconect_screen()
        break
#===================================================================================
    if cn > 9 or cn < 1:
        print("Lỗi !!! Vui lòng nhập lại")
#===================================================================================
    if cn==1:
        ma_tran=start()
#===================================================================================
        print("ma trận đã được nhập là:")
        output_matrix(ma_tran)
        if ma_tran[0,0] == 0:
            print("+---------------------------------------------+")
            print("|                                             |")
            print("|       Không thể đưa ma trận về bậc thang    |")
            print("|                                             |")
            print("+---------------------------------------------+")            
        else:
            result()
            print("             Ma trận bậc thang là :                       ")
            print("+--------------------------------------------------------+")
            Echolon(ma_tran)
            print("+--------------------------------------------------------+")
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break
#===================================================================================
    if cn==2:
        ma_tran=start()
#===================================================================================
        print("ma trận đã được nhập là:")
        output_matrix(ma_tran)
        result()
        print("|                Phương trình Lambda là:                  |")
        print("+--------------------------------------------------------+")
        Find_Lambda_Equ(ma_tran)
        print("+--------------------------------------------------------+")
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break  
#===================================================================================   
    if cn==3:
        ma_tran=start()
        print("ma trận đã được nhập là:")
        output_matrix(ma_tran)
        result()
        print(Dim_matrix(ma_tran))
        print("+--------------------------------------------------------+")
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break
#===================================================================================                
    if cn==4:
        ma_tran=start()
        print("ma trận đã được nhập là:")
        output_matrix(ma_tran)
#===================================================================================
        ma_tran=start()
        print("Hạng của ma trận là: Rank A = ",Rank_matrix(ma_tran))
        print("+--------------------------------------------------------+")  
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break                        
#===================================================================================
    if cn==5:
        ma_tran=start()
        print("ma trận đã được nhập là:")
        output_matrix(ma_tran)
#===================================================================================
        result()
        if is_linearly_independent(ma_tran)== True:
            print("Ma trận độc lập tuyến tính (linearly independent)")
        else:
            print("Ma trận phụ thuộc tuyến tính (linearly dependent)")
        print("+--------------------------------------------------------+")
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break        
    if cn==6:
        n = int(input("Nhập số biến của hệ phương trình: "))
        a = np.zeros((n, n)) 
        b = np.zeros(n)
        for i in range(n):
            for j in range(n):
                a[i][j] = float(input(f"Nhập hệ số ma trận vị trí thứ dòng {i+1} cột {j+1}: "))
            b[i] = float(input(f"Nhập hệ số phương trình {i+1} = "))
#===================================================================================
        result()
        solve_eqs(a,b)   
        while True:
            print("Bấm B để trở về menu")
            b=str(input())
            if b=="b" or b=="B":
                break     
#====================================================================================
    if cn==7:
        choose_menu_vector()
        while True:
                while True:
                    try:
                        a=int(input())
                        break
                    except ValueError:
                        print("Sai cú pháp !!! nhập lại")
                if a==2:
                    so_cot = 2
                    print("=============================================================")
                    print("Nhập tọa độ vector 1")
                    print("=============================================================")
                    vector1 = input_2d_vector(so_cot)
                    vector1 = np.array(vector1)
                    print("=============================================================")
                    print("Nhập tọa độ vector 2")
                    print("=============================================================")
                    vector2 = input_2d_vector(so_cot)
                    vector2 = np.array(vector2)
                    result()
                    vector_2D(vector1,vector2)
                    print("+--------------------------------------------------------+")
                    while True:
                        try:
                            print("Bấm B để trở về menu")
                            b=str(input())
                            if b=="b" or b=="B":
                                break     
                        except ValueError:
                            print("Bấm B để trở về menu")
                    if b=="b" or b=="B":
                        break     
                if a==3:
                    so_cot = 3
                    print("=============================================================")
                    print("Nhập tọa độ vector 1")
                    print("=============================================================")
                    vector1 = input_3d_vector(so_cot)
                    vector1 = np.array(vector1)
                    print("=============================================================")
                    print("Nhập tọa độ vector 2")
                    print("=============================================================")
                    vector2 = input_3d_vector(so_cot)
                    vector2 = np.array(vector2)
                    result()
                    vector_3D(vector1,vector2)
                    print("+--------------------------------------------------------+")
                    while True:
                        try:
                            print("Bấm B để trở về menu")
                            b=str(input())
                            if b=="b" or b=="B":
                                break     
                        except ValueError:
                            print("Bấm B để trở về menu")
                    if b=="b" or b=="B":
                        break                         
                if a!= 2 or a!=3:
                    print("Sai cú pháp !!! nhập lại")