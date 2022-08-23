# -*- coding: utf-8 -*-
"""
Editor Spyder

Este é um arquivo de script temporário.
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def EvolucaoItemA(N, M, dt, dx):
    Mat = []
    for k in range (M+1):
        Mat.append([0]*(N+1))
    for i in range (N+1):
        Mat[0][i] = ((i*dx)**2)*(1-i*dx)**2
    for k in range (M):
        for i in range (1,N):
            f = (10*np.cos(10*k*dt)*(i*dx)**2)*(1-i*dx)**2 - (1 + np.sin(10*k*dt))*(12*(i*dx)**2 - 12*(i*dx) +2)
            Mat[k+1][i] = Mat[k][i] + dt*(((Mat[k][i-1] - 2*Mat[k][i] + Mat[k][i+1])/dx**2) + f)
    return (Mat)
def EvolucaoItemB(N, M, dt, dx):
    Mat = []
    for k in range (M+1):
        Mat.append([0]*(N+1))
    for i in range (N+1):
        Mat[0][i] = (np.exp(-i*dx))
    for k in range (M+1):
        Mat[k][0] = (np.exp(k*dt))
        Mat[k][N] = (np.exp(k*dt - 1))*np.cos(5*k*dt)
    for k in range (M):
        for i in range (1,N):
            t = k*dt
            x = i*dx
            f = 25*(np.exp(t-x))*(t**2)*np.cos(5*t*x)-10*(np.exp(t-x))*t*np.sin(5*t*x)-5*((np.exp(t-x))*x*np.sin(5*t*x))
            Mat[k+1][i] = Mat[k][i] + dt*(((Mat[k][i-1] - 2*Mat[k][i] + Mat[k][i+1])/dx**2) + f)
    return (Mat)
def EvolucaoItemC(N, M, dt, dx):
    Mat = []
    for k in range (M+1):
        Mat.append([0]*(N+1))
    for k in range (M):
        for i in range (1,int(0.25/dx - 1/2)):
            f = 0
            Mat[k+1][i] = Mat[k][i] + dt*(((Mat[k][i-1] - 2*Mat[k][i] + Mat[k][i+1])/dx**2) + f)
        for i in range (int((0.25/dx) - 1/2),int((0.25/dx) + 1/2)):
            f = (10000 * (1-2*(k*dt)**2))/dx
            Mat[k+1][i] = Mat[k][i] + dt*(((Mat[k][i-1] - 2*Mat[k][i] + Mat[k][i+1])/dx**2) + f)
        for i in range (int(0.25/dx + 1/2),N):
            f = 0
            Mat[k+1][i] = Mat[k][i] + dt*(((Mat[k][i-1] - 2*Mat[k][i] + Mat[k][i+1])/dx**2) + f)
    return (Mat)
def ErroEvolucao(N, M, dt, dx, grafico):
    normal = grafico
    funcao = []
    erro = []
    k = M
    for k in range (M+1):
        erro.append([0]*(N+1))
    for k in range (M+1):
        funcao.append([0]*(N+1))
    for i in range (N):
        funcao[M][i] = (1 + np.sin(10*M*dt))*((i*dx)**2)*((1-i*dx)**2)
    for i in range (1,N):
        erro[k][i] = abs(normal[k][i] - funcao[k][i])
    return (erro[M])
def ErroEvolucaoB(N, M, dt, dx, graficasso):
    normal = graficasso
    funcao = []
    erro = []
    k = M
    for k in range (M+1):
        erro.append([0]*(N+1))
    for k in range (M+1):
        funcao.append([0]*(N+1))
    for i in range (N):
        funcao[M][i] = (np.exp(M*dt - i*dx))*np.cos(5*M*dt*i*dx)
    for i in range (1,N):
        erro[k][i] = abs(normal[k][i] - funcao[k][i])
    return (erro[M])
def MaximoErro(N, M, linha):
    maximo = 0
    for i in range(1,N):
        if (abs(linha[i])>maximo):
            maximo = linha[i]
    return maximo
def Decomposicao(N, M, dt, dx, erro, metodo, item):
    Mat = []
    for k in range (N-1):
        Mat.append([0]*(N-1))
    if metodo == 'crank':
        Mat[0][0] = 1 + erro
        for i in range (1,N-1):
            Mat[i-1][i] = -erro/2
            Mat[i][i] = 1 + erro
            Mat[i][i-1] = -erro/2
        L = []
        for k in range (N-1):
            L.append([0]*(N-1))
        D = []
        for k in range (N-1):
            D.append([0]*(N-1))
        D[0][0] = 1 + erro
    if metodo == 'euler':
        Mat[0][0] = 1 + 2*erro
        for i in range (1,N-1):
            Mat[i-1][i] = -erro
            Mat[i][i] = 1 + 2*erro
            Mat[i][i-1] = -erro
        L = []
        for k in range (N-1):
            L.append([0]*(N-1))
        D = []
        for k in range (N-1):
            D.append([0]*(N-1))
        D[0][0] = 1 + 2*erro
    for j in range (N-1):
        for i in range (j):
            h = Mat[j][i]
            L[j][i] = h/D[i][i]
            for k in range(i+1, j+1):
                Mat[j][k] = Mat[j][k] - h*L[k][i]
        D[j][j] = Mat[j][j]
    if item == 'a':
        U = []
        for k in range (M+1):
            U.append([0]*(N+1))
        for i in range (N+1):
            U[0][i] = ((i*dx)**2)*(1-i*dx)**2
        F = []
        for k in range (M+1):
            F.append([0]*(N+1))
        for k in range (M+1):
            for i in range(N+1):
                F[k][i] = (10*np.cos(10*k*dt)*(i*dx)**2)*(1-i*dx)**2 - (1 + np.sin(10*k*dt))*(12*(i*dx)**2 - 12*(i*dx) +2)
    if item == 'b':
        U = []
        for k in range (M+1):
            U.append([0]*(N+1))
        for i in range (N+1):
            U[0][i] = (np.exp(-i*dx))
        for k in range (M+1):
            U[k][0] = (np.exp(k*dt))
            U[k][N] = (np.exp(k*dt - 1))*np.cos(5*k*dt)
        F = []
        for k in range (M+1):
            F.append([0]*(N+1))
        for k in range (M+1):
            for i in range (N+1):
                t = k*dt
                x = i*dx
                F[k][i] = 25*(np.exp(t-x))*(t**2)*np.cos(5*t*x)-10*(np.exp(t-x))*t*np.sin(5*t*x)-5*((np.exp(t-x))*x*np.sin(5*t*x))
    if item == 'c':
        U = []
        for k in range (M+1):
            U.append([0]*(N+1))
        F = []
        for k in range (M+1):
            F.append([0]*(N+1))
        for k in range (M+1):
            for i in range (int(0.25/dx - 1/2)):
                F[k][i]= 0
            for i in range (int((0.25/dx) - 1/2),int((0.25/dx) + 1/2)):
                F[k][i] = (10000 * (1-2*(k*dt)**2))/dx
            for i in range (int(0.25/dx + 1/2),N):
                F[k][i] = 0
    K = []
    for i in range(M+1):
        K.append([0]*(N+1))
    for k in range (M):    
        B = []
        for q in range (N-1):
            B.append(0)
        if metodo == 'crank':
            for i in range(N-1):
                if (i==0):
                    B[i] = (erro/2)*U[k+1][i] + (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
                elif (i==N-2):
                    B[i] = (erro/2)*U[k+1][N] + (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
                else:
                    B[i] = (erro/2)*(U[k][i] - 2*U[k][i+1] + U[k][i+2]) + (dt/2)*(F[k][i+1] + F[k+1][i+1]) + U[k][i+1]
        if metodo == 'euler':
            for i in range(N-1):
                if (i==0):
                    B[i] = U[k][i+1] + dt*(F[k+1][i+1]) + erro*U[k+1][0]
                elif (i==N-2):
                    B[i] = U[k][i+1] + dt*(F[k+1][i+1]) + erro*U[k+1][N]
                else:
                    B[i] = U[k][i+1] + dt*(F[k+1][i+1])
        Z = []
        for v in range (N-1):
            Z.append(0)
        C = []
        for a in range (N-1):
            C.append(0)
        for j in range (N-1):
            Z[j] = B[j]
            for i in range (j):
                Z[j] = Z[j] - L[j][i]*Z[i] 
            C[j] = Z[j]/D[j][j]
        X = []
        for o in range (N-1):
            X.append(0)  
        for t in reversed(range(N-1)):
            X[t] = C[t]
            for i in range (t+1, N-1):
                X[t] = X[t] - L[i][t]*X[i]
        for i in range (1,N):
            U[k+1][i] = X[i-1]
    return U
def Plot3D(N, M, dt, dx, grafico):
    x = []
    t = []
    for i in range(N+1):
        x.append(i*dx)
    for k in range(M+1):
        t.append(k*dt)
    matriz = grafico      
    fig = plt.figure()
    ax = Axes3D(fig)
    x, t = np.meshgrid(x, t)
    u = np.array(matriz)
    surf = ax.plot_surface(x, t, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_xlabel("Posição")
    ax.set_ylabel("Tempo")
    ax.set_zlabel("Temperatura")
def main():
    N = int(input ('N:'))
    parte = (input ('Parte (1 ou 2):'))
    if parte == '1':
        erro = float(input('lambda:'))
        M = int((N*N)/erro)
        T = 1
        dt = T/M
        dx = 1/N
    if parte == '2':
        erro = N
        M = N
        dt = 1/N
        dx = dt
    if parte == '1':
        item = (input ('Item (a, b ou c):'))
        if item == 'a':
            grafico = EvolucaoItemA(N,M,dt,dx)
            Plot3D(N, M, dt, dx, grafico)
            linha = ErroEvolucao(N, M, dt, dx,grafico)
            print('O erro para este N é ')
            print(MaximoErro(N, M, linha))
        if item == 'b':
            grafico = EvolucaoItemB(N,M,dt,dx)
            Plot3D(N,M,dt,dx,grafico)
            linha = ErroEvolucaoB(N, M, dt, dx, grafico)
            print('O erro para este N é ')
            print(MaximoErro(N, M, linha))
        if item == 'c':
            grafico = EvolucaoItemC(N,M,dt,dx)
            Plot3D(N,M,dt,dx,grafico)
    if parte == '2':
        metodo = (input('Metodo (crank ou euler):'))
        item = (input('Item (a, b ou c):'))
        if item == 'a':
            graficasso = Decomposicao(N, M, dt, dx, erro, metodo, item)
            Plot3D(N, M, dt, dx, graficasso)
            linha = ErroEvolucao(N,M,dt,dx,graficasso)
            print('O erro para este N é ')
            print(MaximoErro(N, M, linha))
        if item == 'b':
            graficasso = Decomposicao(N, M, dt, dx, erro, metodo, item)
            Plot3D(N, M, dt, dx, graficasso)
            linha = ErroEvolucaoB(N,M,dt,dx,graficasso)
            print('O erro para este N é ')
            print(MaximoErro(N, M, linha))
        if item == 'c':
            graficasso = Decomposicao(N, M, dt, dx, erro, metodo, item)
            Plot3D(N, M, dt, dx, graficasso)
main()