# ALC-P1
Primeiro Trabalho Disciplinar para COC473 - Álgebra Linear Computacional consiste em criar uma pequena biblioteca de cálculo numérico para matrizes, incluindo métodos para:  
- Task 1 - Solução de um sistema linear de equações AX = B utilizando:  
  - Decomposição LU
  - Decomposição de Cholesky
  - Procedimento Iterativo Jacobi
  - Procedimento Iterativo Gauss-Seidel  
- Task 2 - Calcular autovalores e autoverores de uma matriz utilizando:  
  - Método da Potência
  - Método de Jacobi  
- Task 3 - Fit de pontos de uma função e estimativa de um valor y desconhecido para um dado x por:  
  - Interpolação (Lagrange)
  - Regressão Multilinear (NÃO implementado)

## Compilação
Algumas funções utilizadas nos programas necessitam do C++11, portanto é preciso compilar com a flag adequada:  
- `g++ -std=c++11 task1.cpp -o task1`  
- `g++ -std=c++11 task2.cpp -o task2`  
- `g++ -std=c++11 task3.cpp -o task3`  

Os executáveis (.exe) presentes neste repositório foram compilados em uma máquina com Windows 10 e g++ (tdm64-1) 5.1.0.

## Modo de Uso
TODO: preencher (tem no um README.txt pra cada um por agora)
