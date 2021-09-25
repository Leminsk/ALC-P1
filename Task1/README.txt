#=================================#
|   Modo de uso de "task1.exe"    |
#=================================#

task1.exe tenta resolver o sistema AX = B, no qual X é desconhecido.
A é matriz quadrada e B é um vetor de mesma ordem que A.

Para o funcionamento correto de task1.exe, é necessário fornecer alguns arquivos de entrada que precisam necessariamente estar localizados na mesma pasta/diretório que o executável:

    - "main_input.txt" : cada linha corresponde a um parâmetro de uso (todas as linhas devem ser preenchidas independente do método)
        Linha 1: Ordem da matriz A (N)
        Linha 2: código para seleção do método a ser utilizado (ICOD)
            1: Decomposição LU
            2: Decomposição Cholesky
            3: Jacobi Iterativo
            4: Gauss-Seidel Iterativo
        Linha 3: flag de cálculo do determinante da matriz A
            0: não calcular
            maior que 0: calcular
        Linha 4: o nome do arquivo txt contendo a matriz A (veja abaixo para mais detalhes)
        Linha 5: o nome do arquivo txt contendo o vetor B (veja abaixo para mais detalhes)
        Linha 6: a tolerância a ser utilizada caso o método escolhido seja iterativo

    - "matrix_A.txt" : arquivo de dados da matriz A
        Cada linha corresponde a uma linha (row) da matriz, e cada elemento deve ser separado por vírgula do próximo.
        E.g.: "4, 5.1, 3.99" (a quantia de espaços entre vírgulas e elementos não interfere no programa)

    - "matrix_B.txt" : arquivo de dados do vetor B
        Cada linha corresponde a um elemento de B. (i.e. cada elemento é separado por NEWLINE)

Ao concluir seus cálculos/processos, task1.exe gerará um arquivo de saída chamado "main_output.txt" contendo:
    - solução do sistema X
    - determinante (se requisitado)
    - número de iterações
    - histórico do erro (TOL)

Caso ocorra algum erro, o arquivo de saída conterá informações sobre o motivo da falha do programa.
