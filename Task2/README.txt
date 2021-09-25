+---------------------------------+
|   Modo de uso de "task2.exe"    |
+---------------------------------+

task2.exe tenta calcular os autovalores e autovetores de uma matriz A.

Para o funcionamento correto de task2.exe, é necessário fornecer alguns arquivos de entrada que precisam necessariamente estar localizados na mesma pasta/diretório que o executável:

    - "main_input.txt" : cada linha corresponde a um parâmetro de uso (todas as linhas devem ser preenchidas independente do método)
        Linha 1: Ordem da matriz A (N)
        Linha 2: código para seleção do método a ser utilizado (ICOD)
            1: Método da Potência
            2: Método Jacobi
        Linha 3: flag de cálculo do determinante da matriz A
            0: não calcular
            maior que 0: calcular
        Linha 4: o nome do arquivo txt contendo a matriz A (veja abaixo para mais detalhes)
        Linha 5: a tolerância a ser utilizada caso o método escolhido seja iterativo

    - "matrix_A.txt" : arquivo de dados da matriz A
        Cada linha corresponde a uma linha (row) da matriz, e cada elemento deve ser separado por vírgula do próximo.
        E.g.: "4, 5.1, 3.99" (a quantia de espaços entre vírgulas e elementos não interfere no programa)


Ao concluir seus cálculos/processos, task2.exe gerará um arquivo de saída chamado "main_output.txt" contendo:
    - autovalores de A
    - autovetores respectivos
    - determinante de A (se requisitado)
    - número de iterações
    - histórico do erro (TOL)

Caso ocorra algum erro, o arquivo de saída conterá informações sobre o motivo da falha do programa.