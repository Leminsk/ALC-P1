***********************************
    Modo de uso de "task3.exe"    
***********************************

task3.exe tenta estimar o valor respectivo de uma função desconhecida dado um valor de entrada (x) e um conjunto de pares de coordenadas (x,j).

Para o funcionamento correto de task2.exe, é necessário fornecer alguns arquivos de entrada que precisam necessariamente estar localizados na mesma pasta/diretório que o executável:

    - "main_input.txt" : cada linha corresponde a um parâmetro de uso (todas as linhas devem ser preenchidas independente do método)
        Linha 1: Quantia de pares de coordenadas fornecidas
        Linha 2: código para seleção do método a ser utilizado (ICOD)
            1: Interpolação de Lagrange
            2: Regressão Multilinear (NÃO IMPLEMENTADO)
        Linha 3: o nome do arquivo txt contendo a matriz A (veja abaixo para mais detalhes)
        Linha 4: a coordenada/valor (x) a ser utilizado para a estimativa

    - "coordenadas.txt" : arquivo de dados do conjunto de coordenadas
        Cada linha corresponde a um par de coordenada (x,y), e cada elemento deve ser separado por vírgula do próximo.
        E.g.: "0.10005, 73.99" (a quantia de espaços entre vírgulas e elementos não interfere no programa)


Ao concluir seus cálculos/processos, task3.exe gerará um arquivo de saída chamado "main_output.txt" contendo:
    - o valor estimado da função

Caso ocorra algum erro, o arquivo de saída conterá informações sobre o motivo da falha do programa.