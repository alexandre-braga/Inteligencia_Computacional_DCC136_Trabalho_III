# Trabalho de Grafos
Repositório para armazenar o código fonte do trabalho 3 da disciplina
DCC136 - Inteligência Computacional 2023.3 (UFJF).

Integrantes:
- Alexandre Vitor Silva Braga - 201965501B
- Thais de Jesus Soares       - 202065511B

## Instruções para compilação
Digite `make` no terminal e o bínario `aco_ta` será gerado.

`
Uso: ./aco_ta ARQUIVO_ENTRADA ARQUIVO_SAIDA
`

## Parâmetros
#### Número de iterações
1. Guloso Randomizado: 500
1. Guloso Randomizado Reativo: 5000
1. Colônia de Formigas: 10

#### Número de Formigas
1. Colônia de Formigas: 500

#### Alpha
1. Guloso Randomizado: 0.05
1. Guloso Randomizado Reativo: { 0.05, 0.10, 0.15, 0.30, 0.5 }

#### Número de iterações
Cada algoritmo está sendo executado 10 vezes para achar a média de tempo e do
número de rótulos de cada iteração.

#### Semente de randomização
A semente de randomização será o número de segundos desde 1/1/1970 (Unix Time)
retornado pela chamada de função `time(NULL)`.
