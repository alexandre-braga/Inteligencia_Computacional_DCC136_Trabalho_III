/************************************************************
*   TRABALHO DE GRAFOS - Árvore Geradora De Rotulação Mínima
*       Alexandre Vitor Silva Braga - 201965501B
*       Thais de Jesus Soares - 202065511B
************************************************************/

#include <iostream>

#include "Grafo.hpp"
#include "Menu.hpp"

#define EXPECTED_ARGC(argc) ((argc) == 3)

#define ERR_INVALID_ARGC        1
#define ERR_COULD_NOT_OPEN_FILE 2
#define ERR_READ_FAIL           3

#define OPEN_FILE(fileObj, fileName)                                           \
    do {                                                                       \
        fileObj.open(fileName);                                                \
                                                                               \
        if (!fileObj.is_open()) {                                              \
            std::cerr << programName                                           \
                << ": não foi possível abrir o arquivo `"                      \
                << fileName << "`\n";                                          \
            return ERR_COULD_NOT_OPEN_FILE;                                    \
        }                                                                      \
    } while (0)

int main(int argc, const char *const argv[])
{
    srand(time(NULL));
    const char *programName = *argv;

    if (!EXPECTED_ARGC(argc)) {
        std::cerr << "Uso: " << programName
            << " ARQUIVO_ENTRADA ARQUIVO_SAIDA\n";
        return ERR_INVALID_ARGC;
    }

    std::ifstream arqEntrada;
    std::ofstream arqSaida;

    OPEN_FILE(arqEntrada, argv[1]);
    OPEN_FILE(arqSaida,   argv[2]);

    Grafo g = Grafo::gerarDoArquivo(arqEntrada);

    std::cin.exceptions(std::istream::failbit | std::istream::badbit);

    try {
        Menu::run(g, arqSaida);
    } catch (const std::istream::failure& e) {
        std::cerr << "Falha de leitura\n";
        std::cout << "Saindo...\n";
        return ERR_READ_FAIL;
    }
    return 0;
}
