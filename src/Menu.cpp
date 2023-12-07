#include "Menu.hpp"

#include <limits>
#include <chrono>
#include <iostream>

using namespace std;

using relogio_t = chrono::high_resolution_clock;
using tempo_t   = relogio_t::time_point;

#define DELTA_T(fim, inicio, duration) \
    chrono::duration_cast<duration>((fim) - (inicio)).count()

#define BOILER_PLATE(param, duration, durationFmt, iter, additionalCode)       \
    do {                                                                       \
        double nRotulosTotal = 0;                                              \
        double tempoTotal = 0;                                                 \
                                                                               \
        for (size_t i = 0; i < iter; ++i) {                                    \
            tempo_t inicio = relogio_t::now();                                 \
            const Grafo resultado = param;                                     \
            tempo_t fim = relogio_t::now();                                    \
            double duracao = DELTA_T(fim, inicio, duration);                   \
                                                                               \
            resultado.toDots(out);                                             \
                                                                               \
            std::cerr << "Número de rótulos mínimos encontrados (Execução "    \
                      << i + 1 << "): " << resultado.numeroDeRotulos()         \
                      << "\n";                                                 \
            std::cerr << "Duração: " << duracao << " [" durationFmt "]\n";     \
            nRotulosTotal += resultado.numeroDeRotulos();                      \
            tempoTotal += duracao;                                             \
            additionalCode;                                                    \
        }                                                                      \
        std::cout << "Média do número de rótulos encontrados: "                \
                  << nRotulosTotal / iter << "\n";                             \
        std::cout << "Média da duração de cada iteração: "                     \
                  << tempoTotal / iter << " [" durationFmt "]\n";              \
    } while (0)

static void guloso(const Grafo& g, std::ofstream& out)
{
    BOILER_PLATE(g.algoritmoGuloso(), chrono::microseconds, "µs", 10, {});
}

static void gulosoRandomizado(const Grafo& g, std::ofstream& out)
{
    BOILER_PLATE(g.algoritmoGulosoRandomizado(0.05, 500), chrono::microseconds, "µs", 10, {});
}

static void gulosoRandomizadoReativo(const Grafo& g, std::ofstream& out)
{
    double alfaMaiorProb;

    BOILER_PLATE(g.algoritmoGulosoRandomizadoReativo(
                 { 0.05, 0.10, 0.15, 0.30, 0.5 }, 5000, 10, &alfaMaiorProb),
                 chrono::milliseconds, "ms", 10,
                 std::cout << "Alpha: " << alfaMaiorProb << '\n');
}

static void coloniaFormigas(const Grafo& g, std::ofstream& out)
{
    BOILER_PLATE(g.algoritmoACO(600, 1000, 10, 0.2, 0.001, 20), chrono::milliseconds, "ms", 10, {});
}

static void toDots(const Grafo& g, std::ofstream& out)
{
    g.toDots(out);
}

static struct MenuOption {
    const char *optName;
    void (*action)(const Grafo& g, std::ofstream&);
} menuOpts[] = {
    { "Algoritmo guloso"                    , guloso                   },
    { "Algoritmo guloso randomizado"        , gulosoRandomizado        },
    { "Algoritmo guloso randomizado reativo", gulosoRandomizadoReativo },
    { "Algoritmo Colônia de Formigas"       , coloniaFormigas          },
    { "To Dots"                             , toDots                   }
};

#define OPTS_SIZE (sizeof(menuOpts) / sizeof(MenuOption))

#define IS_VALID(opt) (opt <= OPTS_SIZE)

#define SAIR 0

void Menu::run(const Grafo& g, std::ofstream &out)
{
    std::string tmpBuf;
    size_t selectedOpt;

    while (true) {
        std::cout << "*****************************************\n";
        std::cout << "** Árvore Geradora De Rotulação Mínima **\n";
        std::cout << "*****************************************\n";
        std::cout << "Número de vértices: " << g.numeroDeVertices() << '\n';
        std::cout << "Número de rotulos: " << g.numeroDeRotulos() << '\n';

        for (size_t i = 0; i < OPTS_SIZE; ++i) {
            std::cout << (i + 1) << " - " << menuOpts[i].optName << '\n';
        }
        std::cout << SAIR << " - Sair\n";

        std::cout << "\nR: ";

        std::getline(std::cin, tmpBuf);

        try {
            selectedOpt = std::stoul(tmpBuf);
        } catch (const std::invalid_argument& e) {
            std::cerr << "** Entrada inválida **\n\n";
            continue;
        }

        if (IS_VALID(selectedOpt)) {
            if (selectedOpt == SAIR) {
                return;
            }
            (*menuOpts[selectedOpt - 1].action)(g, out);
            out.flush();
        } else {
            std::cerr << " ** OPÇÃO INVÁLIDA **\n\n";
        }
    }
}
