#ifndef GRAFO_HPP
#define GRAFO_HPP

#include <map>
#include <vector>
#include <ostream>
#include <istream>
#include <memory>

#include "Vertice.hpp"

class Grafo {
    private:
        /**
         * Estrutura auxiliar que será usada para relacionar
         * um rótulo r com as arestas que possuem esse rótulo.
         */
        struct ArestaAux {
            idvertice_t idOrigem;
            idvertice_t idDestino;

            ArestaAux() {};
            ArestaAux(idvertice_t origem, idvertice_t destino):
                idOrigem(origem),
                idDestino(destino)
            {}
        };

        using labelcontainer_t = std::map<rotulo_t, std::list<ArestaAux>>;
        using vertexlist_t = std::vector<Vertice>;

        // Atributos
        vertexlist_t listaVertices;
        labelcontainer_t rotulos;

    public:
        Grafo(size_t numeroDeVertices, size_t numeroDeRotulos);
        static Grafo gerarDoArquivo(std::istream& arqEntrada);
        struct lcp
        {
            rotulo_t cand;
            double prob;
        };

        size_t numeroDeVertices() const;
        size_t numeroDeRotulos() const;

        void fazerAresta(idvertice_t id1, idvertice_t id2, rotulo_t rotulo);
        void desfazerAresta(idvertice_t id1, idvertice_t id2, rotulo_t rotulo);

        void toDots(std::ostream& arqSaida) const;

        Grafo algoritmoGuloso() const;
        Grafo algoritmoGulosoRandomizado(double alfa, size_t nIteracoes) const;
        Grafo algoritmoGulosoRandomizadoReativo(
                const std::vector<double>& alfas, size_t nIteracoes,
                size_t tamanhoBloco, double *alfaMaiorProb = nullptr) const;
        Grafo algoritmoACO(size_t nIteracoes, size_t nFormigas, double tau_min, double tau_max) const;
        Grafo algoritmoACOSmoothing(size_t nIteracoes, size_t nFormigas, size_t bloco, double lambda_max, double zeta_max, double tau_min, double tau_max) const;
        Grafo algoritmoACOSmoothingV2(size_t nIteracoes, size_t nFormigas, size_t bloco, double lambda_max, double zeta_max, double tau_min, double tau_max) const;
        Grafo algoritmoACOBuscaLocal(size_t nIteracoes, size_t nFormigas, double tau_min, double tau_max) const;

    private:
        Grafo algoritmoGulosoHelper(double alfa) const;
        Grafo algoritmoACOHelper(std::unique_ptr<double[]>& vetorProb, std::unique_ptr<double[]>& vetorFeromonio, std::unique_ptr<double[]>& vetorHeuristico, size_t m, size_t minRotulos) const;

        Grafo removeKRotulos(Grafo sol, size_t k, std::unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t>& listCand) const;
        Grafo reparaSolucao(Grafo sol, std::unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t>& listCand, size_t k) const;
        Grafo VNS(Grafo solAtual, std::unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t> listCand) const;
        Grafo algoritmoACOHelperBL(std::unique_ptr<double[]>& vetorFeromonio, std::unique_ptr<double[]>& vetorHeuristico, size_t m, size_t minRotulos, std::unique_ptr<idvertice_t[]>& subArvore,  std::vector<lcp>& LCP) const;

        const Vertice *getVerticeById(idvertice_t id) const;
        Vertice *getVerticeById(idvertice_t id);
};

#endif // GRAFO_HPP
