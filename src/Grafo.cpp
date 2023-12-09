#include "Grafo.hpp"
#include "Aresta.hpp"
#include "utils.hpp"

#include <cmath>
#include <memory>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <iostream>

using namespace std;

Grafo::Grafo(size_t numeroDeVertices, size_t numeroDeRotulos __attribute__((unused)))
{
    // Inicializa a lista de vértice com seus respectivos IDs.
    this->listaVertices.reserve(numeroDeVertices);
    for (idvertice_t id = 0; id < numeroDeVertices; ++id) {
        this->listaVertices.push_back(Vertice(id));
    }
}

Grafo Grafo::gerarDoArquivo(std::istream& arqEntrada)
{
#define INVALID_LABEL(label) ((label) >= nRotulos)

    size_t nVertices, nRotulos;
    arqEntrada >> nVertices >> nRotulos;

    Grafo g(nVertices, nRotulos);

    std::string linha;
    rotulo_t rotulo;
    size_t pos, j;

    arqEntrada.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    for (idvertice_t i = 0; i < nVertices - 1; i++) {
        std::string tmp;

        std::getline(arqEntrada, linha);
        pos = 0;
        j = i + 1;

        while (stringTok(linha, tmp, ' ', &pos)) {
            rotulo = stringToRotulo(tmp);

            if (!INVALID_LABEL(rotulo)) {
                g.fazerAresta(i, j, rotulo);
            }
            ++j;
        }
    }
    return g;
}

Vertice *Grafo::getVerticeById(idvertice_t id)
{
    return &this->listaVertices[id];
}

const Vertice *Grafo::getVerticeById(idvertice_t id) const
{
    return &this->listaVertices[id];
}

size_t Grafo::numeroDeVertices() const
{
    return this->listaVertices.size();
}

size_t Grafo::numeroDeRotulos() const
{
    return this->rotulos.size();
}

void Grafo::fazerAresta(idvertice_t id1, idvertice_t id2,
                        rotulo_t rotulo)
{
    /* Testa se um vértice existe. Caso contrário, lança uma excessão
     * do tipo std::invalid_argument
     */
#define TEST_VERTICE(v, id)                                                    \
    do {                                                                       \
        if ((v) == nullptr) {                                                  \
            throw std::invalid_argument(                                       \
                    std::string("não foi possível encontrar grafo com id ") += \
                    (id));                                                     \
        }                                                                      \
    } while (0)

    if (id1 != id2) {
        Vertice *v1 = this->getVerticeById(id1);
        Vertice *v2 = this->getVerticeById(id2);

        TEST_VERTICE(v1, id1);
        TEST_VERTICE(v2, id2);

        v1->addAdjacente(id2, rotulo);
        v2->addAdjacente(id1, rotulo);
    } else {
        Vertice *v = this->getVerticeById(id1);

        TEST_VERTICE(v, id1);
        v->addAdjacente(id1, rotulo);
    }

    // Tipo de retorno do método std::map::insert()
    // Um pair onde:
    //    - first: iterador que aponda para um par <chave-valor> dentro do std::map
    //    - second: bool: `true` se a chave foi inserida. `false` se a chave já existia dentro do container
    std::pair<Grafo::labelcontainer_t::iterator, bool> ret;

    // Se não existe uma chave `rotulo` em `this->rotulos`, inserir uma lista vazia
    ret = this->rotulos.insert(make_pair(rotulo, std::list<ArestaAux>{}));

    // Adiciona a aresta (id1, id2) à lista de arestas associadas ao rótulo
    ret.first->second.push_back(ArestaAux(id1, id2));
}

void Grafo::toDots(std::ostream& arqSaida) const
{
    arqSaida << "graph {\n";
    for (const Vertice& v : this->listaVertices) {
        for (const Aresta& a : v.listaDeAdjacencia()) {
            if (v.id() < a.idDestino()) {
                arqSaida << "\t" << v.id() << " -- " << a.idDestino()
                    << "[label=" << a.rotulo() << "];\n";
            }
        }
    }
    arqSaida << "}";
}

// Usado pelo algoritmo de kruskal abaixo para testar se dois vértices estão na mesma subárvore
static inline bool estaoNaMesmaSubArvore(const idvertice_t subArvore[], idvertice_t u, idvertice_t v)
{
    return subArvore[u] == subArvore[v];
}

static void juntarSubArvores(idvertice_t subArvore[], idvertice_t u, idvertice_t v, size_t size)
{
    idvertice_t novaSub, subAntiga;

    if (subArvore[u] < subArvore[v]) {
        novaSub = subArvore[u], subAntiga = subArvore[v];
    } else {
        novaSub = subArvore[v], subAntiga = subArvore[u];
    }

    for (idvertice_t id = 0; id < size; ++id) {
        if (subArvore[id] == subAntiga) {
            subArvore[id] = novaSub;
        }
    }
}

/**
 * Algoritmo baseado no algoritmo de kruskal.
 * O algoritmo funciona da seguinte forma:
 * Enquanto (!conexo(grafoSolucao) && !vazia(listaDeCandidatos)) {
 *     r = escolheRotulo(listaDeCandidatos);
 *     Para cada aresta (u, v) que tem rótulo r {
 *         Se (!estaoNaMesmaSubArvore(u, v)) {
 *             adicionar a aresta (u, v) ao grafo solução.
 *         }
 *     }
 * }
 */
Grafo Grafo::algoritmoGulosoHelper(double alfa) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())

    Grafo F(this->numeroDeVertices(), 0);
    std::vector<rotulo_t> L;
    unique_ptr<idvertice_t[]> subArvore(new idvertice_t[this->numeroDeVertices()]);

    L.reserve(this->numeroDeRotulos());

    for (rotulo_t r = 0; r < this->numeroDeRotulos(); ++r) {
        L.push_back(r);
    }

    std::sort(L.begin(), L.end(), [this](rotulo_t r1, rotulo_t r2) {
        return PESO(r1) < PESO(r2);
    });

    for (const Vertice& v : this->listaVertices) {
        subArvore[v.id()] = v.id();
    }

    size_t subs = this->numeroDeVertices();
    while (subs > 1 && !L.empty()) {

        //lista restrita
        size_t aux = L.size() * alfa;
        //decide o rótulo
        size_t idx = (aux < 2) ? 0 : (rand() % aux);
        rotulo_t r = L[idx];

        //adiciona todas arestas do rótulo
        for (const ArestaAux& a : this->rotulos.find(r)->second) {
            //verifica se aresta do rótulo Y já não foi colocada em rótulo X antes
            if (!estaoNaMesmaSubArvore(subArvore.get(), a.idOrigem, a.idDestino)) {
                F.fazerAresta(a.idOrigem, a.idDestino, r);
                juntarSubArvores(subArvore.get(), a.idOrigem, a.idDestino, this->numeroDeVertices());

                if (--subs <= 1) {
                    return F;
                }
            }
        }
        L.erase(L.begin() + idx);
    }
    return F;
}

Grafo Grafo::algoritmoGuloso() const
{
    return this->algoritmoGulosoHelper(0);
}

Grafo Grafo::algoritmoGulosoRandomizado(double alfa, size_t nIteracoes) const
{
    Grafo melhorSol(this->numeroDeVertices(), this->numeroDeRotulos());
    Grafo solAux(this->numeroDeVertices(), 0);

    bool isPrimeiraVez = true;
    while (nIteracoes) {
        solAux = this->algoritmoGulosoHelper(alfa);
        if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
            isPrimeiraVez = false;
            melhorSol = std::move(solAux);
        }
        --nIteracoes;
    }
    return melhorSol;
}

/* Estrutura auxiliar usada no algoritmoGulosoRandomizadoReativo para calcular
 * a média de um alfa
 */
struct Media {
    size_t total;
    size_t nVezes;

    Media():
        total(0),
        nVezes(0)
    {}

    // Calcula a média
    inline operator double() const
    {
        return (this->nVezes == 0) ? 0 : (double) this->total / this->nVezes;
    }
};

static size_t selecionaAlfa(unique_ptr<double[]>& P, size_t m)
{
    double rng = 0.0001 + (double) (rand() % 9900) / 100000;
    double aux = 0;

    for (size_t i = 0; i < m; i++) {
        aux += P[i];
        if (rng <= aux) {
            return i;
        }
    }
    return m - 1;
}

static inline void atualizarMedias(unique_ptr<Media[]>& A, size_t numeroDeRotulos, size_t i)
{
    A[i].total += numeroDeRotulos;
    ++A[i].nVezes;
}

static inline void atualizarProbabilidades(unique_ptr<double[]>& P,
        unique_ptr<Media[]>& A, size_t m, size_t melhorSol)
{
#define SIGMA 10

    unique_ptr<double[]> q(new double[m]);
    double somatorioDeQjs = 0;

    for (size_t i = 0; i < m; ++i) {
        q[i] = std::pow(melhorSol / A[i], SIGMA);
        somatorioDeQjs += q[i];
    }

    for (size_t i = 0; i < m; ++i) {
        P[i] = q[i] / somatorioDeQjs;
    }
}

Grafo Grafo::algoritmoGulosoRandomizadoReativo(const vector<double>& alfas,
        size_t nIteracoes, size_t bloco, double *alfaMaiorProb) const
{
    Grafo melhorSol(this->numeroDeVertices(), this->numeroDeRotulos());
    Grafo solAux(this->numeroDeVertices(), 0);

    unique_ptr<double[]> P(new double[alfas.size()]);
    unique_ptr<Media[]> A(new Media[alfas.size()]);

    bool isPrimeiraVez = true;
    size_t i;

    do {
        double pInicial = (double) 1 / alfas.size();

        for (size_t i = 0; i < alfas.size(); ++i) {
            P[i] = pInicial;
        }
    } while (0);

    for (i = 1; i <= alfas.size() && i <= nIteracoes; ++i) {
        size_t alfaidx = i - 1;

        solAux = this->algoritmoGulosoHelper(alfas[alfaidx]);

        atualizarMedias(A, solAux.numeroDeRotulos(), alfaidx);


        if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
            isPrimeiraVez = false;
            melhorSol = std::move(solAux);
        }
    }

    atualizarProbabilidades(P, A, alfas.size(), melhorSol.numeroDeRotulos());

    for (; i <= nIteracoes; ++i) {
        if (i % bloco == 0) {
            atualizarProbabilidades(P, A, alfas.size(), melhorSol.numeroDeRotulos());
        }

        size_t alfaidx = selecionaAlfa(P, alfas.size());

        solAux = this->algoritmoGulosoHelper(alfas[alfaidx]);

        atualizarMedias(A, solAux.numeroDeRotulos(), alfaidx);


        if (solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
            isPrimeiraVez = false;
            melhorSol = std::move(solAux);
        }
    }

    if (alfaMaiorProb) {
        double maiorProb = P[0];
        *alfaMaiorProb = alfas[0];

        for (size_t i = 1; i < alfas.size(); ++i) {
            if (P[i] > maiorProb) {
                maiorProb = P[i], *alfaMaiorProb = alfas[i];
            }
        }
    }
    return melhorSol;
}

/*-------------------------------------------------------------*/
/*--------------------------- ACO -----------------------------*/
/*-------------------------------------------------------------*/
struct lcp
{
    rotulo_t cand;
    double prob;
};
struct TeN
{
    rotulo_t id;
    double T;
    double N;
    bool content;
};
static inline void inicializarFeromonios(unique_ptr<double[]>& vetorFeromonio, size_t m, double tau_min, double tau_max)
{
    double tau_media = (tau_max + tau_min)/2;
    for (size_t i = 0; i < m; i++) {
        vetorFeromonio[i] = tau_media;
    }
}
static inline void inicializarProbabilidadesACO(unique_ptr<double[]>& vetorProb, unique_ptr<double[]>& vetorHeuristico, size_t tot, size_t m)
{
    do {
        for (size_t i = 0; i < m; ++i) {
            vetorProb[i] = vetorHeuristico[i] / tot;
        }
    } while (0);
}

static inline void atualizarFeromonios(unique_ptr<double[]>& vetorFeromonio, size_t m, std::vector<rotulo_t> rotulosDaMelhorSol)
{
#define RHO 0.1
#define DELTA_L 0.5

    //deposition (Só nas arestas da melhor solução)
    for (size_t i = 0; i < m; ++i) {
        //if i ta no vector de rótulos da melhor solução (criar vector com os rotulos da melhorSol)
        if(std::find(rotulosDaMelhorSol.begin(), rotulosDaMelhorSol.end(), i) != rotulosDaMelhorSol.end()){
            vetorFeromonio[i] = vetorFeromonio[i] + DELTA_L;
        }
    }

    //evaporation
    for (size_t i = 0; i < m; ++i) {
        vetorFeromonio[i] = vetorFeromonio[i] * (1 - RHO);
    }
}
static inline void atualizarProbabilidadesACO(std::vector<lcp>& LCP, unique_ptr<double[]>& vetorFeromonio,
unique_ptr<double[]>& vetorHeuristico, size_t m)
{
#define PESO(r) (-this->rotulos.find(r)->second.size())

    std::vector<TeN> ten(m);
    ten.reserve(m);
    double somatorioDeTNs = 0;

    for (rotulo_t i = 0; i < m; i++ ){
        ten.at(i).id = i;
        ten.at(i).T = 0;
        ten.at(i).N = 0;
        ten.at(i).content = false;
    }

    for (const auto& element : LCP) {
        rotulo_t r = element.cand;
        ten.at(r).id = r;
        ten.at(r).T = vetorFeromonio[r];
        ten.at(r).N = vetorHeuristico[r];
        ten.at(r).content = true;
        somatorioDeTNs += ten.at(r).T * ten.at(r).N;
    }

    // for (const auto& element : ten) {
    //     std::cerr << "ten at r " << element.id << "|" << element.T << "|" << element.N << "|" << element.content << std::endl;
    // }
    // std::cerr << "lcp size:" << LCP.size() << std::endl;
    // std::cerr << "ten size:" << ten.size() << std::endl;

    for (rotulo_t r = 0; r < LCP.size(); r++){

        auto it = std::find_if(ten.begin(), ten.end(), [&LCP, r](const TeN& element) {
            return element.id == LCP.at(r).cand;
        });
        size_t itPos = std::distance(ten.begin(), it);

        // std::cerr << "itPos in Ten: " << itPos << "| Rótulo candidato a atualiza" << LCP.at(r).cand << std::endl;
        // std::cerr << "Verifying Ten at r " << ten.at(itPos).id << "|" << ten.at(itPos).T << "|" << ten.at(itPos).N << "|" << ten.at(itPos).content << std::endl;
        // std::cerr << "LCP at r before " << LCP.at(r).cand << "|" << LCP.at(r).prob << std::endl;

        if(ten.at(itPos).content){
            LCP.at(r).prob = ten.at(r).T * ten.at(r).N / somatorioDeTNs;
        }
        else{
            LCP.at(r).prob = -1;
        }
        // std::cerr << "LCP at r after " << LCP.at(r).cand << "|" << LCP.at(r).prob << std::endl;
    }
}
// static inline void atualizarProbabilidadesACO(std::vector<lcp> LCP, unique_ptr<double[]>& vetorFeromonio,
// unique_ptr<double[]>& vetorHeuristico, size_t m)
// {
// #define PESO(r) (-this->rotulos.find(r)->second.size())

//     std::vector<TeN[]> ten(m);
//     ten.reserve(m);
//     double somatorioDeTNs = 0;

//     for (rotulo_t i = 0; i < m; ++i) {
//         ten.at(i).T = vetorFeromonio[i];
//         ten.at(i).N = vetorHeuristico[i];
//         auto it = std::find_if(LCP.begin(), LCP.end(), [i](const lcp& element) { return element.cand == i; });
//         if(it != LCP.end()){
//             somatorioDeTNs += T[i]*N[i];
//         }
//     }

//     for (size_t i = 0; i < m; ++i) {
//         auto it = std::find_if(LCP.begin(), LCP.end(), [i](const lcp& element) { return element.cand == i; });

//         if(it != LCP.end()){
//             LCP.at(i).prob = T[i]*N[i] / somatorioDeTNs;
//         }
//         else{
//             LCP.at(i).prob = -1;
//         }
//     }
// }
static size_t selecionaRotulo(std::vector<lcp> LCP, size_t tamAtual)
{
    double rng = (double) (rand() % 10000) / 10000;
    double aux = 0;
    LCP.resize(tamAtual);
    std::sort(LCP.begin(),LCP.end(), [](const lcp &x, const lcp &y){ return (x.prob > y.prob);});

    // for (size_t r = 0; r < tamAtual; ++r) {
    //     std::cerr << "Sorted: " << LCP.at(r).cand << "|" << LCP.at(r).prob;
    // }

    for (size_t i = 0; i < tamAtual; i++) {
        aux += LCP.at(i).prob;
        if (rng <= aux) {
            // std::cerr << "\nRetornou Cand: " << LCP.at(i).cand << "| Prob: " << LCP.at(i).prob << std::endl;
            return LCP.at(i).cand;
        }
    }
    return LCP.at(tamAtual-1).cand;
}


Grafo Grafo::algoritmoACOHelper(unique_ptr<double[]>& vetorProb, unique_ptr<double[]>& vetorFeromonio,
unique_ptr<double[]>& vetorHeuristico, size_t m, size_t minRotulos) const
{

    Grafo F(this->numeroDeVertices(), 0);
    std::vector<lcp> LCP(m);
    unique_ptr<idvertice_t[]> subArvore(new idvertice_t[this->numeroDeVertices()]);

    LCP.reserve(m);

    //Lista de rótulos candidatos
    for (rotulo_t r = 0; r < m; ++r) {
        LCP.at(r).cand = r;
        LCP.at(r).prob = vetorProb[r];
    }
    //Coloca todos os vértices na subarvore
    for (const Vertice& v : this->listaVertices) {
        subArvore[v.id()] = v.id();
    }
    size_t removidos = 0;
    std::vector<rotulo_t> rotulosRemovidos;
    size_t subs = this->numeroDeVertices();
    while (subs > 1 && !LCP.empty()) {

        // //----Imprime LCP----//
        // for (size_t r = 0; r < LCP.size(); ++r) {
        //     std::cerr << "LCP: " << LCP.at(r).cand << "|" << LCP.at(r).prob;
        // }
        // std::cerr << std::endl;

        if(F.numeroDeRotulos() >= minRotulos){
            //std::cerr << "Excedeu" << std::endl;
            return F;
        }
        //Atualiza a probabilidade de Cada Rótulo
        atualizarProbabilidadesACO(LCP, vetorFeromonio, vetorHeuristico, m);

        size_t aux = LCP.size();
        // std::cerr << "Tam Atual LCP: " << aux << std::endl;
        // std::cerr << "removidos: " << removidos << std::endl;

        //Decide o Rótulo a ser inserido (Problema)
        //  Melhoria possível - se criar um mapeamento do rotulo com posição no vector de probabiliade
        //  e ordenar esse vector por probabilidade maior, a seleção de rótulo vai ser melhor (FEITO)
        // std::cerr << "Vai dar pau selecionando" << std::endl;
        // size_t idx = (aux < 2) ? 0 : selecionaRotulo(LCP, aux);
        // rotulo_t r = LCP.at(idx).cand;
        rotulo_t r = (aux < 2) ? 0 : selecionaRotulo(LCP, aux);
        // std::cerr << "rótulo r: " << r << std::endl;
        // std::cerr << "Não deu pau selecionando" << std::endl;
        // std::cerr << "Rótulo Escolhido: " << r << std::endl;

        //Adiciona todas arestas do rótulo
        for (const ArestaAux& a : this->rotulos.find(r)->second) {
            //verifica se aresta do rótulo Y já não foi colocada em rótulo X antes
            if (!estaoNaMesmaSubArvore(subArvore.get(), a.idOrigem, a.idDestino)) {
                F.fazerAresta(a.idOrigem, a.idDestino, r);
                juntarSubArvores(subArvore.get(), a.idOrigem, a.idDestino, this->numeroDeVertices());
                if (--subs <= 1) {
                    return F;
                }
            }
        }
        rotulosRemovidos.push_back(r);
        // std::cerr << "Rótulo Sendo Removido: " << r << std::endl;
        // std::cerr << "Rótulo Removidos " << std::endl;
        // for (rotulo_t s = 0; s < rotulosRemovidos.size(); s++){
        //     std::cerr << rotulosRemovidos.at(s) << "-";
        // }
        // std::cerr << std::endl;
        auto it = std::find_if(LCP.begin(), LCP.end(), [r](const lcp& element) {
            return element.cand == r;
        });
        size_t idx = std::distance(LCP.begin(), it);
        LCP.erase(LCP.begin() + idx);
        removidos++;
    }

    return F;
}



Grafo Grafo::algoritmoACO(size_t nIteracoes, size_t nFormigas, size_t bloco, double lambda_max, double tau_min, double tau_max) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
#define CICLOSMAX 0.1*nFormigas

    size_t nRotulos = this->numeroDeRotulos();
    Grafo melhorSol(this->numeroDeVertices(), this->numeroDeRotulos());
    Grafo solAux(this->numeroDeVertices(), 0);

    //Inicializa o Ferômonio
    unique_ptr<double[]> vetorFeromonio(new double[nRotulos]);
    inicializarFeromonios(vetorFeromonio, nRotulos, tau_min, tau_max);

    //Armazena o Peso Heuristico
    size_t numeroArestasDoRotulo = 0;
    size_t numeroTotalArestas = 0;
    unique_ptr<double[]> vetorHeuristico(new double[nRotulos]);
    for (rotulo_t r = 0; r < nRotulos; ++r) {
        numeroArestasDoRotulo = this->rotulos.find(r)->second.size();
        vetorHeuristico[r] = numeroArestasDoRotulo;
        numeroTotalArestas += numeroArestasDoRotulo;
        numeroArestasDoRotulo = 0;
    }

    //   //  ----Imprime vetor Heurístico----
    //     std::cerr << "Vetor Heuristico: ";
    //     for (size_t i = 0; i < nRotulos; i++) {
    //         std::cerr << vetorHeuristico[i] << ' ';
    //     }
    //     std::cerr << "------------------" << std::endl;

    bool isPrimeiraVez = true;
    size_t it;
    size_t formigasSemMelhora = 0;
    size_t iteracoesSemMelhora = 0;
    size_t lambda = 0;

    size_t menorNRotulos = nRotulos;
    std::vector<rotulo_t> rotulosDaMelhorSol;

    //Constroi a solução ao longo das Iterações
    for (it = 0; it <= nIteracoes; ++it) {
        bool itSemMelhora = true;
        //Inicializa as Probabilidades
        unique_ptr<double[]> vetorProb(new double[nRotulos]);
        inicializarProbabilidadesACO(vetorProb, vetorHeuristico, numeroTotalArestas, nRotulos);

        for(size_t fg = 1; fg <= nFormigas; ++fg){

            solAux = this->algoritmoACOHelper(vetorProb, vetorFeromonio, vetorHeuristico, nRotulos, menorNRotulos);
            //std::cout << "Melhor " << melhorSol.numeroDeRotulos() << " | " <<  "Iteração " << it << " | "  << "Formiga " << fg << " | " <<  solAux.numeroDeRotulos() << " |" << std::endl;

            if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
                isPrimeiraVez = false;
                melhorSol = std::move(solAux);
                menorNRotulos = melhorSol.numeroDeRotulos();
                formigasSemMelhora = 0;

                std::cout << "----- Nova Solução encontrada ! ------" << std::endl;
                std::cout << "Melhor " << melhorSol.numeroDeRotulos() << " | " <<  "Iteração " << it << " | "  << "Formiga " << fg << " |" << std::endl;

                size_t posicao;
                rotulosDaMelhorSol.clear();
                for (auto i : melhorSol.rotulos){
		            std::cout << i.first << " ";
                    rotulosDaMelhorSol.push_back(i.first);
                    posicao++;
                }
                std::cout << "--------------------------------------" << std::endl;
                std::cout << std::endl;

                if (itSemMelhora)
                    itSemMelhora = false;
            }
            else{
                formigasSemMelhora++;
            }

            if(formigasSemMelhora >= CICLOSMAX){
                //std::cout << "Sem Melhora - Lista vertices tentados na sol" << std::endl;
                // for (auto i : solAux.rotulos){
		        //     std::cout << i.first << " ";
                // }
                if(formigasSemMelhora >= 1.1*CICLOSMAX){
                    //Busca Local (Opcional)
                    fg = nFormigas;
                }
                //std::cout << std::endl;
            }

        }

        if (itSemMelhora)
            iteracoesSemMelhora++;
        else
            iteracoesSemMelhora = 0;

        formigasSemMelhora = 0;

        //Atualização do Feromônio
        atualizarFeromonios(vetorFeromonio, nRotulos, rotulosDaMelhorSol);

        // if (it % bloco == 0) {
        if (iteracoesSemMelhora >= bloco){
            //std::cout << "Hello Thais" << std::endl;
            //Smoothing (Opcional)
            if (iteracoesSemMelhora < bloco * 2){
                // if (lambda < lambda_max){
                lambda += lambda_max / bloco;
                // }
                //std::cerr << "Hello Thais, não melhora há " << iteracoesSemMelhora << " anos" << std::endl;
                double maxValor = *std::max_element(vetorFeromonio.get(), vetorFeromonio.get() + nRotulos);
                double minValor = *std::min_element(vetorFeromonio.get(), vetorFeromonio.get() + nRotulos);

                double smoothing = (maxValor - minValor) * lambda;

                for (size_t i = 0; i < nRotulos; ++i) {
                    vetorFeromonio[i] = vetorFeromonio[i] + smoothing;
                }
            }

            else {
                // if (iteracoesSemMelhora == bloco * 2){
                inicializarFeromonios(vetorFeromonio, nRotulos, tau_min, tau_max);
                iteracoesSemMelhora = 0;
                lambda = 0;
                // }
            }
            //std::cerr << "Hello Thais, a iteração atual e " << it << " anos" << std::endl;
        }

        //Atualiza Melhor Solução (pós smoothing)
        // if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
        //     isPrimeiraVez = false;
        //     melhorSol = std::move(solAux);
        // }

    }

    // std::cout << "Hello Thais" << std::endl;

    return melhorSol;
}