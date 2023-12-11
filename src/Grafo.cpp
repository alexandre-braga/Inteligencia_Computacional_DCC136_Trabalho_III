#include "Grafo.hpp"
#include "Aresta.hpp"
#include "utils.hpp"

#include <cmath>
#include <memory>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <random>

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

static inline void smoothingDeFeromonios(unique_ptr<double[]>& vetorFeromonio, size_t m, double lambda) {
    double maxValor = *std::max_element(vetorFeromonio.get(), vetorFeromonio.get() + m);
    double minValor = *std::min_element(vetorFeromonio.get(), vetorFeromonio.get() + m);
    double smoothing = (maxValor - minValor) * lambda;

    for (size_t i = 0; i < m; ++i) {
        vetorFeromonio[i] = vetorFeromonio[i] + smoothing;
    }
}

static inline void smoothingDeHeuristicas(unique_ptr<double[]>& vetorHeuristico, size_t m, double zeta) {

    double maxValorHeur = *std::max_element(vetorHeuristico.get(), vetorHeuristico.get() + m);
    double minValorHeur = *std::min_element(vetorHeuristico.get(), vetorHeuristico.get() + m);
    double smoothing2 = (maxValorHeur - minValorHeur) * zeta;

    for (size_t i = 0; i < m; ++i) {
        vetorHeuristico[i] += vetorHeuristico[i] + smoothing2;
    }
}


static inline void atualizarProbabilidadesACO(std::vector<Grafo::lcp>& LCP, unique_ptr<double[]>& vetorFeromonio,
unique_ptr<double[]>& vetorHeuristico, size_t m)
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
    std::vector<TeN> ten(m);
    ten.reserve(m);
    double somatorioDeTNs = 0;

    for (const auto& element : LCP) {
        rotulo_t r = element.cand;
        ten.at(r).id = r;
        ten.at(r).T = vetorFeromonio[r];
        ten.at(r).N = vetorHeuristico[r];
        ten.at(r).content = true;
        somatorioDeTNs += ten.at(r).T * ten.at(r).N;
    }
    for (rotulo_t r = 0; r < LCP.size(); r++){

        auto it = std::find_if(ten.begin(), ten.end(), [&LCP, r](const TeN& element) {
            return element.id == LCP.at(r).cand;
        });
        size_t itPos = std::distance(ten.begin(), it);

        if(ten.at(itPos).content){
            LCP.at(r).prob = ten.at(r).T * ten.at(r).N / somatorioDeTNs;
        }
        else{
            LCP.at(r).prob = -1;
        }
    }
}
static size_t selecionaRotulo(std::vector<Grafo::lcp> LCP, size_t tamAtual)
{
    std::mt19937 rng(std::random_device{}());

    double rngGenerated = 0.01 + (double) (rng() % 2900) / 10000;

    double aux = 0;
    LCP.resize(tamAtual);
    std::sort(LCP.begin(),LCP.end(), [](const Grafo::lcp &x, const Grafo::lcp &y){ return (x.prob > y.prob);});

    for (size_t i = 0; i < tamAtual; i++) {
        aux += LCP.at(i).prob;
        if (rngGenerated <= aux) {
            return LCP.at(i).cand;
        }
    }
    return LCP.at(tamAtual-1).cand;
    // return LCP.at(rngGenerated).cand;
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

    size_t subs = this->numeroDeVertices();
    while (subs > 1 && !LCP.empty()) {

        if(F.numeroDeRotulos() >= minRotulos){
            //std::cerr << "Excedeu" << std::endl;
            return F;
        }

        //Atualiza a probabilidade de Cada Rótulo
        atualizarProbabilidadesACO(LCP, vetorFeromonio, vetorHeuristico, m);

        size_t aux = LCP.size();

        // if(aux == m){
        //     std::cerr << "Final Probs: " << std::endl;
        //     for (size_t i = 0; i < LCP.size(); i++) {
        //         std::cerr << LCP.at(i).prob << " ";
        //     }
        //     std::cerr << std::endl;
        // }

        //Decide o Rótulo a ser inserido
        rotulo_t r = (aux < 2) ? 0 : selecionaRotulo(LCP, aux);

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

        auto it = std::find_if(LCP.begin(), LCP.end(), [r](const lcp& element) {
            return element.cand == r;
        });
        size_t idx = std::distance(LCP.begin(), it);
        LCP.erase(LCP.begin() + idx);

    }


    return F;
}


Grafo Grafo::algoritmoACOSmoothing(size_t nIteracoes, size_t nFormigas, size_t bloco, double lambda_max, double zeta_max, double tau_min, double tau_max) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
#define FORMIGASMAX 0.1*nFormigas

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

    bool isPrimeiraVez = true;
    size_t it;
    size_t formigasSemMelhora = 0;
    size_t iteracoesSemMelhora = 0;
    double lambda = 0;
    double zeta = 0;
    size_t menorNRotulos = nRotulos;
    std::vector<rotulo_t> rotulosDaMelhorSol;

    //Constroi a solução ao longo das Iterações
    for (it = 0; it <= nIteracoes; ++it) {
        bool itSemMelhora = true;

        //Inicializa as Probabilidades
        unique_ptr<double[]> vetorProb(new double[nRotulos]);
        size_t proporcaoZeta = 0;
        for (rotulo_t r = 0; r < nRotulos; ++r) {
            proporcaoZeta += vetorHeuristico[r];
        }
        inicializarProbabilidadesACO(vetorProb, vetorHeuristico, proporcaoZeta, nRotulos);

        for(size_t fg = 1; fg <= nFormigas; ++fg){

            solAux = this->algoritmoACOHelper(vetorProb, vetorFeromonio, vetorHeuristico, nRotulos, menorNRotulos);
            //std::cout << "Iteração " << it << " | "  << "Formiga " << fg << " | " <<  solAux.numeroDeRotulos() << " |" << std::endl;

            if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
                isPrimeiraVez = false;
                melhorSol = std::move(solAux);
                menorNRotulos = melhorSol.numeroDeRotulos();
                formigasSemMelhora = 0;

                // std::cout << "----- Nova Solução encontrada ! ------" << std::endl;
                // std::cout << "Melhor " << melhorSol.numeroDeRotulos() << " | " <<  "Iteração " << it << " | "  << "Formiga " << fg << " |" << std::endl;

                size_t posicao;
                rotulosDaMelhorSol.clear();
                for (auto i : melhorSol.rotulos){
		            // std::cout << i.first << " ";
                    rotulosDaMelhorSol.push_back(i.first);
                    posicao++;
                }
                // std::cout << "--------------------------------------" << std::endl;
                // std::cout << std::endl;

                if (itSemMelhora)
                    itSemMelhora = false;
            }
            else{
                formigasSemMelhora++;
            }

            if(formigasSemMelhora >= FORMIGASMAX){

                // std::cout << "Sem Melhora - Lista vertices tentados na sol" << std::endl;
                // for (auto i : solAux.rotulos){
		        //     std::cout << i.first << " ";
                // }
                // std::cout << std::endl;

                //Busca Local
                if(formigasSemMelhora >= 1.1*FORMIGASMAX){
                    fg = nFormigas;
                }
            }

        }
        formigasSemMelhora = 0;

        if (itSemMelhora)
            iteracoesSemMelhora++;
        else
            iteracoesSemMelhora = 0;

        //Atualização do Feromônio
        atualizarFeromonios(vetorFeromonio, nRotulos, rotulosDaMelhorSol);

        if (iteracoesSemMelhora >= bloco){
            if (iteracoesSemMelhora < bloco * 2){

                // if (lambda < lambda_max){
                lambda += lambda_max / bloco;
                // }
                smoothingDeFeromonios(vetorFeromonio, nRotulos, lambda);
                smoothingDeHeuristicas(vetorHeuristico, nRotulos, zeta);

                // std::cerr << "Hello Thais, LAMBDA atual: " << lambda << std::endl;
                // std::cerr << "Hello Thais, ZETA atual: " << zeta << std::endl;
            }
            else {

                inicializarFeromonios(vetorFeromonio, nRotulos, tau_min, tau_max);
                iteracoesSemMelhora = 0;
                lambda = 0;
                zeta += zeta_max / bloco;

                // std::cerr << "Hello Thais, bora começar de novo ?" << std::endl;
            }
        }
        // std::cerr << "Hello Thais, se passaram " << it << " anos" << std::endl;
        // std::cerr << "------------------------------------------" << std::endl;
        // std::cerr << std::endl;
    }

    // std::cout << "Hello Thais, acabamo a execução" << std::endl;

    return melhorSol;
}

Grafo Grafo::algoritmoACO(size_t nIteracoes, size_t nFormigas, double tau_min, double tau_max) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
#define FORMIGASMAX 0.1*nFormigas

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

    bool isPrimeiraVez = true;
    size_t it;
    size_t menorNRotulos = nRotulos;
    std::vector<rotulo_t> rotulosDaMelhorSol;
    size_t formigasSemMelhora = 0;

    //Constroi a solução ao longo das Iterações
    for (it = 0; it <= nIteracoes; ++it) {

        //Inicializa as Probabilidades
        unique_ptr<double[]> vetorProb(new double[nRotulos]);
        inicializarProbabilidadesACO(vetorProb, vetorHeuristico, numeroTotalArestas, nRotulos);

        for(size_t fg = 1; fg <= nFormigas; ++fg){

            solAux = this->algoritmoACOHelper(vetorProb, vetorFeromonio, vetorHeuristico, nRotulos, menorNRotulos);

            if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
                isPrimeiraVez = false;
                melhorSol = std::move(solAux);
                menorNRotulos = melhorSol.numeroDeRotulos();

                size_t posicao;
                rotulosDaMelhorSol.clear();
                for (auto i : melhorSol.rotulos){
		            rotulosDaMelhorSol.push_back(i.first);
                    posicao++;
                }

            }
            else{
                formigasSemMelhora++;
            }
            if(formigasSemMelhora >= FORMIGASMAX){
                fg = nFormigas;
            }

        }
        formigasSemMelhora = 0;

        atualizarFeromonios(vetorFeromonio, nRotulos, rotulosDaMelhorSol);

    }

    return melhorSol;
}

Grafo Grafo::algoritmoACOV2(size_t nIteracoes, size_t nFormigas, double tau_min, double tau_max) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
#define FORMIGASMAX 0.1*nFormigas

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

    bool isPrimeiraVez = true;
    size_t it;
    size_t menorNRotulos = nRotulos;
    std::vector<rotulo_t> rotulosDaMelhorFormiga;
    size_t formigasSemMelhora = 0;

    //Constroi a solução ao longo das Iterações
    for (it = 0; it <= nIteracoes; ++it) {

        //Inicializa as Probabilidades
        unique_ptr<double[]> vetorProb(new double[nRotulos]);
        inicializarProbabilidadesACO(vetorProb, vetorHeuristico, numeroTotalArestas, nRotulos);

        bool isPrimeiraFormiga = true;
        Grafo melhorFormiga(this->numeroDeVertices(), this->numeroDeRotulos());
        size_t menorNRotulosFormigas = nRotulos;

        for(size_t fg = 1; fg <= nFormigas; ++fg){

            solAux = this->algoritmoACOHelper(vetorProb, vetorFeromonio, vetorHeuristico, nRotulos, menorNRotulosFormigas);

            if (isPrimeiraFormiga || solAux.numeroDeRotulos() < melhorFormiga.numeroDeRotulos()) {
                isPrimeiraFormiga = false;
                melhorFormiga = std::move(solAux);
                menorNRotulosFormigas = melhorFormiga.numeroDeRotulos();
                size_t posicao;
                rotulosDaMelhorFormiga.clear();
                for (auto i : melhorSol.rotulos){
                    rotulosDaMelhorFormiga.push_back(i.first);
                    posicao++;
                }
            }
            else{
                formigasSemMelhora++;
            }
            if(formigasSemMelhora >= FORMIGASMAX){
                fg = nFormigas;
            }

        }
        formigasSemMelhora = 0;
        atualizarFeromonios(vetorFeromonio, nRotulos, rotulosDaMelhorFormiga);

        if (isPrimeiraVez || melhorFormiga.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
            isPrimeiraVez = false;
            melhorSol = std::move(melhorFormiga);
            menorNRotulos = melhorSol.numeroDeRotulos();
        }

    }

    return melhorSol;
}

/*-------------------------------------------------------------*/
/*----------------------- Busca Local -------------------------*/
/*-------------------------------------------------------------*/

// Remoção aleatória de n rótulos da solução
// Inserção de m <= n rótulos que maximizam o número de componentes conexas
// V1: Vizinhança gerada pela remoção de 1 rótulo
// V2: Vizinhança gerada pela remoção de 2 rótulos
// V3: Vizinhança gerada pela remoção de 3 rótulos

static void separarSubArvores(idvertice_t subArvore[], idvertice_t u, idvertice_t v)
{
    subArvore[u] = u;
    subArvore[v] = v;
}

void Grafo::desfazerAresta(idvertice_t id1, idvertice_t id2, rotulo_t rotulo)
{
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

        v1->removeAdjacente(id2, rotulo);
        v2->removeAdjacente(id1, rotulo);
    } else {
        Vertice *v = this->getVerticeById(id1);

        TEST_VERTICE(v, id1);
        v->removeAdjacente(id1, rotulo);
    }
}

Grafo Grafo::removeKRotulos(Grafo s, size_t k, unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t>& listCand) const{
    Grafo s1(this->numeroDeVertices(), 0);
    s1 = s;

    std::vector<rotulo_t> rotulosS;
    std::vector<rotulo_t> rotulosSRemovidos;

    for (auto i : s1.rotulos){
        rotulosS.push_back(i.first);
    }

    size_t i = 0;
    while (i < k) {
        rotulo_t indiceRotuloSorteado = std::rand() % rotulosS.size();
        auto it = std::find(rotulosSRemovidos.begin(), rotulosSRemovidos.end(), rotulosS.at(indiceRotuloSorteado));
        if (it == rotulosSRemovidos.end()){
            i += 1;
            for (const ArestaAux& a : s1.rotulos.find(rotulosS.at(indiceRotuloSorteado))->second) {
                s1.desfazerAresta(a.idOrigem, a.idDestino, rotulosS.at(indiceRotuloSorteado));
                separarSubArvores(subArvore.get(), a.idOrigem, a.idDestino);
            }
            auto itRotulo = s1.rotulos.find(rotulosS.at(indiceRotuloSorteado));
            listCand.push_back(rotulosS.at(indiceRotuloSorteado));
            s1.rotulos.erase(itRotulo);
        }
    }

    return s1;
}

Grafo Grafo::reparaSolucao(Grafo s1, std::unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t>& listCand, size_t k) const{

    Grafo sAux(this->numeroDeVertices(), 0);
    sAux = s1;

    std::vector<idvertice_t> subArvoresDif;
    size_t numSubArvores = 0;

    for (idvertice_t id = 0; id < this->numeroDeVertices(); ++id) {
        auto idAux = std::find(subArvoresDif.begin(), subArvoresDif.end(), subArvore[id]);
        if (idAux == subArvoresDif.end()) {
            subArvoresDif.push_back(subArvore[id]);
            numSubArvores += 1;
        }
    }

    size_t idx;
    size_t tamListCand;
    float qualidade = 1;
    size_t itSemMelhora = 0;
    bool CONTINUE = true;

    while (numSubArvores > 1 && !listCand.empty() && CONTINUE) {
        idx = 0;
        tamListCand = listCand.size();

        std::cout << " oi " << idx << " " << tamListCand << std::endl;

        while (idx < tamListCand){

            std::cout << " hello " << idx << " " << tamListCand << std::endl;

            size_t numSAaux = numSubArvores;
            std::vector<idvertice_t> vert;
            rotulo_t r = listCand[idx];

            // std::cout << "lsatual " << idx << " " << r << std::endl;

            for (const ArestaAux& a : this->rotulos.find(r)->second) {
                // std::cout << "vai entrar?" << std::endl;
                // std::cout << a.idOrigem << " " << a.idDestino << std::endl;
                if (a.idOrigem != a.idDestino && !estaoNaMesmaSubArvore(subArvore.get(), a.idOrigem, a.idDestino)) {
                    // std::cout << a.idOrigem << " " << a.idDestino << std::endl;
                    auto id1 = std::find(vert.begin(), vert.end(), a.idOrigem);
                    auto id2 = std::find(vert.begin(), vert.end(), a.idDestino);

                    // std::cout << "vert " << std::endl;
                    // for (auto l = vert.begin(); l != vert.end(); l++){
                    //     std::cout << *l << " ";
                    // }

                    bool id1Entra = false;
                    bool id2Entra = false;

                    if (id1 == vert.end()){
                        id1Entra = true;
                    }

                    if (id2 == vert.end()){
                        id2Entra = true;
                    }

                    if (id1Entra || id2Entra){
                        numSAaux -= 1;
                    }

                    if (id1Entra){
                        vert.push_back(a.idOrigem);
                    }

                    if (id2Entra){
                        vert.push_back(a.idDestino);
                    }
                }
            }

            // std::cout << "Vert" << std::endl;
            // for (auto l = vert.begin(); l != vert.end(); l++){
            //     std::cout << *l << " ";
            // }

            // std::cout << numSubArvores << " - " << numSAaux << " - " << k << std::endl;

            if (numSAaux <= numSubArvores * qualidade){
                // std::cout << "Entrou" << std::endl;
                for (const ArestaAux& a : this->rotulos.find(r)->second) {
                    if (!estaoNaMesmaSubArvore(subArvore.get(), a.idOrigem, a.idDestino)) {
                        sAux.fazerAresta(a.idOrigem, a.idDestino, r);
                        juntarSubArvores(subArvore.get(), a.idOrigem, a.idDestino, this->numeroDeVertices());
                        if (--numSubArvores <= 1) {
                            // std::cout << "Achou (o que eu n garanto)" << std::endl;
                            // for (const Vertice& v : this->listaVertices) {
                            //     std::cout << subArvore[v.id()] << " | ";
                            // }
                            return sAux;
                        }
                    }
                }

                // std::cout << "nsa " << numSubArvores << std::endl;

                // std::cout << "\nSubarvore de cada vertice" << std::endl;
                // for (const Vertice& v : this->listaVertices) {
                //     std::cout << v.id() << " - " << subArvore[v.id()] << " | ";
                // }

                // std::cout << "\nLC" << std::endl;
                // for (auto l = listCand.begin(); l != listCand.end(); l++){
                //     std::cout << *l << " ";
                // }

                listCand.erase(listCand.begin() + idx);

                // std::cout << "\nLC" << std::endl;
                // for (auto l = listCand.begin(); l != listCand.end(); l++){
                //     std::cout << *l << " ";
                // }

                break;
            }

            idx++;

            if (itSemMelhora < qualidade * 10 && idx >= tamListCand){
                itSemMelhora += 1;
                // qualidade += 0.1;
            }
            else {
                if (itSemMelhora >= qualidade * 10){
                    CONTINUE = false;
                }
            }
        }
    }

    // if (listCand.empty()){
    //     std::cout << "com ctz deu ruim" << std::endl;
    // }
    // else{
    //     std::cout << "\nLC" << std::endl;
    //     for (auto l = listCand.begin(); l != listCand.end(); l++){
    //         std::cout << *l << " ";
    //     }
    // }

    // std::cout << "num " << numSubArvores << " " << idx << std::endl;

    // std::cout << "Acho que deu ruim" << std::endl;
    // for (const Vertice& v : this->listaVertices) {
    //     std::cout << subArvore[v.id()] << " | ";
    // }

    std::cout << "\nSubArvore: " << std::endl;
    for (const Vertice& v : this->listaVertices) {
        std::cout << subArvore[v.id()] << " ";
    }

    return sAux;
}

Grafo Grafo::VNS(Grafo solAtual, unique_ptr<idvertice_t[]>& subArvore, std::vector<rotulo_t> listCand) const
{
#define NUM_VIZINHANCAS_VNS 4
#define IT_MAX_SEM_MELHORA 10
    Grafo s(this->numeroDeVertices(), 0);
    Grafo s1(this->numeroDeVertices(), 0);

    // s -> melhor solução
    s = solAtual;

    size_t itSemMelhora = 0;
    size_t k;

    do {
        bool melhoraIt = false;
        k = 1;
        while (k <= NUM_VIZINHANCAS_VNS) {

            std::cout << "Remove" << std::endl;
            // Pertubação
            s1 = this->removeKRotulos(s, k, subArvore, listCand);
            std::cout << "Repara" << std::endl;
            // Reparo
            s1 = this->reparaSolucao(s1, subArvore, listCand, k);

            // Busca Local (talvez depois (ou não :()))
            std::cout << "Move" << std::endl;
            // Mover

            std::cout << "\nSubArvore: " << std::endl;
            for (const Vertice& v : this->listaVertices) {
                std::cout << subArvore[v.id()] << " ";
            }

            std::cout << "s1 " << s1.numeroDeRotulos() << " s " << s.numeroDeRotulos() << std::endl;
            if (s1.numeroDeRotulos() < s.numeroDeRotulos()) {
                s = std::move(s1);
                k = 1;
                melhoraIt = true;
            }
            else {
                k += 1;
            }

            std::cout << itSemMelhora << " " << k << std::endl;
        }

        if (melhoraIt){
            itSemMelhora = 0;
        }
        else{
            itSemMelhora++;
        }
    } while (itSemMelhora < IT_MAX_SEM_MELHORA);

    return s;
}

Grafo Grafo::algoritmoACOHelperBL(unique_ptr<double[]>& vetorFeromonio, unique_ptr<double[]>& vetorHeuristico, size_t m, size_t minRotulos, unique_ptr<idvertice_t[]>& subArvore, std::vector<lcp>& LCP) const
{

    Grafo F(this->numeroDeVertices(), 0);

    size_t subs = this->numeroDeVertices();
    while (subs > 1 && !LCP.empty()) {

        if(F.numeroDeRotulos() >= minRotulos){
            return F;
        }

        atualizarProbabilidadesACO(LCP, vetorFeromonio, vetorHeuristico, m);

        size_t aux = LCP.size();
        rotulo_t r = (aux < 2) ? 0 : selecionaRotulo(LCP, aux);

        for (const ArestaAux& a : this->rotulos.find(r)->second) {
            if (!estaoNaMesmaSubArvore(subArvore.get(), a.idOrigem, a.idDestino)) {
                F.fazerAresta(a.idOrigem, a.idDestino, r);
                juntarSubArvores(subArvore.get(), a.idOrigem, a.idDestino, this->numeroDeVertices());
                if (--subs <= 1) {
                    return F;
                }
            }
        }

        auto it = std::find_if(LCP.begin(), LCP.end(), [r](const lcp& element) {
            return element.cand == r;
        });
        size_t idx = std::distance(LCP.begin(), it);
        LCP.erase(LCP.begin() + idx);
    }

    return F;
}

Grafo Grafo::algoritmoACOBuscaLocal(size_t nIteracoes, size_t nFormigas, double tau_min, double tau_max) const
{
#define PESO(r) (-this->rotulos.find(r)->second.size())
#define FORMIGASMAX 0.1*nFormigas

    size_t nRotulos = this->numeroDeRotulos();
    Grafo melhorSol(this->numeroDeVertices(), this->numeroDeRotulos());
    Grafo solAux(this->numeroDeVertices(), 0);
    Grafo solBL(this->numeroDeVertices(), 0);

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

    bool isPrimeiraVez = true;
    size_t it;
    size_t menorNRotulos = nRotulos;
    std::vector<rotulo_t> rotulosDaMelhorSol;
    size_t formigasSemMelhora = 0;

    //Constroi a solução ao longo das Iterações
    for (it = 0; it <= nIteracoes; ++it) {

        //Inicializa as Probabilidades
        unique_ptr<double[]> vetorProb(new double[nRotulos]);
        inicializarProbabilidadesACO(vetorProb, vetorHeuristico, numeroTotalArestas, nRotulos);

        for(size_t fg = 1; fg <= nFormigas; ++fg){

            std::vector<lcp> LCP(nRotulos);
            unique_ptr<idvertice_t[]> subArvore(new idvertice_t[this->numeroDeVertices()]);

            LCP.reserve(nRotulos);

            //Lista de rótulos candidatos
            for (rotulo_t r = 0; r < nRotulos; ++r) {
                LCP.at(r).cand = r;
                LCP.at(r).prob = vetorProb[r];
            }
            //Coloca todos os vértices na subarvore
            for (const Vertice& v : this->listaVertices) {
                subArvore[v.id()] = v.id();
            }

            solAux = this->algoritmoACOHelperBL(vetorFeromonio, vetorHeuristico, nRotulos, menorNRotulos, subArvore, LCP);

            std::vector<rotulo_t> ListaRestante;
            ListaRestante.reserve(nRotulos);
            for (const auto& element : LCP) {
                rotulo_t r = element.cand;
                ListaRestante.push_back(r);
            }

            if (isPrimeiraVez || solAux.numeroDeRotulos() < melhorSol.numeroDeRotulos()) {
                isPrimeiraVez = false;
                melhorSol = std::move(solAux);
                menorNRotulos = melhorSol.numeroDeRotulos();
            }
            else{
                formigasSemMelhora++;
            }
            if(formigasSemMelhora >= FORMIGASMAX){
                std::cout << "\nVish, vamo tentar busca local agr alexandre " << std::endl;
                solBL = this->VNS(solAux, subArvore, ListaRestante);
                if (solBL.numeroDeRotulos() < melhorSol.numeroDeRotulos()){
                    melhorSol = std::move(solBL);
                    formigasSemMelhora = 0;
                }
                fg = nFormigas;
            }

        }
        formigasSemMelhora = 0;

        atualizarFeromonios(vetorFeromonio, nRotulos, rotulosDaMelhorSol);

    }

    return melhorSol;
}