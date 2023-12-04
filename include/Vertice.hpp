#ifndef VERTICE_HPP
#define VERTICE_HPP

#include <list>
#include <string>

#include "Aresta.hpp"

class Vertice {
    private:
        std::list<Aresta> listaAresta;
        idvertice_t _id;

    public:
        Vertice(idvertice_t id);

        idvertice_t id() const;
        const std::list<Aresta>& listaDeAdjacencia() const;

        void addAdjacente(idvertice_t id, rotulo_t rotulo);
};

#endif // VERTICE_HPP
