#include "Vertice.hpp"

#include <algorithm>

Vertice::Vertice(idvertice_t id):
    _id(id)
{
}

idvertice_t Vertice::id() const
{
    return this->_id;
}

void Vertice::addAdjacente(idvertice_t id, rotulo_t rotulo)
{
    for (const Aresta& a : this->listaAresta) {
        if (a.idDestino() == id) {
            return;
        }
    }
    this->listaAresta.push_back(Aresta(id, rotulo));
}

void Vertice::removeAdjacente(idvertice_t id, rotulo_t rotulo)
{
    auto it = std::find(this->listaAresta.begin(), this->listaAresta.end(), Aresta(id, rotulo));
    if (it != this->listaAresta.end())
        this->listaAresta.erase(it);
}

const std::list<Aresta>& Vertice::listaDeAdjacencia() const
{
    return this->listaAresta;
}
