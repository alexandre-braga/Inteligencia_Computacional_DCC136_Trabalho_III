#include "Vertice.hpp"

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

const std::list<Aresta>& Vertice::listaDeAdjacencia() const
{
    return this->listaAresta;
}
