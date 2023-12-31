#include "Aresta.hpp"

Aresta::Aresta(idvertice_t idVertice, rotulo_t rotulo):
    _idDestino(idVertice),
    _rotulo(rotulo)
{
}

rotulo_t Aresta::rotulo() const
{
    return this->_rotulo;
}

idvertice_t Aresta::idDestino() const
{
    return this->_idDestino;
}

bool Aresta::operator==(const Aresta& other) const {
    return (_idDestino == other._idDestino) && (_rotulo == other._rotulo);
}