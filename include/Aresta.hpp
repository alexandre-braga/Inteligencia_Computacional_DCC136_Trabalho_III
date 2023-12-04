#ifndef ARESTA_HPP
#define ARESTA_HPP

#include <string>

#define stringToRotulo stoul

using rotulo_t = unsigned;
using idvertice_t = size_t;

class Aresta {
    private:
        idvertice_t _idDestino;
        rotulo_t _rotulo;

    public:
        Aresta(idvertice_t idVertice, rotulo_t rotulo);

        rotulo_t rotulo() const;
        idvertice_t idDestino() const;
};

#endif // ARESTA_HPP
