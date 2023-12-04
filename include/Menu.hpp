#ifndef MENU_HPP
#define MENU_HPP

#include <fstream>

#include "Grafo.hpp"

class Menu {
    public:
        static void run(const Grafo& g, std::ofstream& out);
};

#endif // MENU_HPP
