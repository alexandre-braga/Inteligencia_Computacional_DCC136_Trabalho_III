#include "utils.hpp"

bool stringTok(const std::string& s, std::string& dest, char delim, size_t *posptr)
{
#define SET_POSPTR(val) do { if (posptr) { *posptr = (val); } } while (0)
#define GET_SUBSTR() ((pos == 0) ? s : s.substr(pos))

    size_t pos = (posptr) ? (*posptr) : 0;

    if (pos == s.length()) {
        return false;
    }

    if (delim == '\0') {
        dest = GET_SUBSTR();
        SET_POSPTR(s.length());
    } else {
        while (s.at(pos) == delim) {
            ++pos;
        }

        size_t end = s.find(delim, pos);

        /* Se não encontrou 'delim' */
        if (end == std::string::npos) {
            dest = GET_SUBSTR();
            SET_POSPTR(s.length());
        } else {
            dest = s.substr(pos, end - pos);
            SET_POSPTR(end + 1);
        }
    }
    return true;
}
