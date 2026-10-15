// Read back a BDF with ExtFEMIO's own reader to validate the format.
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include "bdf/cards.h"
#include "bdf/file.h"

namespace bc = dnvgl::extfem::bdf::cards;

int main(int argc, char **argv) {
    if (argc < 2) return 2;
    std::ifstream fs(argv[1]);
    dnvgl::extfem::bdf::input::bdf_file inp(fs);
    std::list<std::string> lines, entries;
    std::map<std::string, size_t> counts;
    size_t cards = 0, unknown = 0, errors = 0;
    while (!inp.eof()) {
        lines.clear();
        inp.get(lines);
        if (lines.empty()) break;
        std::unique_ptr<bc::__base::card> cur;
        entries.clear();
        bc::__base::card::card_split(lines, entries);
        try {
            bc::dispatch(entries, cur);
        } catch (const std::exception &e) {
            ++errors;
            std::cerr << "PARSE ERROR: " << e.what() << "\n  at: "
                      << lines.front() << "\n";
            continue;
        }
        if (!cur) continue;
        ++cards;
        auto t = cur->card_type();
        if (t == bc::types::UNKNOWN) {
            ++unknown;
            std::string name = lines.front().substr(0, 8);
            while (!name.empty() && name.back() == ' ') name.pop_back();
            counts["UNKNOWN:" + name]++;
        } else {
            std::ostringstream os; os << static_cast<int>(t);
            counts["ok:" + os.str()]++;
        }
    }
    std::cout << "cards=" << cards << " unknown=" << unknown
              << " errors=" << errors << "\n";
    for (auto &kv : counts) std::cout << "  " << kv.first << ": " << kv.second << "\n";
    return errors ? 1 : 0;
}
