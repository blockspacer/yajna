// Yajna
//
// Written in 2015 by Martinho Fernandes <martinho.fernandes@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related
// and neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along with this software.
// If not, see <http://creativecommons.org/publicdomain/zero/1.0/>

// Hashlife implementation

#ifndef HLIFE_HPP
#define HLIFE_HPP

namespace hlife {
    struct world;

    union cell {
    private:
        enum { nw, ne, sw, se };

    public:
        cell const& result(world& w, int level) const {
            if(level == 2) {
                auto nnw = q[nw]->q[nw] + q[nw]->q[ne] + q[ne]->q[nw]
                         + q[nw]->q[sw]                + q[ne]->q[sw]
                         + q[sw]->q[nw] + q[sw]->q[ne] + q[se]->q[nw]
                         ;
                auto nne = q[nw]->q[ne] + q[ne]->q[nw] + q[ne]->q[ne]
                         + q[nw]->q[se]                + q[ne]->q[se]
                         + q[sw]->q[ne] + q[se]->q[nw] + q[se]->q[ne]
                         ;
                auto nsw = q[nw]->q[sw] + q[nw]->q[se] + q[ne]->q[sw]
                         + q[sw]->q[nw]                + q[se]->q[nw]
                         + q[sw]->q[sw] + q[sw]->q[se] + q[se]->q[sw]
                         ;
                auto nse = q[nw]->q[se] + q[ne]->q[sw] + q[ne]->q[se]
                         + q[sw]->q[ne]                + q[se]->q[ne]
                         + q[sw]->q[se] + q[se]->q[sw] + q[se]->q[se]
                         ;
                return get_cell(w,
                        get_leaf(nnw, q[nw]->q[se]),
                        get_leaf(nne, q[ne]->q[sw]),
                        get_leaf(nsw, q[sw]->q[ne]),
                        get_leaf(nse, q[se]->q[nw]));
            } else if(future) return *future;
            else throw nullptr; // TODO recurse
        }
    private:
        cell const& get_leaf(int neighbours, bool alive) {
            static cell live = cell(true);
            static cell dead = cell(false);
            if(neighbours == 2 && alive) return live;
            if(neighbours == 3) return live;
            return dead;
        }

        cell const& get_cell(world&, cell const& cnw, cell const& cne, cell const& csw, cell const& cse) {
            // TODO find existing in world, or add to world
            return cell(cnw, cne, csw, cse);
        }

        operator bool() const { return alive; }

        cell(bool status) : alive(status) {}
        cell(cell const& cnw, cell const& cne, cell const& csw, cell const& cse)
        : q{ &cnw, &cne, &csw, &cse }, future(nullptr) {}

        struct {
            cell const* q[4];
            mutable cell* future;
        };
        bool alive;
    };
}

#endif // HLIFE_HPP

