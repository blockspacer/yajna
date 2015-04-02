// Yajna
//
// Written in 2015 by Martinho Fernandes <martinho.fernandes@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related
// and neighboring rights to this software to the public domain cellspacewide. This software is
// distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along with this software.
// If not, see <http://creativecommons.org/publicdomain/zero/1.0/>

// Hashlife implementation

#ifndef HLIFE_HPP
#define HLIFE_HPP

#include <array>
#include <functional>
#include <unordered_set>
#include <tuple>
#include <memory>
#include <utility>
#include <cstddef>
#include <cassert>

namespace hlife {
    // Cellspace is the set of all possible cells.
    struct cellspace;
    // A single cell (can be a macro-cell or a leaf)
    union cell;
    using cell_ptr = cell const*;
    using cell_ref = cell const&;

    // A point in space-time.
    struct point {
        int x, y, t;

        bool operator==(point const& other) { return tied() == other.tied(); }
        bool operator!=(point const& other) { return tied() != other.tied(); }

    private:
        std::tuple<int const&, int const&, int const&> tied() const { return std::tie(x, y, t); }
    };

    union cell {
    public:
        // Passkey for ctors
        struct key {
        private:
            key() = default;
            friend struct cellspace;
        };

        // Except for the two leaves, ctors should only be used by cellspace to
        // generate cells uniquely and as needed. The passkey idiom is used to
        // enforce that.
        cell(key const&, bool status) : alive(status) {}
        cell(key const&, cell_ref nw, cell_ref ne, cell_ref sw, cell_ref se)
        : q{ &nw, &ne, &sw, &se, nullptr } {}

        // Tests if a point is in a cell's light cone The only properties of
        // the cell that are needed are the extrinsic ones: its level (aka
        // size), and the coordinates of its center. This could be a non-static
        // member, and should probably be, but since it uses no intrinsic
        // properties, we can avoid the unnecessary evaluation of some cells.
        static bool in_light_cone(int level, point center, point p) {
            // Leaves have a degenerate light cone, so the point has to be at
            // the leaf's position to be in the cone.
            if(level == 0) return center == p;

            // If the point is in the past it's not in the light cone.
            if(p.t < center.t) return false;

            // We'll need the offsets of the quadrant centers soon...
            int center_offset = 1 << (level-2);

            // If the point is in the future and lies within this cell's future
            // space (the central pseudo-quadrant), we search the future cell.
            if(p.t > center.t
               && p.x >= (center.x - offset)
               && p.x <  (center.x + offset)
               && p.y >= (center.y - offset)
               && p.y <  (center.y + offset)) {
                return in_light_cone(level-1, qcenter, p);
            }

            // Otherwise we need to search for the point in a specific quadrant
            // The quadrant's center is at the same time, but at different
            // spatial coordinates that differ from this cell's center by a
            // quarter the side of this cell in both axes.
            bool north = p.x < center.x;
            bool west  = p.y < center.y;
            point qcenter {
                center.x + (west?  -1 : 1) * center_offset,
                center.y + (north? -1 : 1) * center_offset,
                center.t
            };

            // Recursively search the quadrant for this point
            return in_light_cone(level-1, qcenter, p);
        }

    private:
        // Retrieves a matching cell from the cellspace.
        // This should be the only way to obtain macro-cells.
        // We simply assume that all cells exist in cellspace
        // and none needs to be created through the ctors.
        static cell_ref get_cell(cellspace& space, cell_ref nw, cell_ref ne, cell_ref sw, cell_ref se);

        // Live and dead leaf cells from the cellspace.
        static cell_ref live(cellspace& space);
        static cell_ref dead(cellspace& space);

        // Evaluates this cell, effectively computing the future of this cell's
        // quadrants. This result is a cell one size smaller as the rest of the
        // cell depends on neighboring cells.
        cell_ref result(cellspace& space, int level) const {
            assert(level > 1); // 1-cells cannot be evaluated

            // Early exit for memoised results.
            if(q.future) return *q.future;

            // Recursive evaluation bottoms out at 2-cells since 1-cells cannot
            // be evaluated.
            if(level == 2) {
                // 2-cells are evaluated normally by counting live neighbours.
                int nnw = *q.nw->q.nw + *q.nw->q.ne + *q.ne->q.nw
                        + *q.nw->q.sw               + *q.ne->q.sw
                        + *q.sw->q.nw + *q.sw->q.ne + *q.se->q.nw;
                int nne = *q.nw->q.ne + *q.ne->q.nw + *q.ne->q.ne
                        + *q.nw->q.se               + *q.ne->q.se
                        + *q.sw->q.ne + *q.se->q.nw + *q.se->q.ne;
                int nsw = *q.nw->q.sw + *q.nw->q.se + *q.ne->q.sw
                        + *q.sw->q.nw               + *q.se->q.nw
                        + *q.sw->q.sw + *q.sw->q.se + *q.se->q.sw;
                int nse = *q.nw->q.se + *q.ne->q.sw + *q.ne->q.se
                        + *q.sw->q.ne               + *q.se->q.ne
                        + *q.sw->q.se + *q.se->q.sw + *q.se->q.se;
                q.future = &get_cell(space,
                        future_leaf(space, nnw, *q.nw->q.se),
                        future_leaf(space, nne, *q.ne->q.sw),
                        future_leaf(space, nsw, *q.sw->q.ne),
                        future_leaf(space, nse, *q.se->q.nw));
            } else {
                // n-cells are evaluated by combining the results of nine n-2-cells...
                cell_ref inw = q.nw->result(space, level-1);
                cell_ref in  = result_horizontal(space, level-1, *q.nw, *q.ne);
                cell_ref ine = q.ne->result(space, level-1);
                cell_ref iw  = result_vertical(space, level-1, *q.nw, *q.sw);
                cell_ref ix  = result_center(space, level-1);
                cell_ref ie  = result_vertical(space, level-1, *q.ne, *q.se);
                cell_ref isw = q.sw->result(space, level-1);
                cell_ref is  = result_horizontal(space, level-1, *q.sw, *q.se);
                cell_ref ise = q.se->result(space, level-1);

                // ... into four n-1-cells...
                cell_ref gnw = get_cell(space, inw, in, iw, ix);
                cell_ref gne = get_cell(space, in, ine, ix, ie);
                cell_ref gsw = get_cell(space, iw, ix, isw, is);
                cell_ref gse = get_cell(space, ix, ie, is, ise);

                // ... and then into a single n-1-cell
                q.future = &get_cell(space,
                        gnw.result(space, level-1),
                        gne.result(space, level-1),
                        gsw.result(space, level-1),
                        gse.result(space, level-1));
            }
            return *q.future;
        }

        // Evaluates the pseudo-quadrant that straddles the four quadrants in
        // the center.
        cell_ref result_center(cellspace& space, int level) const {
            return get_cell(space, *q.nw->q.se, *q.ne->q.sw, *q.sw->q.ne, *q.se->q.nw).result(space, level);
        }
        // Evaluates the pseudo-quadrant that straddles the two given quadrants
        // horizontally.
        friend cell_ref result_horizontal(cellspace& space, int level, cell_ref w, cell_ref e) {
            return cell::get_cell(space, *w.q.ne, *e.q.nw, *w.q.se, *e.q.sw).result(space, level);
        }
        // Evaluates the pseudo-quadrant that straddles the two given quadrants
        // vertically.
        friend cell_ref result_vertical(cellspace& space, int level, cell_ref n, cell_ref s) {
            return cell::get_cell(space, *n.q.sw, *n.q.se, *s.q.nw, *s.q.ne).result(space, level);
        }

        // Chooses the right next generation 0-cell for a given cell state and
        // number of live neighbours.
        static cell_ref future_leaf(cellspace& space, int neighbours, bool alive) {
            if(neighbours == 2 && alive) return live(space);
            if(neighbours == 3) return live(space);
            return dead(space);
        }

        operator bool() const { return alive; }

        // This is an immovable brick as its identity somewhat relies on its
        // address
        cell(cell const&) = delete;
        cell& operator=(cell const&) = delete;
        ~cell() = default;

        // Definition of cell equivalent for use in a hash set.
        friend struct equivalence;
        struct equivalence {
            static std::size_t combine() { return 0; }
            template <typename... Args>
            static std::size_t combine(std::size_t n, Args... args) {
                std::size_t seed = combine(args...);
                return seed ^ (n + 0x9e3779b9 + (seed<<6) + (seed>>2));
            }
            std::size_t operator()(cell_ref c) const { // hash
                std::hash<cell_ptr> h;
                return combine(h(c.q.nw), h(c.q.ne), h(c.q.sw), h(c.q.se));
            }
            bool operator()(cell_ref a, cell_ref b) const { // equivalence
                return a.q.nw == b.q.nw && a.q.ne == b.q.ne && a.q.sw == b.q.sw && a.q.se == b.q.se;
            }
        };

        // A cell is either a leaf which can be alive or not...
        bool alive;
        // ... or it is a macro-cell that links to other cells.
        struct {
            // The four cell quadrants.
            cell_ptr nw, ne, sw, se;
            // The cell obtained by evaluating the future of this cell.
            // This is mutable because it is evaluated lazily.
            mutable cell_ptr future;
        } q;
        // No tagging is used to tell the two apart; all operations know which
        // kind of cell they work on from context.
    };

    // A cellspace is the set of all possible cells.
    struct cellspace {
    public:
        cellspace() {
            // 0-cells and 1-cells cannot be generated through the recursive
            // process, so they need to be pre-seeded.

            // 0-cells (leaves) are static locals obtained from factories.
            // They don't live in the hash set as they cannot be hashed.
            auto leaf_of = [this](bool b) -> cell_ref { return b? live_cell() : dead_cell(); };

            // Pre-generate all 16 1-cells from bit patterns 0-15.
            for(int i = 0; i < 16; ++i)
                cell_with(leaf_of(i&8), leaf_of(i&4), leaf_of(i&2), leaf_of(i&1));
        }

        // The one live cell.
        cell_ref live_cell() const {
            static cell_ref live = cell(cell::key(), true);
            return live;
        }
        // The one dead cell.
        cell_ref dead_cell() const {
            static cell_ref dead = cell(cell::key(), false);
            return dead;
        }

        // Obtains a cell with the given quadrants. Cells are created lazily
        // when requested for the first time.
        cell_ref cell_with(cell_ref nw, cell_ref ne, cell_ref sw, cell_ref se) const {
            return *cells.emplace(cell::key(), nw, ne, sw, se).first;
        }

    private:
        mutable std::unordered_set<cell, cell::equivalence, cell::equivalence> cells;
    };

    cell_ref cell::live(cellspace& space) { return space.live_cell(); }
    cell_ref cell::dead(cellspace& space) { return space.dead_cell(); }
    cell_ref cell::get_cell(cellspace& space, cell_ref nw, cell_ref ne, cell_ref sw, cell_ref se) {
        return space.cell_with(nw, ne, sw, se);
    }

    struct world {
    public:
        // Generates an empty square world with sides 2^level using the given cellspace.
        world(std::shared_ptr<cellspace> space, int level)
        : space(std::move(space))
        , level(level)
        , root([this]() -> cell_ref {
                // pre-generate emptiness for all cells up to level
                cell_ptr empty = &this->space->dead_cell();
                for(int n = 1; n < this->level; ++n) {
                    // make an empty higher level cell with all four quadrants
                    // consisting of the previous level's empty cell.
                    empty = &this->space->cell_with(*empty, *empty, *empty, *empty);
                }
                // return the empty cell at the top
                return *empty;
            }()) {}

        // TODO growth
        // TODO mutation
        // TODO display

    private:
        std::shared_ptr<cellspace> space;
        int level;
        cell_ref root;
    };
}

#endif // HLIFE_HPP
