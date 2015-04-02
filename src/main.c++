#define NONIUS_RUNNER
#include <nonius.h++>
#define private public
#include <hlife/hlife.h++>

NONIUS_BENCHMARK("generate-16384", []{
    auto space = std::make_shared<hlife::cellspace>();
    hlife::world w(space, 16384);
    w.root.result(*space, 16383);
})
