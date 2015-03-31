#define NONIUS_RUNNER
#include <nonius.h++>
#define private public
#include <hlife/hlife.h++>

NONIUS_BENCHMARK("generate-16384", []{
    hlife::cellspace L(16384);
    L.root.result(L, 16383);
})
