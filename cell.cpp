#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"

// This file demonstrates how to implement custom ORAM wide cell types.
// Such types can be structures of arbitrary numbers of RegAS and RegXS
// fields.  The example here imagines a cell of a binary search tree,
// where you would want the key to be additively shared (so that you can
// easily do comparisons), the pointers field to be XOR shared (so that
// you can easily do bit operations to pack two pointers and maybe some
// tree balancing information into one field) and the value doesn't
// really matter, but XOR shared is usually slightly more efficient.

// Now we use the cell in various ways.  This function is called by
// online.cpp.

void cell(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=4;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);
    run_coroutines(tio, [&tio, depth] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        Duoram<Cell> oram(tio.player(), size);
        auto A = oram.flat(tio, yield);
        Cell init;
        init.key.set(0xffffffffffffffff);
        init.pointers.set(0xeeeeeeeeeeeeeeee);
        init.value.set(0xdddddddddddddddd);
        A.init(init);
        Cell c;
        c.key.set(0x0102030405060708);
        c.pointers.set(0x1112131415161718);
        c.value.set(0x2122232425262728);
        // Explicit write
        A[0] = c;
        RegAS idx;
        // Explicit read
        Cell expl_read_c = A[0];
        printf("expl_read_c = ");
        expl_read_c.dump();
        printf("\n");
        // ORAM read
        Cell oram_read_c = A[idx];
        printf("oram_read_c = ");
        oram_read_c.dump();
        printf("\n");

        RegXS valueupdate;
        valueupdate.set(0x4040404040404040 * tio.player());
        RegXS pointersset;
        pointersset.set(0x123456789abcdef0 * tio.player());
        // Explicit update and write of individual fields
        A[1].CELL_VALUE += valueupdate;
        A[3].CELL_POINTERS = pointersset;
        // Explicit read of individual field
        RegXS pointval = A[0].CELL_POINTERS;
        printf("pointval = ");
        pointval.dump();
        printf("\n");

        idx.set(1 * tio.player());
        // ORAM read of individual field
        RegXS oram_value_read = A[idx].CELL_VALUE;
        printf("oram_value_read = ");
        oram_value_read.dump();
        printf("\n");
        valueupdate.set(0x8080808080808080 * tio.player());
        // ORAM update of individual field
        A[idx].CELL_VALUE += valueupdate;
        idx.set(2 * tio.player());
        // ORAM write of individual field
        A[idx].CELL_VALUE = valueupdate;

        c.key.set(0x0102030405060708 * tio.player());
        c.pointers.set(0x1112131415161718 * tio.player());
        c.value.set(0x2122232425262728 * tio.player());
        // ORAM update of full Cell
        A[idx] += c;
        idx.set(3 * tio.player());
        // ORAM write of full Cell
        A[idx] = c;

        printf("\n");

        if (depth < 10) {
            oram.dump();
            auto R = A.reconstruct();
            if (tio.player() == 0) {
                for(size_t i=0;i<R.size();++i) {
                    printf("\n%04lx ", i);
                    R[i].dump();
                }
                printf("\n");
            }
        }
    });
}
