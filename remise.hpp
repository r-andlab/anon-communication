#ifndef __REMISE_HPP__
#define __REMISE_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"
#include "duoram.hpp"

 
class Writer{

    public:
    
    uint64_t * access_list;
    uint64_t * proofs;

    Writer(size_t list_size) {
        access_list = new uint64_t[list_size];
        proofs = new uint64_t[list_size];
    } 
};




class RemiseSAM {
private:
    Duoram < RegAS > oram;
    Duoram < RegAS > tokendb;
    Duoram < RegAS > proofdb;

    size_t MAX_SIZE;
    size_t num_items;

    // Basic restore heap property at a secret shared index
    // Takes in as an input the XOR shares of the index at which
    // the heap property has to be restored
    // Returns the XOR shares of the index of the smaller child
    RegXS restore_heap_property(MPCTIO &tio, yield_t & yield, RegXS index);

    // Optimized restore heap property at a secret shared index
    // Takes in as an input the XOR shares of the index at which
    // the heap property has to be restored
    // Returns the XOR shares of the index of the smaller child and
    // comparison between the left and right child
    std::pair<RegXS, RegBS> restore_heap_property_optimized(MPCTIO &tio, yield_t & yield, RegXS index, size_t layer, typename Duoram<RegAS>::template OblivIndex<RegXS,3> oidx);

    // Restore heap property at an index in clear
    // Takes in as an input the index (in clear) at which
    // the heap property has to be restored
    // Returns the XOR shares of the index of the smaller child and
    // comparison between the left and right child
    std::pair<RegXS, RegBS> restore_heap_property_at_explicit_index(MPCTIO &tio, yield_t & yield,  size_t index);

public:

    RemiseSAM(int player_num, size_t size) : oram(player_num, size), tokendb(player_num, size), proofdb(player_num, size) {};
    
    void initialize_proofDB(MPCTIO &tio, yield_t & yield, size_t dbsize, Writer writer);

    void post_message(MPCTIO &tio, yield_t & yield, RegAS ind, RegAS val, RegAS proof);

    RegAS read_message(MPCTIO& tio, yield_t& yield, RegAS ind, RegAS proof_submitted);

    void provide_or_remove_access(MPCTIO& tio, yield_t& yield, RegAS ind, RegAS proof_submitted);

    void publish_bulletin_board(MPCTIO &tio, yield_t & yield, size_t dbsize);
};
     
void Remise(MPCIO &mpcio, const PRACOptions &opts, char **args);

#endif
