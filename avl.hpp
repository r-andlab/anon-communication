#ifndef __AVL_HPP__
#define __AVL_HPP__

#include <optional>
#include <math.h>
#include <stdio.h>
#include <string>
#include "types.hpp"
#include "duoram.hpp"
#include "cdpf.hpp"
#include "mpcio.hpp"
#include "options.hpp"
#include "bst.hpp"

/*
    Macro definitions:

    AVL_RANDOMIZE_INSERTS: Randomize keys of items inserted in the unit
    tests. When turned off, items with incremental keys are inserted

    AVL_DEBUG: General debug flag

    AVL_DEBUG_BB: Debug flag for balance bit computations
*/

// #define AVL_RANDOMIZE_INSERTS
// #define AVL_DEBUG
// #define AVL_DEBUG_BB

/*
  For AVL tree we'll treat the pointers fields as:
  < L_ptr (31 bits), R_ptr (31 bits), bal_L (1 bit), bal_R (1 bit)>
  Where L_ptr and R_ptr are pointers to the left and right child respectively,
  and bal_L and bal_R are the balance bits.

  Consequently AVL has its own versions of extract and set pointers for its children.
*/

#define AVL_PTR_SIZE 31

inline int AVL_TTL(size_t n) {
    if(n==0) {
        return 0;
    } else if (n==1) {
        return 1;
    } else {
        double logn = log2(n);
        double TTL = 1.44 * logn;
        return (int(ceil(TTL)));
    }
}

inline RegXS getAVLLeftPtr(RegXS pointer){
    return (pointer>>33);
}

inline RegXS getAVLRightPtr(RegXS pointer){
    return ((pointer&(0x00000001FFFFFFFC))>>2);
}

inline void setAVLLeftPtr(RegXS &pointer, RegXS new_ptr){
    pointer&=(0x00000001FFFFFFFF);
    pointer+=(new_ptr<<33);
}

inline void setAVLRightPtr(RegXS &pointer, RegXS new_ptr){
    pointer&=(0xFFFFFFFE00000003);
    pointer+=(new_ptr<<2);
}

inline RegBS getLeftBal(RegXS pointer){
    RegBS bal_l;
    bool bal_l_bit = ((pointer.share() & (0x0000000000000002))>>1);
    bal_l.set(bal_l_bit);
    return bal_l;
}

inline RegBS getRightBal(RegXS pointer){
    RegBS bal_r;
    bool bal_r_bit = (pointer.share() & (0x0000000000000001));
    bal_r.set(bal_r_bit);
    return bal_r;
}

inline void setLeftBal(RegXS &pointer, RegBS bal_l){
    value_t temp_ptr = pointer.share();
    temp_ptr&=(0xFFFFFFFFFFFFFFFD);
    temp_ptr^=((value_t)(bal_l.share()<<1));
    pointer.set(temp_ptr);
}

inline void setRightBal(RegXS &pointer, RegBS bal_r){
    value_t temp_ptr = pointer.share();
    temp_ptr&=(0xFFFFFFFFFFFFFFFE);
    temp_ptr^=((value_t)(bal_r.share()));
    pointer.set(temp_ptr);
}

inline void dumpAVL(Node n) {
    RegBS left_bal, right_bal;
    left_bal = getLeftBal(n.pointers);
    right_bal = getRightBal(n.pointers);
    printf("[%016lx %016lx(L:%ld, R:%ld) %d %d %016lx]", n.key.share(), n.pointers.share(),
          getAVLLeftPtr(n.pointers).xshare, getAVLRightPtr(n.pointers).xshare,
          left_bal.share(), right_bal.share(), n.value.share());
}

struct avl_del_return {
    // Flag to indicate if the key this deletion targets requires a successor swap
    RegBS F_ss;
    // Pointers to node to be deleted that would be replaced by successor node
    RegXS N_d;
    // Pointers to successor node that would replace deleted node
    RegXS N_s;
    // F_r: Flag for updating child pointer with returned pointer
    RegBS F_r;
    RegXS ret_ptr;
};

struct avl_insert_return {
  RegXS gp_node; // grandparent node
  RegXS p_node; // parent node
  RegXS c_node; // child node

  // Direction bits: 0 = Left, 1 = Right
  RegBS dir_gpp; // Direction bit from grandparent to parent node
  RegBS dir_pc; // Direction bit from p_node to c_node
  RegBS dir_cn; // Direction bit from c_node to new_node

  RegBS imbalance;
};

class AVL {
  private:
    Duoram<Node> oram;
    RegXS root;

    size_t num_items = 0;
    size_t cur_max_index = 0;
    size_t MAX_SIZE;
    int MAX_DEPTH;
    bool OPTIMIZED;

    std::vector<RegXS> empty_locations;

    std::tuple<RegBS, RegBS, RegXS, RegBS> insert(MPCTIO &tio, yield_t &yield, RegXS ptr,
        RegXS ins_addr, RegAS ins_key, Duoram<Node>::Flat &A, int TTL, RegBS isDummy,
        avl_insert_return &ret);

    void rotate(MPCTIO &tio, yield_t &yield, RegXS &gp_pointers, RegXS p_ptr,
        RegXS &p_pointers, RegXS c_ptr, RegXS &c_pointers, RegBS dir_gpp,
        RegBS dir_pc, RegBS isNotDummy, RegBS F_gp);

    std::tuple<RegBS, RegBS, RegBS, RegBS> updateBalanceIns(MPCTIO &tio, yield_t &yield,
        RegBS bal_l, RegBS bal_r, RegBS bal_upd, RegBS child_dir);

    void updateChildPointers(MPCTIO &tio, yield_t &yield, RegXS &left, RegXS &right,
          RegBS c_prime, const avl_del_return &ret_struct);

    void fixImbalance(MPCTIO &tio, yield_t &yield, Duoram<Node>::Flat &A,
        Duoram<Node>::OblivIndex<RegXS,1> oidx, RegXS oidx_oldptrs, RegXS ptr,
        RegXS nodeptrs, RegBS p_bal_l, RegBS p_bal_r, RegBS &bal_upd, RegBS c_prime,
        RegXS cs_ptr, RegBS imb, RegBS &F_ri, avl_del_return &ret_struct);

    void updateRetStruct(MPCTIO &tio, yield_t &yield, RegXS ptr, RegBS F_rs,
        RegBS F_dh, RegBS F_ri, RegBS &bal_upd, avl_del_return &ret_struct);

    std::tuple<bool, RegBS> del(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS del_key,
        Duoram<Node>::Flat &A, RegBS F_af, RegBS F_fs, int TTL,
        avl_del_return &ret_struct);

    std::tuple<RegBS, RegBS, RegBS, RegBS> updateBalanceDel(MPCTIO &tio, yield_t &yield,
        RegBS bal_l, RegBS bal_r, RegBS bal_upd, RegBS child_dir);

    bool lookup(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS key,
        Duoram<Node>::Flat &A, int TTL, RegBS isDummy, Node *ret_node);

    void pretty_print(const std::vector<Node> &R, value_t node,
        const std::string &prefix, bool is_left_child, bool is_right_child);

    std::tuple<bool, bool, bool, address_t> check_avl(const std::vector<Node> &R,
        value_t node, value_t min_key, value_t max_key);

  public:
    AVL(int num_players, size_t size, bool opt_flag = true) :
            oram(num_players, size), OPTIMIZED(opt_flag) {
        this->MAX_SIZE = size;
        MAX_DEPTH = 0;
        while(size>0) {
          MAX_DEPTH+=1;
          size=size>>1;
        }
    };

    void init(){
        num_items=0;
        cur_max_index=0;
        empty_locations.clear();
    }

    void insert(MPCTIO &tio, yield_t &yield, const Node &node);

    // Deletes the first node that matches del_key. If an item with del_key
    // does not exist in the tree, it results in an explicit (non-oblivious)
    // failure.
    bool del(MPCTIO &tio, yield_t &yield, RegAS del_key);

    // Returns the first node that matches key
    bool lookup(MPCTIO &tio, yield_t &yield, RegAS key, Node *ret_node);

    // Non-obliviously initialize an AVL tree of a particular size
    void initialize(MPCTIO &tio, yield_t &yield, size_t depth);

    // Display and correctness check functions
    void pretty_print(MPCTIO &tio, yield_t &yield);
    bool check_avl(MPCTIO &tio, yield_t &yield);
    void print_oram(MPCTIO &tio, yield_t &yield);

    // For test functions ONLY:
    Duoram<Node>* get_oram() {
        return &oram;
    };

    RegXS get_root() {
        return root;
    };
};

void avl(MPCIO &mpcio, const PRACOptions &opts, char **args);
void avl_tests(MPCIO &mpcio, const PRACOptions &opts, char **args);

#endif
