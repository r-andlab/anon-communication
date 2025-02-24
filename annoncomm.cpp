#include <functional>
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"


#include "annoncomm.hpp"   


void SAM::post_message(MPCTIO &tio, yield_t & yield, RegAS ind, RegAS val, RegAS proof) {
    
    auto RemiseDB = oram.flat(tio, yield);
    auto ProofDB = proofdb.flat(tio, yield);
 
    CDPF cdpf = tio.cdpf(yield);

    RegAS zero_share;

    zero_share.ashare = 0;
    
    RegBS audit_share, final_audit_check; 

    //uint64_t ind_rec = mpc_reconstruct(tio, yield, ind);

    //std::cout << "ind_rec = " << ind_rec << std::endl;
    
    RegAS ProofShare = ProofDB[ind];

    RegBS lt_audit, eq_audit, gt_audit, lt_zero_msg, eq_zero_msg, gt_zero_msg;

    run_coroutines(tio,
        [&cdpf, &tio, &yield, &proof, &ProofShare, &lt_audit, &eq_audit, &gt_audit](yield_t &inner_yield) {
            std::tie(lt_audit, eq_audit, gt_audit) = cdpf.compare(tio, inner_yield, proof - ProofShare, tio.aes_ops());
        },
        [&cdpf, &tio, &yield, &val, &zero_share, &lt_zero_msg, &eq_zero_msg, &gt_zero_msg](yield_t &inner_yield) {
            std::tie(lt_zero_msg, eq_zero_msg, gt_zero_msg) = cdpf.compare(tio, inner_yield, val - zero_share, tio.aes_ops());
        }
    );

    RegBS tmp1 =  eq_zero_msg ^ eq_audit;    
    RegBS tmp2; 
   
    mpc_and(tio, yield, tmp2, eq_zero_msg, eq_audit);
 
    final_audit_check = tmp1 ^ tmp2;
 
    bool final_audit_check_reconstruction = mpc_reconstruct(tio, yield, final_audit_check);

    #ifdef VERBOSE
        std::cout << "final_audit_check_reconstruction = " << (int) final_audit_check_reconstruction << std::endl;
    #endif
   
    if(final_audit_check_reconstruction) 
    {
        RemiseDB[ind] = val;
        
        #ifdef VERBOSE
        std::cout << "WRITE REQUEST IS SUCCESS" << std::endl;
        #endif

    }
    else
    {
        RemiseDB[ind] = zero_share;
        
        #ifdef VERBOSE
        std::cout << "write request failed!" << std::endl;
        #endif
    }
}



RegAS SAM::read_message(MPCTIO &tio, yield_t & yield, RegAS ind, RegAS proof_submitted) {
      
    auto TokenDB = tokendb.flat(tio, yield);

    auto RemiseDB = oram.flat(tio, yield);
    
    RegAS proof_share = TokenDB[ind];

    CDPF cdpf = tio.cdpf(yield);

    auto[lt_read_auth, eq_read_auth, gt_read_auth] = cdpf.compare(tio, yield, proof_share - proof_submitted, tio.aes_ops());   
    
    bool read_auth = mpc_reconstruct(tio, yield, eq_read_auth);

    RegAS read_out;

    if(read_auth) 
    {
        read_out = RemiseDB[ind];
    }
    else
    {
        #ifdef VERBOSE
         //  std::cout << "Read Failed" << std::endl;
        #endif
    }

    return read_out;
}
 

 void SAM::provide_or_remove_access(MPCTIO &tio, yield_t & yield, RegAS ind, RegAS proof_submitted) {

  //  std::cout << "Provide Access" << std::endl;
    
    CDPF cdpf = tio.cdpf(yield);
    
    auto ProofDB = proofdb.flat(tio, yield);
    auto TokenDB = tokendb.flat(tio, yield);

    RegAS ProofShare = ProofDB[ind];

    auto[lt_audit, eq_audit, gt_audit] = cdpf.compare(tio, yield, proof_submitted - ProofShare, tio.aes_ops()); // Checks if the message submitted is a zero message

    bool auth_reconstruction = mpc_reconstruct(tio, yield, eq_audit);

    if(auth_reconstruction) {
        TokenDB[ind] = proof_submitted;
    } 
    else
    {
     // std::cout << "read access update failed" << std::endl;  
    }
}
 


 
void SAM::initialize_proofDB(MPCTIO &tio, yield_t & yield, size_t dbsize, Writer writer) {
   
    auto ProofDB = proofdb.flat(tio, yield);
    size_t count = 0;

    for(size_t j = 0; j < dbsize; ++j)
    {
       uint64_t proof;
       arc4random_buf(&proof, sizeof(uint64_t));
       RegAS proof_j;
      
       if(tio.player() == 0) proof_j.ashare = 0;
       if(tio.player() == 1) proof_j.ashare = proof;

       ProofDB[j] = proof_j;
        
       if(j < 1000) 
       {
          writer.access_list[count] = j;
          writer.proofs[count] = proof;
          ++count;
       }
       
    }
}


void SAM::publish_bulletin_board(MPCTIO &tio, yield_t & yield, size_t dbsize) {

auto RemiseDB = oram.flat(tio, yield);

auto published = RemiseDB.reconstruct();

    for(size_t j = 0; j < dbsize; ++j)
    {
        std::cout << published[j].share() << std::endl; 
    }
}
 
 


void AnnonComm(MPCIO & mpcio,  const PRACOptions & opts, char ** args) {

    std::cout << "Remise code is being run" << std::endl;

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    int nargs = 0;

    while (args[nargs] != nullptr) {
        ++nargs;
    }

    int dbsize = 0;
 
    size_t n_correct_posts = 10;
 
    size_t n_failed_posts = 10;

    size_t n_cover_traffic = 10;

    size_t n_provide_access = 10;
    

    size_t n_reads = 10;
    
    for (int i = 0; i < nargs; i += 2) {
        std::string option = args[i];
        if (option == "-m" && i + 1 < nargs) {
            dbsize = std::atoi(args[i + 1]);
        } else if (option == "-v" && i + 1 < nargs) {
            n_correct_posts = std::atoi(args[i + 1]);
        } else if (option == "-f" && i + 1 < nargs) {
            n_failed_posts = std::atoi(args[i + 1]);
        } else if (option == "-t" && i + 1 < nargs) {
            n_cover_traffic = std::atoi(args[i + 1]);
        } else if (option == "-g" && i + 1 < nargs) {
            n_provide_access = std::atoi(args[i + 1]);
        } else if (option == "-r" && i + 1 < nargs) {
            n_reads = std::atoi(args[i + 1]);
        }
    }

    run_coroutines(tio, [ & tio, dbsize, n_correct_posts, n_failed_posts, n_cover_traffic, n_provide_access, n_reads, &mpcio](yield_t & yield) {
 
        size_t size = size_t(1) << dbsize;

        SAM bulletinboard(tio.player(), size);
 
        Writer writer(1001);

        
        bulletinboard.initialize_proofDB(tio, yield, size, writer);

        // std::cout << "\n===== Init Stats =====\n";
        // tio.sync_lamport();
        // mpcio.dump_stats(std::cout);
        
        mpcio.reset_stats();
        tio.reset_lamport();
 
        RegAS ind, inserted_val, proof_submitted;
        inserted_val.randomize(8);
        
        //uint64_t index_to_which_posted = mpc_reconstruct(tio, yield, ind);
        //std::cout << "index_to_which_posted = " << index_to_which_posted << std::endl;
        
        for(size_t j = 0; j < n_correct_posts; ++ j)
        {
            if(tio.player() == 0) proof_submitted.ashare = 0;       
            if(tio.player() == 1) proof_submitted.ashare = writer.proofs[j];
            
            if(tio.player() == 0) ind.ashare = 0;       
            if(tio.player() == 1) ind.ashare = j;

            bulletinboard.post_message(tio, yield, ind,  inserted_val, proof_submitted);
            
            #ifdef VERBOSE
            std::cout << std::endl << std::endl << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            #endif
        }

        std::cout << "\n===== Correct Post Stats =====\n";
        std::cout << "n_correct_posts = " << n_correct_posts << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();
 
        
        #ifdef VERBOSE
            bulletinboard.publish_bulletin_board(tio, yield, size);
        #endif

        for(size_t j = 0; j < n_failed_posts; ++ j)
        {
            ind.randomize();
            bulletinboard.post_message(tio, yield, ind,  inserted_val, proof_submitted);

            #ifdef VERBOSE
            std::cout << std::endl << std::endl << std::endl;
            std::cout << std::endl << std::endl << std::endl;
            #endif
        }

        std::cout << "\n\n\n\n\n\n===== Failed Post Stats =====\n";
        std::cout << "n_failed_posts = " << n_failed_posts << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();
    
    
        for(size_t j = 0; j < n_cover_traffic; ++j)
        {
         ind.randomize();
         inserted_val.ashare = 0;
         bulletinboard.post_message(tio, yield, ind,  inserted_val, proof_submitted);
        }

        std::cout << "\n\n\n\n\n\n===== Cover Traffic Stats =====\n";
        std::cout << "n_cover_traffic = " << n_cover_traffic << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();



        for(size_t j = 0; j < n_reads; ++j)
        {
         bulletinboard.read_message(tio, yield, ind, proof_submitted);
        }

        std::cout << "\n\n\n\n\n\n===== Reading Stats =====\n";
        std::cout << "n_reads = " << n_reads << std::endl;
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();
    }
    );
}
