 
 

#ifndef LOWMC_LOWMC_H__
#define LOWMC_LOWMC_H__

#include <type_traits>          // std::conditional
#include <cstring>              // std::memset
#include <bitset>               // std::bitset
#include <array>                // std::array
#include <cstddef>              // std::size_t
#include <thread>               // std::thread
#include <future>               // std::async

#include "../block.h"
 #include "../mpcio.hpp"
#include <boost/utility.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/thread.hpp>


namespace lowmc
{

/// the default block length in bytes (must be 128 or 256)
constexpr size_t default_block_len = 128;

/// the default number of rounds
constexpr size_t default_rounds    = 29;

/// the default number of s-boxes per round (must be between 1 and block_len/3)
constexpr size_t default_sboxes    = 11;

/// default number of slices in the bitsliced implementation
 

/// mask used in the evaluation of s-boxes
static const std::string sbox_mask = "01001001001001001001001001001001001001001"
	"00100100100100100100100100100100100100100100100100100100100100100100100100"
	"10010010010010010010010010010010010010010010010010010010010010010010010010"
	"0100100100100100100100100100100100100100100100100100100100100100100";

template <size_t block_len = default_block_len,
          size_t rounds    = default_rounds,
          size_t sboxes    = default_sboxes>
class lowmc
{
  public:
 
	static_assert(block_len == 128 || block_len == 256,
		"block_len must be 128 or 256");
	static_assert(rounds > 0,
		"rounds must be positive");
	static_assert(sboxes > 0,
		"sboxes must be positive");
	static_assert(sboxes < block_len/3,
		"sboxes must be less than 3*block_len");
	static_assert(rounds > 0,
		"rounds must be positive");

	/// number of bits in a block (Contraints: @ref block_size = 128 or @ref block_size = 256)
	static constexpr auto block_size = block_len;

	/// number of rounds (Constraints: 0 < @ref num_rounds)
	static constexpr auto num_rounds = rounds;

	/// number of s-boxes per round (Contraints: 0 < @ref sboxes_per_round ≤ @ref block_size/3)
	static constexpr auto sboxes_per_round = sboxes;

	/// number of bits per block unaffected by the s-boxes in each round
	static constexpr auto identity_len = block_len - 3 * sboxes;

	/// total number of s-boxes across all rounds
	static constexpr auto sboxes_total = sboxes * rounds;

	/// type representing a single lowmc block
	using block_t = std::conditional_t<block_len == 128,
		block<__m128i>, block<__m256i>>;


  //private:
	/// mask for the highest-order bit in each s-box
	block_t mask = sbox_mask;
    block_t maska = mask;
 //    block_t maska = block_t(std::string(sbox_mask, 0, identity_len  - 1));
 	
	// /// mask for middle-order bit in each s-box
	   block_t maskb = maska >> 1;
	// /// mask for low-order bit in each s-box
	   block_t maskc = maska >> 2;
	// /// mask for the all-but-the-highest-order bit in each s-box
	   block_t maskbc = maskb | maskc;
    static const block_t* round_constants; // Declare as static
    static const uint64_t matrices[29][2*128];//      [rounds][(block_len/64) * block_len];
 
/** @defgroup group1 Regular encryption
 *  This is the first group
 *  @{
 */

	/// Plain-ol' ECB-mode encryption of a 1-block message
	inline auto encrypt(const block_t & msg) const
	{

		//std::cout << "Plain Old Encrypt!" << std::endl;
		auto c = msg;
		for (size_t i = 0; i < rounds; ++i)
		{
			c = substitute(c);
		 	auto mat = reinterpret_cast<const block_t *>(matrices[i]);
		    c = mul(mat, c, round_constants[i]);
		}
		return c;
	} 
 

	inline auto encrypt2_p0p1(block_t msg_share, MPCTIO& tio, yield_t &yield) const
	{

		//std::cout << "encrypt2_p0p1 called by P0 and P1" << std::endl << std::endl;
		//printf("This function encrypts the message\n");
 
		block_t c0 = msg_share;
		
		for (size_t i = 0; i < rounds; ++i)
		{
			
			//std::cout << "i = " << i << std::endl;
			block_t gamma, blind0, blinded_c1;
		
			block_t blind_recv;
			
			tio.recv_server(&blind_recv, sizeof(block_t));
		
			//std::cout << "blind_recv = " << blind_recv[0] << " <> " << blind_recv[1] << std::endl; 

			tio.recv_server(&gamma, sizeof(gamma));
	 
		    block_t blinded_c = c0 ^ blind_recv;

			tio.queue_peer(&blinded_c, sizeof(blinded_c));

			yield();

			block_t blinded_c_recv;

			tio.recv_peer(&blinded_c_recv, sizeof(blinded_c_recv));
			
			//std::cout << "gamma_recv = " << gamma.mX[0] << " <> " << gamma.mX[1] << std::endl; 
			c0 = substitute2_p0p1(c0, blind_recv, blinded_c_recv, gamma);

			auto mat = reinterpret_cast<const block_t *>(matrices[i]);

			bool p = false;
			
			if(tio.player() == 1) p = true;
			
		    c0 = mul(mat, c0, p ? round_constants[i] : 0);

			// std::cout << std::endl << " --------------------------------- \n\n";
 
		}
		return c0;
	} // encrypt2_p0p1_new

 	/// s-boxes
	inline auto & substitute(block_t & msg) const
	{
		auto srli1 = (msg >> 1) & maskbc;
		auto srli2 = (msg >> 2) & maskc;

		auto tmp = msg & srli1;

		auto bc = (tmp << 2) & maska;
		auto ac = (msg & srli2) << 1;
		auto ab = (tmp >> 1) & maskc;

		msg ^= (bc | ac | ab) ^ srli1 ^ srli2;

		return msg;
	} // substitute

	// inline auto substitute2_p2(randomness & rand0, randomness & rand1,
	// 		randomness & rand2) const0
	inline auto substitute2_p2(block_t blind0, block_t blind1, block_t blind2, MPCTIO& tio, yield_t& yield) const
	{
		 	
		auto tmp1 = ((blind0 >> 1) & blind1) ^ ((blind1 >> 1) & blind0);
		auto tmp2 = ((blind0 >> 2) & blind1) ^ ((blind1 >> 2) & blind0);

		auto bc = (tmp1 << 2) & maska;
		auto ac = (tmp2 << 1) & maskb;
		auto ab = (tmp1 >> 1) & maskc;

		auto gamma0 = (bc | ac | ab) ^ blind2;
		auto gamma1 = blind2;

		// printf("gamma0 = %llu <> %llu\n", gamma0.mX[0], gamma0.mX[1]);
		// printf("gamma1 = %llu <> %llu\n", gamma1.mX[0], gamma1.mX[1]);
	
		return std::make_pair(gamma0, gamma1);
	} // substitute2_p2

	// inline void encrypt2_p2(randomness & rand0,	randomness & rand1, randomness & rand2) const
	inline void encrypt2_p2(MPCTIO& tio, yield_t & yield) const
	{
		//std::cout << "encrypt2_p2 called by P2" << std::endl << std::endl;

		for (size_t i = 0; i < rounds; ++i)
		{

			block_t blind0, blind1, blind2;

			arc4random_buf(&blind0, sizeof(block_t));
			arc4random_buf(&blind1, sizeof(block_t));
			arc4random_buf(&blind2, sizeof(block_t));
 
			tio.queue_p0(&blind0, sizeof(block_t));
			tio.queue_p1(&blind1, sizeof(block_t));
			yield();
			
			//auto [gamma0, gamma1] = substitute2_p2(rand0, rand1, rand2);
			auto [gamma0, gamma1] = substitute2_p2(blind0, blind1, blind2, tio, yield);
	
			// std::cout << "gamma0 = " << gamma0.mX[0] << " " << gamma0.mX[1] << std::endl; 
			// std::cout << "gamma1 = " << gamma1.mX[0] << " " << gamma1.mX[1] << std::endl; 
			tio.queue_p0(&gamma0, sizeof(gamma0));
			tio.queue_p1(&gamma1, sizeof(gamma1));

			yield();

 
		}
	} // encrypt2_p2
 
 




	inline auto mul(const block_t * matrix, const block_t & msg, const block_t & constant) const
	{
		auto result = constant;		

		for (size_t i = 0; i < sizeof(block_t)/sizeof(uint64_t); ++i)
		{
			uint64_t bitset = static_cast<typename block_t::value_type>(msg)[i];
			while (bitset != 0)
			{
				result ^= matrix[i*64 + __builtin_ctzl(bitset)];
				bitset ^= bitset & -bitset;
			}
		}

		return result; 
	} // mul

  
	inline auto & substitute2_p0p1(block_t & share0, const block_t & blind0,
		const block_t & blinded_share1, const block_t & gamma0) const
	{
		auto blinded_msg = share0 ^ blinded_share1;

		auto srli1 = (share0 >> 1) & maskbc;
		auto srli2 = (share0 >> 2) & maskc;

		auto tmp = (blinded_msg & srli1) ^ (blind0 & (blinded_share1 >> 1));

		auto bc = (tmp << 2) & maska;
		auto ac = (((blinded_msg & srli2) ^ (blind0 & (blinded_share1 >> 2))) << 1) & maskb;
		auto ab = (tmp >> 1) & maskc;

		share0 ^= (bc | ac | ab) ^ srli1 ^ srli2 ^ gamma0;
		return share0;
	} // substitute2_p0p1
	
}; // class lowmc

} // namespace lowmc
 
#endif // LOWMC_LOWMC_H__


 