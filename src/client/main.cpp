#include <iostream>
#include <rpc/client.h>
#include <rpc/rpc_error.h>
#include <rdpf.hpp>
#include <cell.hpp>
#include <duoram.hpp>

bool post_message(const char *ip, RegAS index, RegAS val, RegAS proof) {
  try {
    rpc::client client(ip, 8080);
    
    auto result = client.call("post", index.ashare, val.ashare, proof.ashare).as<int>();

    if (result != 0) {
      std::cout << "Protocol Error: " << result << std::endl;
      return false;
    }
  } catch (rpc::system_error &e) {
    std::cout << "System Error: " << e.what() << std::endl;
    return false;
  } catch (rpc::rpc_error &e) {
    std::cout << "RPC Error: " << e.what() << std::endl;
    return false;
  }

  return true;
}

struct TokenDB {
  std::vector<uint64_t> access_list;
  std::vector<uint64_t> proofs;
};

TokenDB get_tokenDB(const char *ip) {
  TokenDB db;

  rpc::client client(ip, 8080);
  
  db.access_list = client.call("get_access_list").as<std::vector<uint64_t>>();
  db.proofs = client.call("get_proofdb").as<std::vector<uint64_t>>();

  return db;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage: anoncomm_client <server1_ip> <server2_ip>" << std::endl;
    return -1;
  }

  const char *server1_ip = argv[1];
  const char *server2_ip = argv[2];

  TokenDB tdb = get_tokenDB(server1_ip);

  RegAS msg, proof;

  msg.randomize(8);
  proof.ashare = tdb.proofs[0];

  RegAS idx1;
  idx1.ashare = 0;

  if (!post_message(server1_ip, idx1, msg, proof)) return -1;

  RegAS idx2;
  idx2.ashare = 1;
  
  if (!post_message(server2_ip, idx2, msg, proof)) return -1;
  
  return 0;
}
