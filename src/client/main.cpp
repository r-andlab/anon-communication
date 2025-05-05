#include <iostream>
#include <rpc/client.h>
#include <rpc/rpc_error.h>
#include <rdpf.hpp>
#include <cell.hpp>
#include <duoram.hpp>

bool post_message(const char *ip, RegAS ind, RegAS val, RegAS proof) {
  try {
    auto result = client.call("post", ind, val, proof).as<int>();

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

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage: anoncomm_client <server1_ip> <server2_ip>" << std::endl;
    return -1;
  }

  const char *server1_ip = argv[1];
  const char *server2_ip = argv[2];

  Writer writer(1001);

  RegAS ind, val, proof;
  inserted_val.randomize(8);

  if (!post_message(server1_ip, ind, val, proof)) return -1;
  if (!post_message(server2_ip, ind, val, proof)) return -1;
  
  return 0;
}
