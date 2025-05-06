import random
import grpc
import bulletin_board_pb2 as pb2
import bulletin_board_pb2_grpc as pb2_grpc

def secret_share_message(message: str):
    """Split a message into two additive shares such that share0 + share1 = message (byte-wise mod 256)."""
    data = message.encode('utf-8')
    # Generate a random byte for each byte of the message
    share0_bytes = bytes(random.randrange(0, 256) for _ in range(len(data)))
    # Compute share1 as the difference (mod 256) between message bytes and share0 bytes
    share1_bytes = bytes((m - s) % 256 for m, s in zip(data, share0_bytes))
    return share0_bytes, share1_bytes

if __name__ == "__main__":
    # Example message to post
    message = "Hello, anonymous world!"
    index = 1  # Bulletin board index for this message

    # Create two additive shares of the message
    share0, share1 = secret_share_message(message)
    # Set up gRPC channels to the two database servers
    channel0 = grpc.insecure_channel("localhost:50051")
    channel1 = grpc.insecure_channel("localhost:50052")
    stub0 = pb2_grpc.DatabaseServiceStub(channel0)
    stub1 = pb2_grpc.DatabaseServiceStub(channel1)

    print(f"Client: Submitting message '{message}' at index {index}")
    # Submit each share to its respective server
    response0 = stub0.SubmitShare(pb2.ShareRequest(index=index, share=share0))
    response1 = stub1.SubmitShare(pb2.ShareRequest(index=index, share=share1))
    # Print out the responses from both servers
    print(f"Client: Server0 response -> {response0.status}")
    print(f"Client: Server1 response -> {response1.status}")

    # Verify that combining the shares reproduces the original message (for demo purposes)
    reconstructed = bytes((a + b) % 256 for a, b in zip(share0, share1))
    try:
        recon_msg = reconstructed.decode('utf-8')
    except UnicodeDecodeError:
        recon_msg = reconstructed
    print(f"Client: Reconstructed message from shares -> {recon_msg}")
    print("Client: Submission complete. You can check server0_storage.txt and server1_storage.txt for stored shares.")
