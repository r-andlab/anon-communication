import logging
from concurrent import futures

import grpc
import bulletin_board_pb2 as pb2
import bulletin_board_pb2_grpc as pb2_grpc

# Addresses of the two database servers (for callback)
SERVER0_ADDRESS = "localhost:50051"
SERVER1_ADDRESS = "localhost:50052"

class AuditServiceHandler(pb2_grpc.AuditServiceServicer):
    """Audit server logic: verify matching shares from two database servers."""
    def __init__(self):
        # Store one share until the matching share arrives: { index: {server_id: share_bytes} }
        self.pending_shares = {}
        # Prepare stubs to call back to the two database servers for commit
        self.stub0 = pb2_grpc.DatabaseServiceStub(grpc.insecure_channel(SERVER0_ADDRESS))
        self.stub1 = pb2_grpc.DatabaseServiceStub(grpc.insecure_channel(SERVER1_ADDRESS))

    def VerifyShare(self, request, context):
        """RPC method called by a database server to verify a received share."""
        idx = request.index
        share_bytes = request.share
        server_id = request.server_id
        logging.info(f"Audit: Received share for index {idx} from {server_id}")

        # If this is the first share for the given index, store it and wait for the other share
        if idx not in self.pending_shares:
            self.pending_shares[idx] = {server_id: share_bytes}
            logging.info(f"Audit: Stored share from {server_id} for index {idx}, waiting for matching share.")
            return pb2.VerifyResponse(status="PENDING")
        else:
            # Another share for this index is already pending. Retrieve it.
            other_entries = self.pending_shares[idx]
            if server_id in other_entries:
                # This scenario shouldn’t normally happen (duplicate from same server)
                logging.warning(f"Audit: Duplicate share from {server_id} for index {idx} ignored.")
                return pb2.VerifyResponse(status="ERROR")
            # Combine this share with the pending share from the other server
            (other_server_id, other_share) = next(iter(other_entries.items()))
            other_entries[server_id] = share_bytes  # add the new share
            logging.info(f"Audit: Both shares for index {idx} received (from {other_server_id} and {server_id}).")

            # (Optional) Reconstruct the original message by adding byte values mod 256
            combined_bytes = bytes((a + b) % 256 for a, b in zip(other_share, share_bytes))
            try:
                # Try decoding as UTF-8 text for logging
                msg_text = combined_bytes.decode('utf-8')
            except UnicodeDecodeError:
                msg_text = combined_bytes  # if not text, keep as bytes
            logging.info(f"Audit: Reconstructed message for index {idx}: {msg_text}")

            # Instruct both database servers to commit their shares for this index
            logging.info(f"Audit: Instructing servers to commit index {idx}.")
            self.stub0.CommitShare(pb2.CommitRequest(index=idx))
            self.stub1.CommitShare(pb2.CommitRequest(index=idx))
            # Remove the pending entry since verification is done
            del self.pending_shares[idx]
            return pb2.VerifyResponse(status="VERIFIED")

def serve():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    pb2_grpc.add_AuditServiceServicer_to_server(AuditServiceHandler(), server)
    server.add_insecure_port("[::]:50053")  # Listen on port 50053
    logging.info("Starting Audit server on port 50053...")
    server.start()
    server.wait_for_termination()

if __name__ == "__main__":
    serve()
