import logging
from concurrent import futures

import grpc
import bulletin_board_pb2 as pb2
import bulletin_board_pb2_grpc as pb2_grpc

# Audit server address for verification requests
AUDIT_ADDRESS = "localhost:50053"

class DatabaseServiceHandler(pb2_grpc.DatabaseServiceServicer):
    """Database server logic for storing shares (Server 0)."""
    def __init__(self, server_id):
        self.server_id = server_id
        self.pending_shares = {}  # Pending shares awaiting verification, {index: share_bytes}
        # Prepare a stub to communicate with the audit server
        self.audit_stub = pb2_grpc.AuditServiceStub(grpc.insecure_channel(AUDIT_ADDRESS))

    def SubmitShare(self, request, context):
        """RPC method for client to submit a share of a message."""
        idx = request.index
        share_bytes = request.share
        logging.info(f"{self.server_id}: Received share for index {idx}")
        # Store the share in pending buffer
        self.pending_shares[idx] = share_bytes
        # Send the share to the audit server for verification
        response = self.audit_stub.VerifyShare(pb2.VerifyRequest(index=idx, share=share_bytes, server_id=self.server_id))
        if response.status == "PENDING":
            # Audit is waiting for the other share
            logging.info(f"{self.server_id}: Audit status PENDING for index {idx} (awaiting other share).")
            return pb2.ShareResponse(status="PENDING")
        elif response.status == "VERIFIED":
            # Audit found the matching share and triggered commit
            logging.info(f"{self.server_id}: Audit VERIFIED index {idx}. Share stored successfully.")
            return pb2.ShareResponse(status="STORED")
        else:
            logging.error(f"{self.server_id}: Audit returned error for index {idx}.")
            return pb2.ShareResponse(status="ERROR")

    def CommitShare(self, request, context):
        """RPC method called by audit server to commit (store) a verified share."""
        idx = request.index
        share_bytes = self.pending_shares.get(idx)
        if share_bytes is None:
            logging.warning(f"{self.server_id}: No pending share for index {idx} to commit.")
            return pb2.CommitResponse(status="NOT_FOUND")
        # Write the share to local storage (append to a text file for persistence)
        filename = f"{self.server_id}_storage.txt"
        try:
            share_text = share_bytes.decode('utf-8')
        except UnicodeDecodeError:
            # If not valid text, store as hex string for readability
            share_text = share_bytes.hex()
        with open(filename, "a") as f:
            f.write(f"{idx}: {share_text}\n")
        logging.info(f"{self.server_id}: Committed share for index {idx} to {filename}.")
        # Remove from pending buffer since it's now stored
        del self.pending_shares[idx]
        return pb2.CommitResponse(status="OK")

def serve():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=5))
    # Initialize handler with unique server_id (e.g., "server0")
    pb2_grpc.add_DatabaseServiceServicer_to_server(DatabaseServiceHandler("server0"), server)
    server.add_insecure_port("[::]:50051")  # Listen on port 50051 for client
    logging.info("Starting Database Server 0 on port 50051...")
    server.start()
    server.wait_for_termination()

if __name__ == "__main__":
    serve()
