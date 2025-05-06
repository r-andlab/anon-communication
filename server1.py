import logging
from concurrent import futures

import grpc
import bulletin_board_pb2 as pb2
import bulletin_board_pb2_grpc as pb2_grpc

AUDIT_ADDRESS = "localhost:50053"

class DatabaseServiceHandler(pb2_grpc.DatabaseServiceServicer):
    """Database server logic for storing shares (Server 1)."""
    def __init__(self, server_id):
        self.server_id = server_id
        self.pending_shares = {}
        self.audit_stub = pb2_grpc.AuditServiceStub(grpc.insecure_channel(AUDIT_ADDRESS))

    def SubmitShare(self, request, context):
        idx = request.index
        share_bytes = request.share
        logging.info(f"{self.server_id}: Received share for index {idx}")
        self.pending_shares[idx] = share_bytes
        response = self.audit_stub.VerifyShare(pb2.VerifyRequest(index=idx, share=share_bytes, server_id=self.server_id))
        if response.status == "PENDING":
            logging.info(f"{self.server_id}: Audit status PENDING for index {idx} (awaiting other share).")
            return pb2.ShareResponse(status="PENDING")
        elif response.status == "VERIFIED":
            logging.info(f"{self.server_id}: Audit VERIFIED index {idx}. Share stored successfully.")
            return pb2.ShareResponse(status="STORED")
        else:
            logging.error(f"{self.server_id}: Audit returned error for index {idx}.")
            return pb2.ShareResponse(status="ERROR")

    def CommitShare(self, request, context):
        idx = request.index
        share_bytes = self.pending_shares.get(idx)
        if share_bytes is None:
            logging.warning(f"{self.server_id}: No pending share for index {idx} to commit.")
            return pb2.CommitResponse(status="NOT_FOUND")
        filename = f"{self.server_id}_storage.txt"
        try:
            share_text = share_bytes.decode('utf-8')
        except UnicodeDecodeError:
            share_text = share_bytes.hex()
        with open(filename, "a") as f:
            f.write(f"{idx}: {share_text}\n")
        logging.info(f"{self.server_id}: Committed share for index {idx} to {filename}.")
        del self.pending_shares[idx]
        return pb2.CommitResponse(status="OK")

def serve():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=5))
    pb2_grpc.add_DatabaseServiceServicer_to_server(DatabaseServiceHandler("server1"), server)
    server.add_insecure_port("[::]:50052")  # Listen on port 50052
    logging.info("Starting Database Server 1 on port 50052...")
    server.start()
    server.wait_for_termination()

if __name__ == "__main__":
    serve()
