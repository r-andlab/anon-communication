# Anonymous Bulletin Board with Secret Sharing

This project implements a secure, anonymous bulletin board system using secret sharing techniques to protect message confidentiality. The system distributes each message across multiple servers, ensuring that no single server can read the complete content.

## System Architecture

The system consists of three server components and a client:

1. **Database Server 0** (server0.py)
2. **Database Server 1** (server1.py)
3. **Audit Server** (audit_server.py)
4. **Client** (client.py)

### Core Concept: Secret Sharing

Messages are split into two "shares" using additive secret sharing. Each share individually appears random, but when combined they reconstruct the original message. This system uses a simple modular addition scheme where:

- For each byte in the message, a random byte is generated for Share 0
- Share 1 is calculated so that `(Share 0 + Share 1) mod 256 = Original Byte`
- To reconstruct: `(Share 0 + Share 1) mod 256` recovers the original message

## Component Details

### Client (client.py)

The client is responsible for:
- Taking a message to be posted on the bulletin board
- Splitting it into two shares using the secret sharing algorithm
- Sending each share to a different database server
- Optionally verifying reconstruction for demonstration purposes

Usage example:
```python
# Message "Hello, anonymous world!" is split and sent to servers
```

### Database Servers (server0.py, server1.py)

These two servers function identically but run on different ports (50051 and 50052):
- Receive shares from clients via gRPC
- Store pending shares temporarily
- Forward shares to the audit server for verification
- Commit verified shares to local storage files (server0_storage.txt, server1_storage.txt)
- Implement two RPC methods:
  - `SubmitShare`: Called by clients to submit message shares
  - `CommitShare`: Called by the audit server to commit verified shares

### Audit Server (audit_server.py)

The audit server is the central coordinator:
- Receives shares from both database servers
- Pairs matching shares by index
- Verifies both shares are from the same original message
- Instructs database servers to commit their shares when verified
- Optionally reconstructs and logs the complete message for demonstration/debugging
- Implements the RPC method:
  - `VerifyShare`: Called by database servers to verify share integrity

### Protocol Buffer Files

- **bulletin_board.proto**: Defines the protocol buffer messages and service interfaces
- **bulletin_board_pb2.py**: Generated Python code for protocol messages
- **bulletin_board_pb2_grpc.py**: Generated Python code for gRPC service stubs and servers

## Storage Files

- **server0_storage.txt**: Persistent storage for shares on Server 0
- **server1_storage.txt**: Persistent storage for shares on Server 1

Each line follows the format: `<index>: <share_data>`

## Security Properties

This system provides:
- **Confidentiality**: No single server can read complete messages
- **Integrity**: The audit server ensures shares come from the same message
- **Availability**: Messages can be reconstructed as long as both shares exist
- **Auditability**: The audit server logs message reconstructions

## Running the System

1. Start the audit server:
   ```
   python audit_server.py
   ```

2. Start both database servers:
   ```
   python server0.py
   python server1.py
   ```

3. Run the client to post a message:
   ```
   python client.py
   ```

## Limitations

- This implementation uses insecure gRPC channels for simplicity
- In a production system, secure channels and authentication would be required
- The system requires all three servers to be honest (no collusion)
- Basic error handling with minimal recovery mechanisms

## Extensions

Possible extensions to this system:
- Support for more advanced secret sharing schemes (Shamir's Secret Sharing)
- Threshold-based reconstruction (requiring only k of n shares)
- Message retrieval API for authorized users
- Digital signatures for message authentication