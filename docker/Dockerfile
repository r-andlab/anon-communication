FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y wget git build-essential net-tools iproute2 iperf iputils-ping libbsd-dev libboost-all-dev numactl time
WORKDIR /root
COPY . prac
WORKDIR /root/prac
RUN make -j `nproc`
