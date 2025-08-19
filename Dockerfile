FROM alpine:latest

# Install Python, pip, and bash
RUN apk add --no-cache bash python3 py3-pip

# Set working directory
WORKDIR /app

# Copy all files into container
COPY . .

# Create virtual environment & install dependencies
RUN python3 -m venv /opt/venv && \
    . /opt/venv/bin/activate && \
    pip install --no-cache-dir -r requirements.txt

# Ensure the venv's Python & pip are used by default
ENV PATH="/opt/venv/bin:$PATH"

# Default to bash shell
CMD ["bash"]
