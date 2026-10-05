FROM python:3.12-slim
# RDKit's drawing module needs X11 render libraries the slim image omits; OPSIN (IUPAC numbering) needs Java.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 libxext6 libexpat1 default-jre-headless curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*
# OPSIN 2.9.0, verified against the SHA-256 published with the GitHub release.
ARG OPSIN=opsin-cli-2.9.0-jar-with-dependencies.jar
ARG OPSIN_SHA256=c2e29326c281f87b59a05d934d8589adac6e9d17b95b984931b3e739111b360f
RUN mkdir -p /opt/opsin \
    && curl -fsSL -o /opt/opsin/$OPSIN https://github.com/dan2097/opsin/releases/download/2.9.0/$OPSIN \
    && echo "$OPSIN_SHA256  /opt/opsin/$OPSIN" | sha256sum -c -
WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY moltalk ./moltalk
RUN pip install --no-cache-dir . && useradd --create-home chem
USER chem
# 0.0.0.0 is only the in-container bind address. Host/Origin validation still only accepts localhost unless
# MOLTALK_ALLOWED_HOSTS lists the public hostname. Cloud Run sets PORT.
ENV HOST=0.0.0.0 PORT=8080 OPENBLAS_NUM_THREADS=1 \
    MOLTALK_OPSIN_JAR=/opt/opsin/opsin-cli-2.9.0-jar-with-dependencies.jar \
    MOLTALK_WORKERS=2 MOLTALK_TIMEOUT_S=25
EXPOSE 8080
CMD ["moltalk", "--transport", "streamable-http"]
