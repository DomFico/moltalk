FROM python:3.12-slim
# RDKit's drawing module links against X11 render libraries that the slim image omits.
RUN apt-get update && apt-get install -y --no-install-recommends libxrender1 libxext6 libexpat1 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY pyproject.toml ./
COPY moltalk ./moltalk
RUN pip install --no-cache-dir . && useradd --create-home chem
USER chem
# 0.0.0.0 is only the in-container bind address. Publish the port on 127.0.0.1 (see README);
# Host/Origin validation still only accepts localhost unless MOLTALK_ALLOWED_HOSTS is set.
ENV HOST=0.0.0.0 PORT=8000 OPENBLAS_NUM_THREADS=1
EXPOSE 8000
CMD ["moltalk", "--transport", "streamable-http"]
