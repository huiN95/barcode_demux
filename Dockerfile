# syntax=docker/dockerfile:1.6
# ========= builder =========
FROM 192.168.3.38:5000/algo/adacus_dev:1.11.0 AS builder

WORKDIR /app
COPY . .

RUN bash -lc 'set -eux; \
    . "$HOME/.cargo/env"; \
    make build-and-install; \
    strip -s /usr/bin/adapter_demux || true'

# 3) 再复制真实源码（这一步变更最频繁，放后面）
# COPY . .


# ========= runtime =========
# 方案A：更小的 Debian slim（推荐先试这个）
FROM debian:bookworm-slim AS runtime

ENV RUST_LOG=info \
    RUST_BACKTRACE=1

# 常用运行时基础：证书 + 时区数据（可选）+ 最小依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 非 root 运行：同时创建一个可写的工作目录
# RUN useradd -m -u 10001 -s /usr/sbin/nologin appuser \
#     && mkdir -p /work \
#     && chown -R appuser:appuser /work

# 可选：如果你程序会写日志/输出到当前目录，强烈建议把 WORKDIR 指到可写目录
WORKDIR /work

COPY --from=builder /usr/bin/barcode_demux /usr/bin/barcode_demux

USER root
ENTRYPOINT ["/usr/bin/barcode_demux"]