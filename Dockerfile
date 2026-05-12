# syntax=docker/dockerfile:1.6

# ========= builder =========
FROM 192.168.3.38:5000/algo/adacus_dev:1.11.0 AS builder

WORKDIR /app

# 1) 先拷贝依赖描述文件，利用缓存
COPY Cargo.toml Cargo.lock ./

# 若是 workspace，可按需补充
# COPY crates/*/Cargo.toml crates/*/

# 预热依赖缓存
RUN mkdir -p src && echo 'fn main(){}' > src/main.rs && \
    . "$HOME/.cargo/env" && \
    cargo build --release && \
    rm -rf target/release/deps && \
    rm -rf src

# 2) 再拷贝真实源码
COPY . .

# 3) 真正构建与安装
RUN . "$HOME/.cargo/env" && \
    cargo build --release && \
    cp target/release/barcode_demux /usr/bin/

# ========= runtime =========
FROM debian:bookworm-slim AS runtime

ENV RUST_LOG=info \
    RUST_BACKTRACE=1

RUN sed -i 's|http://deb.debian.org/debian|http://mirrors.tuna.tsinghua.edu.cn/debian|g' /etc/apt/sources.list \
    --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root

COPY --from=builder /usr/bin/barcode_demux /usr/bin/barcode_demux

# ENTRYPOINT ["/usr/bin/barcode_demux"]


