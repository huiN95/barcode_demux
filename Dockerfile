# syntax=docker/dockerfile:1.6

# # ========= builder =========
# FROM 192.168.3.38:5000/algo/adacus_dev:1.11.0 AS builder
# WORKDIR /app

# # 如果仓库里没有 Cargo.lock，就只 copy toml，然后生成 lock（可选）
# COPY Cargo.toml ./
# RUN bash -lc 'set -eux; . "$HOME/.cargo/env"; cargo generate-lockfile'

# # 预编译依赖（缓存用）
# RUN --mount=type=cache,target=/root/.cargo/registry,sharing=locked \
#     --mount=type=cache,target=/root/.cargo/git,sharing=locked \
#     --mount=type=cache,target=/app/target,sharing=locked \
#     bash -lc 'set -eux; \
#     mkdir -p src && echo "fn main(){}" > src/main.rs; \
#     . "$HOME/.cargo/env"; \
#     cargo build --release'

# # 复制完整源码并构建
# COPY . .
# RUN --mount=type=cache,target=/root/.cargo/registry,sharing=locked \
#     --mount=type=cache,target=/root/.cargo/git,sharing=locked \
#     --mount=type=cache,target=/app/target,sharing=locked \
#     bash -lc 'set -eux; \
#     . "$HOME/.cargo/env"; \
#     make build-and-install; \
#     strip -s /usr/bin/primer_demux || true'


# # ========= runtime =========
# FROM debian:bookworm AS runtime

# ENV RUST_LOG=info \
#     RUST_BACKTRACE=1

# # apt 加速：BuildKit cache（确保 DOCKER_BUILDKIT=1）
# RUN sed -i 's|http://deb.debian.org/debian|http://mirrors.tuna.tsinghua.edu.cn/debian|g' /etc/apt/sources.list \
#     --mount=type=cache,target=/var/cache/apt,sharing=locked \
#     --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
#     set -eux; \
#     apt-get update; \
#     apt-get install -y --no-install-recommends ca-certificates; \
#     rm -rf /var/lib/apt/lists/*

# # root 运行：默认工作目录
# WORKDIR /work

# COPY --from=builder /usr/bin/primer_demux /usr/bin/primer_demux

# ENTRYPOINT ["/usr/bin/primer_demux"]


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

ENTRYPOINT ["/usr/bin/barcode_demux"]


