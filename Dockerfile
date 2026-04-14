# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm AS builder

WORKDIR /app

COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release --all-features --bin genomemask --locked && \
    strip target/release/genomemask

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/genomemask /usr/local/bin/genomemask

# Set up non-root user
RUN useradd -m -u 1000 cuser && \
    chmod +x /usr/local/bin/genomemask

USER cuser
WORKDIR /data

RUN genomemask --help
