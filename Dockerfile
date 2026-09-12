FROM rust:1-bookworm

WORKDIR /app

COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src

RUN cargo build --release --locked

WORKDIR /data
ENTRYPOINT ["/app/target/release/carbon_nanotube"]
CMD ["run", "--sp", "/data"]
