FROM rust:1-bookworm

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3 \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src

RUN cargo build --release --locked

COPY validation/run_particles_scaling.py /app/validation/run_particles_scaling.py

WORKDIR /data
ENTRYPOINT ["/app/target/release/carbon_nanotube"]
CMD ["run", "--sp", "/data"]
