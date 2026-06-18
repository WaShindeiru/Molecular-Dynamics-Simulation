# Molecular-Dynamics-Simulation

Run with info logs, from defined particles, parameters.json and particles_initial.json are in one directory:
```bash
RUST_LOG=info cargo run --profile release-checks -- run --sp <path>
```

Run with info logs, with particles generator, parameters.json and generator_config.json are in one directory:
```bash
RUST_LOG=info cargo run --profile release-checks -- run --sg <path>
```

Run with info logs, with particles generator, specify parameters.json and generator_config.json separately:
```bash
RUST_LOG=info cargo run --profile release-checks -- run -s <simulation-configuration-path> -g <generator-configuration-path>
```

Run with info logs, from defined partilces, specify parameters.json and particles_initial.json separately:
```bash
RUST_LOG=info cargo run --profile release-checks -- run -s <simulation-configuration-path> -p <particles_initial-path>
```
