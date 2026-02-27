# Contributing

## Development Setup

```bash
# Clone repository
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell

# Build Docker environment
docker-compose build dev

# Verify environment
docker-compose run --rm dev bash scripts/verify_environment.sh
```

## Running Tests

```bash
# All tests
docker-compose run --rm dev cargo test

# Specific test
docker-compose run --rm dev cargo test test_name

# With output
docker-compose run --rm dev cargo test -- --nocapture
```

## Code Style

```bash
# Format code
docker-compose run --rm dev cargo fmt

# Lint
docker-compose run --rm dev cargo clippy
```

## Project Structure

```
src/
├── anndata/      # AnnData (.h5ad) reader/writer
├── seurat/       # Seurat object handling
├── rds/          # RDS format reader/writer
├── ir/           # Intermediate representation
├── sparse/       # Sparse matrix operations
├── backend/      # Storage backend abstraction
├── validation/   # Validation engine
└── bin/          # CLI entry point
```

## Adding Features

1. Create feature branch
2. Implement changes
3. Add tests
4. Update documentation
5. Submit pull request

## Testing Guidelines

- Write unit tests for new functions
- Add integration tests for new features
- Ensure all tests pass before submitting

## Documentation

- Update relevant docs in `docs/`
- Keep both English and Chinese versions in sync
- Use clear, concise language

## Commit Messages

```
feat: add new feature
fix: fix bug
docs: update documentation
test: add tests
refactor: code refactoring
```

## Questions?

Open an issue on GitHub.
