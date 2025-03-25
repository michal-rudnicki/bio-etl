# === Ścieżki ===
APP_DIR=bio_worker

# === Docker ===
up:
	docker-compose up --build

down:
	docker-compose down

restart:
	docker-compose down -v
	docker-compose --env-file .env up --build

logs:
	docker-compose logs -f

# === Formatowanie ===
format:
	black $(APP_DIR)
	isort $(APP_DIR)

# === Lintowanie ===
lint:
	flake8 $(APP_DIR)
	pylint $(APP_DIR)

# === Testy lokalne ===
test:
	set -a && source .env && set +a && pytest bio_worker/bio_worker/tests

# === Wszystko razem ===
check: format lint test

# === Czyszczenie (opcjonalnie) ===
clean:
	docker system prune -af