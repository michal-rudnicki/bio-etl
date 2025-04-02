# === Ścieżki ===
APP_DIR=bio_worker
ENV_FILE=bio_worker/.env

# === Docker ===
up:
	d	test -f $(ENV_FILE)|| (echo "Brak pliku .env! Utwórz go najpierw."; exit 1)
	docker-compose --env-file $(ENV_FILE) up --build

down:
	test -f $(ENV_FILE) || (echo "Brak pliku .env! Utwórz go najpierw."; exit 1)
	docker-compose --env-file $(ENV_FILE) down -v

restart:
	test -f $(ENV_FILE) || (echo "Brak pliku .env! Utwórz go najpierw."; exit 1)
	docker-compose --env-file $(ENV_FILE) down -v
	docker-compose --env-file $(ENV_FILE) up --build

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
	set -a && source $(ENV_FILE) && set +a && pytest bio_worker/bio_worker/tests

# === Wszystko razem ===
check: format lint test

# === Czyszczenie (opcjonalnie) ===
clean:
	docker system prune -af

pre-commit:
	pre-commit run --all-files