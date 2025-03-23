FROM python:3.11

WORKDIR /app

COPY . .

# Instalujemy wymagane biblioteki
RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "-m", "main"]