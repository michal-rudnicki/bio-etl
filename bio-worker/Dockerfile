FROM python:3.11

WORKDIR /app

COPY worker.py .

# Instalujemy wymagane biblioteki
RUN pip install biopython psycopg2-binary

CMD ["python", "worker.py"]