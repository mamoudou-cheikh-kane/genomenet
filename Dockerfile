FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY app.py ./
COPY assets ./assets
COPY data ./data
# 7860 : port attendu par Hugging Face Spaces ; Render fournit $PORT
EXPOSE 7860
CMD ["sh", "-c", "gunicorn app:server --bind 0.0.0.0:${PORT:-7860} --workers 2 --timeout 120"]
