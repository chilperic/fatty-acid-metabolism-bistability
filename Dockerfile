FROM python:3.10-slim
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir -e .[dev]
CMD ["bash", "-lc", "pytest -q && python scripts/01_simulate_timecourse.py && python scripts/02_scan_acetyl_influx.py"]
